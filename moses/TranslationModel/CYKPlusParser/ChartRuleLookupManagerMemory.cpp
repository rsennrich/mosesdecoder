/***********************************************************************
  Moses - factored phrase-based language decoder
  Copyright (C) 2011 University of Edinburgh

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 ***********************************************************************/

#include <iostream>
#include "ChartRuleLookupManagerMemory.h"
#include "DotChartInMemory.h"
#include "CompletedRuleCollection.h"

#include "moses/ChartParser.h"
#include "moses/InputType.h"
#include "moses/ChartParserCallback.h"
#include "moses/StaticData.h"
#include "moses/NonTerminal.h"
#include "moses/ChartCellCollection.h"
#include "moses/TranslationModel/PhraseDictionaryMemory.h"

using namespace std;

namespace Moses
{

ChartRuleLookupManagerMemory::ChartRuleLookupManagerMemory(
  const ChartParser &parser,
  const ChartCellCollectionBase &cellColl,
  const PhraseDictionaryMemory &ruleTable)
  : ChartRuleLookupManagerCYKPlus(parser, cellColl)
  , m_ruleTable(ruleTable)
  , m_softMatchingMap(StaticData::Instance().GetSoftMatches())
{
  UTIL_THROW_IF2(m_dottedRuleColls.size() != 0,
		  "Dotted rule collection not correctly initialized");

  size_t sourceSize = parser.GetSize();
  m_dottedRuleColls.resize(sourceSize);

  const PhraseDictionaryNodeMemory &rootNode = m_ruleTable.GetRootNode();

  m_isSoftMatching = !m_softMatchingMap.empty();

  for (size_t ind = 0; ind < m_dottedRuleColls.size(); ++ind) {
#ifdef USE_BOOST_POOL
    DottedRuleInMemory *initDottedRule = m_dottedRulePool.malloc();
    new (initDottedRule) DottedRuleInMemory(rootNode);
#else
    DottedRuleInMemory *initDottedRule = new DottedRuleInMemory(rootNode);
#endif

    DottedRuleColl *dottedRuleColl = new DottedRuleColl(sourceSize - ind + 1);
    dottedRuleColl->Add(0, initDottedRule); // init rule. stores the top node in tree

    m_dottedRuleColls[ind] = dottedRuleColl;
  }
}

ChartRuleLookupManagerMemory::~ChartRuleLookupManagerMemory()
{
  RemoveAllInColl(m_dottedRuleColls);
}

void ChartRuleLookupManagerMemory::GetChartRuleCollection(
  const WordsRange &range,
  size_t lastPos, //unused
  ChartParserCallback &outColl)
{
  size_t relEndPos = range.GetEndPos() - range.GetStartPos();
  size_t absEndPos = range.GetEndPos();

  // not necessary for normal search, but implements rule pruning for incremental search
  CompletedRuleCollection tmpColl = CompletedRuleCollection();

  // MAIN LOOP. create list of nodes of target phrases

  // get list of all rules that apply to spans at same starting position
  DottedRuleColl &dottedRuleCol = *m_dottedRuleColls[range.GetStartPos()];
  const DottedRuleList &expandableDottedRuleList = dottedRuleCol.GetExpandableDottedRuleList();

  const ChartCellLabel &sourceWordLabel = GetSourceAt(absEndPos);

  // loop through the rules
  // (note that expandableDottedRuleList can be expanded as the loop runs
  //  through calls to ExtendPartialRuleApplication())
  for (size_t ind = 0; ind < expandableDottedRuleList.size(); ++ind) {
    // rule we are about to extend
    const DottedRuleInMemory &prevDottedRule = *expandableDottedRuleList[ind];
    // we will now try to extend it, starting after where it ended
    size_t startPos = prevDottedRule.IsRoot()
                      ? range.GetStartPos()
                      : prevDottedRule.GetWordsRange().GetEndPos() + 1;

    // search for terminal symbol
    // (if only one more word position needs to be covered)
    if (startPos == absEndPos) {

      // look up in rule dictionary, if the current rule can be extended
      // with the source word in the last position
      const Word &sourceWord = sourceWordLabel.GetLabel();
      const PhraseDictionaryNodeMemory::TerminalMap & terminals = prevDottedRule.GetLastNode().GetTerminalMap();

      // if node has small number of terminal edges, test word equality for each.
      if (terminals.size() < 5) {
        for (PhraseDictionaryNodeMemory::TerminalMap::const_iterator iter = terminals.begin(); iter != terminals.end(); ++iter) {
          const Word & word = iter->first;
          if (TerminalEqualityPred()(word, sourceWord)) {
            const PhraseDictionaryNodeMemory *node = & iter->second;
#ifdef USE_BOOST_POOL
            DottedRuleInMemory *dottedRule = m_dottedRulePool.malloc();
            new (dottedRule) DottedRuleInMemory(*node, sourceWordLabel,
                                            prevDottedRule);
#else
            DottedRuleInMemory *dottedRule = new DottedRuleInMemory(*node,
              sourceWordLabel,
              prevDottedRule);
#endif
            dottedRuleCol.Add(relEndPos+1, dottedRule);
          }
        }
      }
      // else, do hash lookup
      else {
        const PhraseDictionaryNodeMemory *node = prevDottedRule.GetLastNode().GetChild(sourceWord);

        // if we found a new rule -> create it and add it to the list
        if (node != NULL) {
          // create the rule
#ifdef USE_BOOST_POOL
          DottedRuleInMemory *dottedRule = m_dottedRulePool.malloc();
          new (dottedRule) DottedRuleInMemory(*node, sourceWordLabel,
                                            prevDottedRule);
#else
          DottedRuleInMemory *dottedRule = new DottedRuleInMemory(*node,
              sourceWordLabel,
              prevDottedRule);
#endif
          dottedRuleCol.Add(relEndPos+1, dottedRule);
        }
      }
    }

    // search for non-terminals
    size_t endPos, stackInd;

    // span is already complete covered? nothing can be done
    if (startPos > absEndPos)
      continue;

    else if (startPos == range.GetStartPos() && range.GetEndPos() > range.GetStartPos()) {
      // We're at the root of the prefix tree so won't try to cover the full
      // span (i.e. we don't allow non-lexical unary rules).  However, we need
      // to match non-unary rules that begin with a non-terminal child, so we
      // do that in two steps: during this iteration we search for non-terminals
      // that cover all but the last source word in the span (there won't
      // already be running nodes for these because that would have required a
      // non-lexical unary rule match for an earlier span).  Any matches will
      // result in running nodes being appended to the list and on subsequent
      // iterations (for this same span), we'll extend them to cover the final
      // word.
      endPos = absEndPos - 1;
      stackInd = relEndPos;
    } else {
      endPos = absEndPos;
      stackInd = relEndPos + 1;
    }


    ExtendPartialRuleApplication(prevDottedRule, startPos, endPos, stackInd,
                                 dottedRuleCol);
  }

  // list of rules that that cover the entire span
  DottedRuleList &rules = dottedRuleCol.Get(relEndPos + 1);

  // look up target sides for the rules
  DottedRuleList::const_iterator iterRule;
  for (iterRule = rules.begin(); iterRule != rules.end(); ++iterRule) {
    const DottedRuleInMemory &dottedRule = **iterRule;
    const PhraseDictionaryNodeMemory &node = dottedRule.GetLastNode();

    // look up target sides
    const TargetPhraseCollection &tpc = node.GetTargetPhraseCollection();

    // add the fully expanded rule (with lexical target side)
    AddCompletedRule(dottedRule, tpc, tmpColl, outColl);
  }

  for (vector<CompletedRule*>::const_iterator iter = tmpColl.begin(); iter != tmpColl.end(); ++iter) {
    outColl.Add((*iter)->GetTPC(), (*iter)->GetStackVector(), range);
  }

  tmpColl.Clear();
  dottedRuleCol.Clear(relEndPos+1);
}

// Given a partial rule application ending at startPos-1 and given the sets of
// source and target non-terminals covering the span [startPos, endPos],
// determines the full or partial rule applications that can be produced through
// extending the current rule application by a single non-terminal.
void ChartRuleLookupManagerMemory::ExtendPartialRuleApplication(
  const DottedRuleInMemory &prevDottedRule,
  size_t startPos,
  size_t endPos,
  size_t stackInd,
  DottedRuleColl & dottedRuleColl)
{

  // target non-terminal labels for the span
  const ChartCellLabelSet &targetNonTerms = GetTargetLabelSet(startPos, endPos);

    if (targetNonTerms.GetSize() == 0) {
      return;
    }

#if !defined(UNLABELLED_SOURCE)
    // source non-terminal labels for the span
    const InputPath &inputPath = GetParser().GetInputPath(startPos, endPos);
    const std::vector<bool> &sourceNonTermArray = inputPath.GetNonTerminalArray();

    // can this ever be true? Moses seems to pad the non-terminal set of the input with [X]
    if (inputPath.GetNonTerminalSet().size() == 0) {
      return;
    }
#endif

    // note where it was found in the prefix tree of the rule dictionary
    const PhraseDictionaryNodeMemory &node = prevDottedRule.GetLastNode();

    // non-terminal labels in phrase dictionary node
    const PhraseDictionaryNodeMemory::NonTerminalMap & nonTermMap = node.GetNonTerminalMap();


    if (nonTermMap.empty()) {
      return;
    }

    // loop over possible expansions of the rule
    PhraseDictionaryNodeMemory::NonTerminalMap::const_iterator p;
    PhraseDictionaryNodeMemory::NonTerminalMap::const_iterator end = nonTermMap.end();
    for (p = nonTermMap.begin(); p != end; ++p) {
      // does it match possible source and target non-terminals?
#if defined(UNLABELLED_SOURCE)
      const Word &targetNonTerm = p->first;
#else
      const PhraseDictionaryNodeMemory::NonTerminalMapKey &key = p->first;
      const Word &sourceNonTerm = key.first;
      // check if source label matches
      if (! sourceNonTermArray[sourceNonTerm[0]->GetId()]) {
        continue;
      }
      const Word &targetNonTerm = key.second;
#endif
      //soft matching of NTs
      if (m_isSoftMatching && !m_softMatchingMap[targetNonTerm[0]->GetId()].empty()) {
        const std::vector<Word>& softMatches = m_softMatchingMap[targetNonTerm[0]->GetId()];
        for (std::vector<Word>::const_iterator softMatch = softMatches.begin(); softMatch != softMatches.end(); ++softMatch) {
          const ChartCellLabel *cellLabel = targetNonTerms.Find(*softMatch);
          if (cellLabel == NULL) {
            continue;
          }
          // create new rule
          const PhraseDictionaryNodeMemory &child = p->second;
#ifdef USE_BOOST_POOL
          DottedRuleInMemory *rule = m_dottedRulePool.malloc();
          new (rule) DottedRuleInMemory(child, *cellLabel, prevDottedRule);
#else
          DottedRuleInMemory *rule = new DottedRuleInMemory(child, *cellLabel,
              prevDottedRule);
#endif
          dottedRuleColl.Add(stackInd, rule);
        }
      } // end of soft matches lookup

      const ChartCellLabel *cellLabel = targetNonTerms.Find(targetNonTerm);
      if (cellLabel == NULL) {
        continue;
      }
      // create new rule
      const PhraseDictionaryNodeMemory &child = p->second;
#ifdef USE_BOOST_POOL
      DottedRuleInMemory *rule = m_dottedRulePool.malloc();
      new (rule) DottedRuleInMemory(child, *cellLabel, prevDottedRule);
#else
      DottedRuleInMemory *rule = new DottedRuleInMemory(child, *cellLabel,
          prevDottedRule);
#endif
      dottedRuleColl.Add(stackInd, rule);
    }
}

}  // namespace Moses
