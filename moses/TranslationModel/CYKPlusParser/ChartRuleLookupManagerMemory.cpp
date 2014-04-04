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

#include "moses/ChartParser.h"
#include "moses/InputType.h"
#include "moses/ChartParserCallback.h"
#include "moses/StaticData.h"
#include "moses/NonTerminal.h"
#include "moses/ChartCellCollection.h"
#include "moses/FactorCollection.h"
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

  size_t sourceSize = parser.GetSize();

  m_completedRules.resize(sourceSize);

  m_isSoftMatching = !m_softMatchingMap.empty();
}

void ChartRuleLookupManagerMemory::GetChartRuleCollection(
  const WordsRange &range,
  size_t lastPos,
  ChartParserCallback &outColl)
{
  size_t startPos = range.GetStartPos();
  size_t absEndPos = range.GetEndPos();

  m_lastPos = lastPos;
  m_stackVec.clear();
  m_outColl = &outColl;
  m_unaryPos = absEndPos-1; // rules ending in this position are unary and should not be added to collection

  CreateFastLookupVectors(startPos, lastPos);

  const PhraseDictionaryNodeMemory &rootNode = m_ruleTable.GetRootNode();

  // size-1 terminal rules
  if (startPos == absEndPos) {
    const Word &sourceWord = GetSourceAt(absEndPos).GetLabel();
    const PhraseDictionaryNodeMemory *child = rootNode.GetChild(sourceWord);

    // if we found a new rule -> directly add it to the out collection
    if (child != NULL) {
        const TargetPhraseCollection &tpc = child->GetTargetPhraseCollection();
        outColl.Add(tpc, m_stackVec, range);
    }
  }
  // all rules starting with nonterminal
  else if (absEndPos > startPos) {
    GetNonTerminalExtensionFixedSpan(&rootNode, startPos, absEndPos-1);
    // all (non-unary) rules starting with terminal
    if (absEndPos == startPos+1) {
      GetTerminalExtension(&rootNode, absEndPos-1);
    }
  }

  // copy temporarily stored rules to out collection
  CompletedRuleCollection rules = m_completedRules[absEndPos];
  for (vector<CompletedRule*>::const_iterator iter = rules.begin(); iter != rules.end(); ++iter) {
    outColl.Add((*iter)->GetTPC(), (*iter)->GetStackVector(), range);
  }

  m_completedRules[absEndPos].Clear();

}

// Create vector 
void ChartRuleLookupManagerMemory::CreateFastLookupVectors(size_t firstPos,
    size_t lastPos) {

    size_t numNonTerms = FactorCollection::Instance().GetNumNonTerminals();

    //re-use data structure from last step, but remove chart cells that would break max-chart-span
    m_fastLookupVector.resize(lastPos+1);
    for (size_t startPos = firstPos+1; startPos <= lastPos; startPos++) {
        ChartCellMatrix & cellMatrix = m_fastLookupVector[startPos];
        cellMatrix.resize(numNonTerms);
        for (size_t i = 0; i < numNonTerms; i++) {
            ChartCellVector &matches = cellMatrix[i];
            for (size_t j = 0; j < matches.size(); j++) {
                if (matches[j].first > lastPos) {
                    matches.resize(j);
                    break;
                }
            }
        }
    }

    // populate vector with all chart cells that start at current start position (may still have some redundancy)
    ChartCellMatrix & cellMatrix = m_fastLookupVector[firstPos];
    cellMatrix.clear();
    cellMatrix.resize(numNonTerms);
    for (size_t endPos = firstPos; endPos <= lastPos; endPos++) {

        // target non-terminal labels for the span
        const ChartCellLabelSet &targetNonTerms = GetTargetLabelSet(firstPos, endPos);

        if (targetNonTerms.GetSize() == 0) {
            continue;
        }

#if !defined(UNLABELLED_SOURCE)
        // source non-terminal labels for the span
        const InputPath &inputPath = GetParser().GetInputPath(firstPos, endPos);
        const std::vector<bool> &sourceNonTermArray = inputPath.GetNonTerminalArray();

        // can this ever be true? Moses seems to pad the non-terminal set of the input with [X]
        if (inputPath.GetNonTerminalSet().size() == 0) {
            continue;
        }
#endif

        for (size_t i = 0; i < numNonTerms; i++) {
            const ChartCellLabel *cellLabel = targetNonTerms.Find(i);
            if (cellLabel != NULL) {
                cellMatrix[i].push_back(ChartCellPos(endPos, cellLabel));
            }
        }
    }

}

// if a (partial) rule matches, add it to list completed rules (if non-unary and non-empty), and try find expansions that have this partial rule as prefix.
void ChartRuleLookupManagerMemory::AddAndExtend(
    const PhraseDictionaryNodeMemory *node,
    size_t endPos) {

    const TargetPhraseCollection &tpc = node->GetTargetPhraseCollection();
    // add target phrase collection (except if rule is empty or unary)
    if (!tpc.IsEmpty() && endPos != m_unaryPos) {
      m_completedRules[endPos].Add(tpc, m_stackVec, *m_outColl);
    }

    // get all further extensions of rule (until reaching end of sentence or max-chart-span)
    if (endPos < m_lastPos) {
      if (!node->GetTerminalMap().empty()) {
        GetTerminalExtension(node, endPos+1);
      }
      if (!node->GetNonTerminalMap().empty()) {
          GetNonTerminalExtension(node, endPos+1);
      }
    }
}


// search all possible terminal extensions of a partial rule (pointed at by node) at a given position
// recursively try to expand partial rules into full rules up to m_lastPos.
void ChartRuleLookupManagerMemory::GetTerminalExtension(
    const PhraseDictionaryNodeMemory *node,
    size_t pos) {

    const Word &sourceWord = GetSourceAt(pos).GetLabel();
    const PhraseDictionaryNodeMemory::TerminalMap & terminals = node->GetTerminalMap();

    // if node has small number of terminal edges, test word equality for each.
    if (terminals.size() < 5) {
      for (PhraseDictionaryNodeMemory::TerminalMap::const_iterator iter = terminals.begin(); iter != terminals.end(); ++iter) {
        const Word & word = iter->first;
        if (word == sourceWord) {
          const PhraseDictionaryNodeMemory *child = & iter->second;
          AddAndExtend(child, pos);
          break;
        }
      }
    }
    // else, do hash lookup
    else {
      const PhraseDictionaryNodeMemory *child = node->GetChild(sourceWord);
      if (child != NULL) {
        AddAndExtend(child, pos);
      }
    }
}

// search all nonterminal possible nonterminal extensions of a partial rule (pointed at by node) for a given span (StartPos, endPos).
// recursively try to expand partial rules into full rules up to m_lastPos.
void ChartRuleLookupManagerMemory::GetNonTerminalExtensionFixedSpan(
    const PhraseDictionaryNodeMemory *node,
    size_t startPos,
    size_t endPos) {

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

    // non-terminal labels in phrase dictionary node
    const PhraseDictionaryNodeMemory::NonTerminalMap & nonTermMap = node->GetNonTerminalMap();

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
          m_stackVec.push_back(cellLabel);
          AddAndExtend(&child, endPos);
          m_stackVec.pop_back();
        }
      } // end of soft matches lookup

      const ChartCellLabel *cellLabel = targetNonTerms.Find(targetNonTerm);
      if (cellLabel == NULL) {
        continue;
      }
      // create new rule
      const PhraseDictionaryNodeMemory &child = p->second;
      m_stackVec.push_back(cellLabel);
      AddAndExtend(&child, endPos);
      m_stackVec.pop_back();
    }
}

// search all nonterminal possible nonterminal extensions of a partial rule (pointed at by node) for a variable span (starting from startPos).
// recursively try to expand partial rules into full rules up to m_lastPos.
void ChartRuleLookupManagerMemory::GetNonTerminalExtension(
    const PhraseDictionaryNodeMemory *node,
    size_t startPos) {

    const ChartCellMatrix &fastLookup = m_fastLookupVector.at(startPos);

    // non-terminal labels in phrase dictionary node
    const PhraseDictionaryNodeMemory::NonTerminalMap & nonTermMap = node->GetNonTerminalMap();

    // loop over possible expansions of the rule
    PhraseDictionaryNodeMemory::NonTerminalMap::const_iterator p;
    PhraseDictionaryNodeMemory::NonTerminalMap::const_iterator end = nonTermMap.end();
    for (p = nonTermMap.begin(); p != end; ++p) {
      // does it match possible source and target non-terminals?
#if defined(UNLABELLED_SOURCE)
      const Word &targetNonTerm = p->first;
#else
      const Word &targetNonTerm = p->first.second;
#endif
      //soft matching of NTs
      if (m_isSoftMatching && !m_softMatchingMap[targetNonTerm[0]->GetId()].empty()) {
        const std::vector<Word>& softMatches = m_softMatchingMap[targetNonTerm[0]->GetId()];
        const PhraseDictionaryNodeMemory &child = p->second;
        for (std::vector<Word>::const_iterator softMatch = softMatches.begin(); softMatch != softMatches.end(); ++softMatch) {
          const ChartCellVector &matches = fastLookup.at((*softMatch)[0]->GetId());
          for (ChartCellVector::const_iterator match = matches.begin(); match != matches.end(); ++match) {
            m_stackVec.push_back(match->second);
            AddAndExtend(&child, match->first);
            m_stackVec.pop_back();
          }
        }
      } // end of soft matches lookup

      const PhraseDictionaryNodeMemory &child = p->second;
      const ChartCellVector &matches = fastLookup[targetNonTerm[0]->GetId()];
      for (ChartCellVector::const_iterator match = matches.begin(); match != matches.end(); ++match) {
        m_stackVec.push_back(match->second);
        AddAndExtend(&child, match->first);
        m_stackVec.pop_back();
      }
    }
}

}  // namespace Moses
