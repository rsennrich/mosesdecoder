#pragma once
#include <string>
#include "StatelessFeatureFunction.h"

namespace Moses
{

class SetSourcePhrase : public StatelessFeatureFunction
{
public:
  SetSourcePhrase(const std::string &line);

  virtual bool IsUseable(const FactorMask &mask) const
  { return true; }

  virtual void Evaluate(const Phrase &source
						, const TargetPhrase &targetPhrase
						, ScoreComponentCollection &scoreBreakdown
						, ScoreComponentCollection &estimatedFutureScore) const;

  virtual void Evaluate(const InputType &input
                         , const InputPath &inputPath
                         , const TargetPhrase &targetPhrase
                         , ScoreComponentCollection &scoreBreakdown
                         , ScoreComponentCollection *estimatedFutureScore = NULL) const
  {}

  virtual void Evaluate(const Hypothesis& hypo,
                        ScoreComponentCollection* accumulator) const
  {}

  virtual void EvaluateChart(const ChartHypothesis &hypo,
                             ScoreComponentCollection* accumulator) const
  {}

  std::vector<float> DefaultWeights() const
  { return std::vector<float>(); }

};

}

