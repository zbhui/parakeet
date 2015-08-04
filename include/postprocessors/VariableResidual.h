
#pragma once
#include "GeneralPostprocessor.h"

class VariableResidual : public GeneralPostprocessor
{
public:
	VariableResidual(const InputParameters &parameters);

  virtual void initialize() {}
  virtual void execute() {}

  virtual Real getValue();
};

template<>
InputParameters validParams<VariableResidual>();
