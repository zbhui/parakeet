
#pragma once

#include "GeneralPostprocessor.h"

class NumTimeStep: public GeneralPostprocessor
{
public:
	NumTimeStep(const InputParameters &parameters);

  virtual void initialize() {}
  virtual void execute() {}

  virtual Real getValue();

protected:
  FEProblem & _feproblem;
};

template<>
InputParameters validParams<NumTimeStep>();
