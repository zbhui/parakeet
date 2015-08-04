
#pragma once
#include "GeneralPostprocessor.h"

class CFDResidual : public GeneralPostprocessor
{
public:
	CFDResidual(const InputParameters &parameters);

  virtual void initialize() {}
  virtual void execute() {}

  virtual Real getValue();
};

template<>
InputParameters validParams<CFDResidual>();
