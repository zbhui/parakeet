
#include "CFDResidual.h"

#include "FEProblem.h"
#include "SubProblem.h"

template<>
InputParameters validParams<CFDResidual>()
{
  InputParameters params = validParams<GeneralPostprocessor>();
  return params;
}

CFDResidual::CFDResidual(const InputParameters &parameters) :
    GeneralPostprocessor(parameters)
{}

Real
CFDResidual::getValue()
{
  return _fe_problem.getNonlinearSystem()._initial_residual_before_preset_bcs;

//  return _fe_problem.getNonlinearSystem().residualVector(Moose::KT_ALL).l2_norm();
}

