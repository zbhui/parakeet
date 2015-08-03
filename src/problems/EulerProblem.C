
#include "EulerProblem.h"

using Eigen::Vector3d;

template<>
InputParameters validParams<EulerProblem>()
{
  InputParameters params = validParams<CFDProblem>();

  return params;
}

EulerProblem::EulerProblem(const InputParameters &params) :
	CFDProblem(params)
{
}

