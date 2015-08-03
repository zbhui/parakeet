
#include "NavierStokesProblem.h"

using Eigen::Vector3d;

template<>
InputParameters validParams<NavierStokesProblem>()
{
  InputParameters params = validParams<CFDProblem>();

  return params;
}

NavierStokesProblem::NavierStokesProblem(const InputParameters &params) :
	CFDProblem(params)
{
}



