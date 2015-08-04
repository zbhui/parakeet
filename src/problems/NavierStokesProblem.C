
#include "NavierStokesProblem.h"

using Eigen::Vector3d;

template<>
InputParameters validParams<NavierStokesProblem>()
{
  InputParameters params = validParams<CFDProblem>();
  params.set<MooseEnum>("vis_type") = "CONSTANT";
  return params;
}

NavierStokesProblem::NavierStokesProblem(const InputParameters &params) :
	CFDProblem(params)
{
}

