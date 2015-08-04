
#include "EulerProblem.h"

using Eigen::Vector3d;

template<>
InputParameters validParams<EulerProblem>()
{
  InputParameters params = validParams<CFDProblem>();
  params.set<MooseEnum>("vis_type") = "INVISCOUS";
  return params;
}

EulerProblem::EulerProblem(const InputParameters &params) :
	CFDProblem(params)
{
}

