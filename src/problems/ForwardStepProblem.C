
#include "ForwardStepProblem.h"
#include "boost/assign.hpp"
using namespace boost::assign;

template<>
InputParameters validParams<ForwardStepProblem>()
{
  InputParameters params = validParams<Riemann2DProblem>();

  return params;
}

ForwardStepProblem::ForwardStepProblem(const InputParameters &params) :
	Riemann2DProblem(params)
{
	_mach = 3;
	Real mach2 = _mach*_mach;
	Real r1 =1.4;
	Real p1 = 1;
	Real c = sqrt(_gamma*p1/r1);

	_shock_depart[0] += r1,3, 0, 0, p1;
}

int ForwardStepProblem::pointLocator(const Point& p)
{
	return 0;
}

