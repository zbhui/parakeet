
#include "CFDInitialCondition.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDInitialCondition>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<unsigned int>("component", "分量");
  return params;
}

CFDInitialCondition::CFDInitialCondition(const InputParameters & parameters) :
    InitialCondition(parameters),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
    _component(getParam<unsigned int>("component"))
{}

Real CFDInitialCondition::value(const Point & p)
{
	  return value(_component, p);
}

Real CFDInitialCondition::value(int component, const Point & p)
{
	return _cfd_problem.initialCondition(p, component);
}
