#include "ProblemElementalL2Error.h"

#include "CFDProblem.h"

template<>
InputParameters validParams<ProblemElementalL2Error>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  return params;
}

ProblemElementalL2Error::ProblemElementalL2Error(const InputParameters &parameters) :
		ElementIntegralPostprocessor(parameters),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
	_nl(_cfd_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
	_variables(_nl.getVariableNames()),
	_n_equations(_variables.size())
{
	for (int eq = 0; eq < _nl.getVariableNames().size(); ++eq)
	{
		MooseVariable &val = _cfd_problem.getVariable(_tid, _variables[eq]);
		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		addMooseVariableDependency(&val);
	}
}

Real ProblemElementalL2Error::getValue()
{
  return std::sqrt(ElementIntegralPostprocessor::getValue());
}

Real ProblemElementalL2Error::computeQpIntegral()
{
	Real err =  (*_uh[0])[_qp] - _cfd_problem.valueExact(_t, _q_point[_qp], 0);
	return err*err;
}
