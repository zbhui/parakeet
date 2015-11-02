#include "ProblemElementalL1Error.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<ProblemElementalL1Error>()
{
	InputParameters param = validParams<ElementIntegralPostprocessor>();
//	param.addParam<Real>("test", 0, "");
//	param.addRequiredParam<bool>("something", "ddd");
//	MooseEnum error_type(MooseEnum("L1 L2 LINF", "L2"));
//	params.addParam<MooseEnum>("error_type", error_type, "Error Type");
	return param;
}

ProblemElementalL1Error::ProblemElementalL1Error(const InputParameters & parameters) :
		ElementIntegralPostprocessor(parameters),
//		_test(getParam<Real>("test")),
//		_is_something(getParam<bool>("something")),
//		_max(0.0),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
		_nl(_cfd_problem.getNonlinearSystem()),
		_tid(parameters.get<THREAD_ID>("_tid")),
		_variables(_nl.getVariableNames()),
		_n_equations(_variables.size())

//		_error_type(getParam<MooseEnum>("error_type"))
	{
		for (int eq = 0; eq < _nl.getVariableNames().size(); ++eq)
		{
			MooseVariable &val = _cfd_problem.getVariable(_tid, _variables[eq]);
			_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
			addMooseVariableDependency(&val);
		}
	}

Real ProblemElementalL1Error::getValue()
{

  return ElementIntegralPostprocessor::getValue();
}

Real ProblemElementalL1Error::computeQpIntegral()
{
	Real err =  (*_uh[0])[_qp] - _cfd_problem.valueExact(_t, _q_point[_qp], 0);
	return std::abs(err);
}
