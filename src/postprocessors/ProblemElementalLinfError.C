/*
 * ProblemElementalLinfError.C
 *
 *	    无穷范数的求解
 *	    L_inf = max(|a1|,|a2|,|a3|,...)
 *
 *  Created on: 2015年10月29日
 *      Author: hui
 */

#include "ProblemElementalLinfError.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<ProblemElementalLinfError>()
{
	InputParameters param = validParams<ElementIntegralPostprocessor>();
	param.addParam<Real>("h",1.0,"value of the mesh cell");
	return param;
}

ProblemElementalLinfError::ProblemElementalLinfError(const InputParameters & parameters) :
		ElementIntegralPostprocessor(parameters),
		_max(0),
		_h(getParam<Real>("h")),
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

void
ProblemElementalLinfError::initialize()
{
  _max = 0;
}

Real ProblemElementalLinfError::getValue()
{
	gatherMax(_max);
	return _max;
}

void ProblemElementalLinfError::threadJoin(const UserObject & y)
{
  const ProblemElementalLinfError & pps = static_cast<const ProblemElementalLinfError &>(y);
  _max = std::max(_max, pps._max);
}

void ProblemElementalLinfError::execute()
{
  Real max(0);
  max = computeIntegral();
  _max = std::max(_max, max);
}

Real ProblemElementalLinfError::computeIntegral()
{
  Real sum = 0;

  for (_qp=0; _qp<_qrule->n_points(); _qp++)
    sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();
  return sum/(_h * _h);
}

Real ProblemElementalLinfError::computeQpIntegral()
{
	Real error(0.0);
	error = (*_uh[0])[_qp] - _cfd_problem.valueExact(_t, _q_point[_qp], 0);
	return std::abs(error);
}






