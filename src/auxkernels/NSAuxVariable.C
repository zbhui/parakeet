
#include "NSAuxVariable.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<NSAuxVariable>()
{
  InputParameters params = validParams<MultiAuxKernel>();
  params.addRequiredCoupledVar("variables", "守恒变量");

  return params;
}

NSAuxVariable::NSAuxVariable(const InputParameters & parameters) :
	MultiAuxKernel(parameters),
    _fe_problem(*parameters.get<FEProblem *>("_fe_problem")),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
	_cfd_data(_cfd_problem)
{
	size_t n_equation = coupledComponents("variables");
	for (size_t i = 0; i < n_equation; ++i)
	{
		_uh.push_back(&coupledValue("variables", i));
	}

}

Real NSAuxVariable::computeValue()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
	}

	_cfd_data.reinit();

	if(_ivar == 0) return _cfd_data.p;
	else return _cfd_data.m;

}

