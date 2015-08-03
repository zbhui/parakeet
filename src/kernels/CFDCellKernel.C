#include "CFDCellKernel.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDCellKernel>()
{
  InputParameters params = validParams<MultiKernel>();
  params.addParam<Real>("perturbation", 1E-08, "有限差分求Jacobian矩阵的变量增量");

  return params;
}

CFDCellKernel::CFDCellKernel(const InputParameters & parameters):
		MultiKernel(parameters),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
		_cfd_data(_cfd_problem),
		_perturbation(getParam<Real>("perturbation"))

{

}

void CFDCellKernel::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data.duh[i] = (*_grad_uh[i])[_qp];
	}
	_cfd_data.reinit();
	fluxTerm();
}

Real CFDCellKernel::computeQpResidual(unsigned int p)
{
	return -_flux[p]*_grad_test[_i][_qp];
}

void CFDCellKernel::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
		_flux_old[q] =  _flux[q];

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _perturbation;
		_cfd_data.reinit();
		fluxTerm();
		for (int p = 0; p < _n_equation; ++p)
			_flux_jacobi_variable[p][q] = (_flux[p] - _flux_old[p])/_perturbation;

		_cfd_data.uh[q] -= _perturbation;
	}
//	_cfd_data.reinit();
}

void CFDCellKernel::fluxTerm()
{
	for (int q = 0; q < _n_equation; ++q)
		_flux[q] =  _cfd_data.invis_flux[q];
}

Real CFDCellKernel::computeQpJacobian(unsigned int p, unsigned int q)
{
	return -(_flux_jacobi_variable[p][q]*_phi[_j][_qp])*_grad_test[_i][_qp];
}




