#include "EulerCellKernel.h"

template<>
InputParameters validParams<EulerCellKernel>()
{
  InputParameters params = validParams<MultiKernel>();
  params += validParams<CFDBase>();

  return params;
}

EulerCellKernel::EulerCellKernel(const InputParameters & parameters):
		MultiKernel(parameters),
		CFDBase(parameters),
		_cfd_data(_mach, _reynolds, _gamma, _prandtl)
{

}

void EulerCellKernel::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
		_cfd_data.uh[i] = (*_uh[i])[_qp];

	_cfd_data.reinit();
	fluxTerm();
}

Real EulerCellKernel::computeQpResidual(unsigned int p)
{
	return -_flux[p]*_grad_test[_i][_qp];
}

void EulerCellKernel::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
		_flux_old[q] =  _flux[q];

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _ds;
		_cfd_data.reinit();
		fluxTerm();
		for (int p = 0; p < _n_equation; ++p)
			_jacobi_variable[p][q] = (_flux[p] - _flux_old[p])/_ds;

		_cfd_data.uh[q] -= _ds;
	}
//	_cfd_data.reinit();
}

void EulerCellKernel::fluxTerm()
{
	for (int q = 0; q < _n_equation; ++q)
		_flux[q] =  _cfd_data.invis_flux[q];
}

Real EulerCellKernel::computeQpJacobian(unsigned int p, unsigned int q)
{
	return -(_jacobi_variable[p][q]*_phi[_j][_qp])*_grad_test[_i][_qp];
}




