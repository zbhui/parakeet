#include "CFDCellKernel.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDCellKernel>()
{
  InputParameters params = validParams<MultiKernel>();
  params.addParam<Real>("perturbation", 1E-08, "有限差分求Jacobian矩阵的变量增量");
  params.addParam<IndicatorName>("shock_indicator", "The name of the Indicator that this Marker uses.");
  return params;
}

CFDCellKernel::CFDCellKernel(const InputParameters & parameters):
	MultiKernel(parameters),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
	_cfd_data(_cfd_problem),
	_perturbation(getParam<Real>("perturbation")),
	_has_artificial_vis(isParamValid("shock_indicator"))
{
	if(_has_artificial_vis)
		_error_vector = &(getErrorVector(getParam<IndicatorName>("shock_indicator")));
}

void CFDCellKernel::reinit()
{
	_cfd_data.reinit();
	reinitArtificialViscous();
	fluxTerm();
}

void CFDCellKernel::reinitArtificialViscous()
{
	Real h = _current_elem->hmax();
	if(_has_artificial_vis)
		for (int q = 0; q < _n_equation; ++q)
		{
			_artificial_viscous[q] = h*(*_error_vector)[_current_elem->id()]*_cfd_data.duh[q];
		}

	else
		for (int q = 0; q < _n_equation; ++q)
		{
			_artificial_viscous[q] = 0*_cfd_data.duh[q];
		}
}

void CFDCellKernel::reinitViscous()
{
	_cfd_data.reinitViscous();
	reinitArtificialViscous();
	for (int q = 0; q < _n_equation; ++q)
		_viscous[q] =  -_cfd_data.vis_flux[q] - _artificial_viscous[q];
}

void CFDCellKernel::fluxTerm()
{
	for (int q = 0; q < _n_equation; ++q)
	{
		_viscous[q] =  -_cfd_data.vis_flux[q] - _artificial_viscous[q];
		_flux[q] =  _cfd_data.invis_flux[q] + _viscous[q];
	}
}

void CFDCellKernel::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data.duh[i] = (*_grad_uh[i])[_qp];
	}
	reinit();
}

Real CFDCellKernel::computeQpResidual(unsigned int p)
{
	return -_flux[p]*_grad_test[_i][_qp];
}

void CFDCellKernel::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
	{
		_flux_old[q] =  _flux[q];
		_viscous_old[q] =  _viscous[q];
	}

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
			_flux_jacobi_variable[p][q] = (_flux[p] - _flux_old[p])/_perturbation;

		_cfd_data.uh[q] -= _perturbation;
	}

	_cfd_data.reinitInviscous();
	for (int q = 0; q < _n_equation; ++q)
	for (int beta = 0; beta < 3; ++beta)
	{
		_cfd_data.duh[q](beta) += _perturbation;
		reinitViscous();
		for (int p = 0; p < _n_equation; ++p)
		for (int alpha = 0; alpha< 3; ++alpha)
//			_flux_jacobi_grad_variable[p][q](alpha, beta) = (_flux[p](alpha) - _flux_old[p](alpha))/_perturbation;
			_flux_jacobi_grad_variable[p][q](alpha, beta) = (_viscous[p](alpha) - _viscous_old[p](alpha))/_perturbation;

		_cfd_data.duh[q](beta) -= _perturbation;
	}
}

Real CFDCellKernel::computeQpJacobian(unsigned int p, unsigned int q)
{
	Real r = (_flux_jacobi_variable[p][q]*_phi[_j][_qp]+_flux_jacobi_grad_variable[p][q]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
	return -r ;
}

ErrorVector & CFDCellKernel::getErrorVector(std::string indicator)
{
  return _fe_problem.adaptivity().getErrorVector(indicator);
}


