
#include "CFDBC.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDBC>()
{
	InputParameters params = validParams<MultiIntegratedBC>();
	params.addParam<Real>("perturbation", 1E-08, "有限差分求Jacobian矩阵的变量增量");
	return params;
}

CFDBC::CFDBC(const InputParameters & parameters):
		MultiIntegratedBC(parameters),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
		_cfd_data(_cfd_problem),
		_cfd_data_neighbor(_cfd_problem),
		_perturbation(getParam<Real>("perturbation"))
{

}

void CFDBC::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data.duh[i] = (*_grad_uh[i])[_qp];
		_cfd_data_neighbor.duh[i] = (*_grad_uh[i])[_qp];
	}

	_cfd_data.reinit();
	boundaryCondition();
	_cfd_data_neighbor.reinit();

	fluxRiemann();
}

Real CFDBC::computeQpResidual(unsigned int p)
{
	return _flux[p] * _test[_i][_qp];
}

Real CFDBC::computeQpJacobian(unsigned int p, unsigned int q)
{
	return _jacobi_variable[p][q]*_phi[_j][_qp]*_test[_i][_qp];
}

void CFDBC::fluxRiemann()
{
	Real lam = fabs(_cfd_data.vel*_normals[_qp]) + _cfd_data.c;
	lam += fabs(_cfd_data_neighbor.vel*_normals[_qp]) + _cfd_data_neighbor.c;
	lam /= 2.;
	for (int p = 0; p < _n_equation; ++p)
	{
		_flux[p] = 0.5*(_cfd_data.invis_flux[p] + _cfd_data_neighbor.invis_flux[p])*_normals[_qp] +
				    lam*(_cfd_data.uh[p] - _cfd_data_neighbor.uh[p]);
	}
}

void CFDBC::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
		_flux_old[q] = _flux[q];

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _perturbation;
		_cfd_data.reinit();
		boundaryCondition();
		_cfd_data_neighbor.reinit();
		fluxRiemann();
		for (int p = 0; p < _n_equation; ++p)
			_jacobi_variable[p][q] = (_flux[p] - _flux_old[p])/_perturbation;

		_cfd_data.uh[q] -= _perturbation;
	}
//	_cfd_data.reinit();
//	_cfd_data_neighbor.reinit();
}

void CFDBC::boundaryCondition()
{
	for(int p = 0; p < _n_equation; ++p)
		_cfd_data_neighbor.uh[p] = _cfd_problem.boundaryCondition(_t, _q_point[_qp], p);
}
