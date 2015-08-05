
#include "CFDFaceKernel.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDFaceKernel>()
{
	  InputParameters params = validParams<MultiDGKernel>();
	  params.addParam<Real>("perturbation", 1E-08, "有限差分求Jacobian矩阵的变量增量");

	  return params;
}
CFDFaceKernel::CFDFaceKernel(const InputParameters & parameters):
		MultiDGKernel(parameters),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
		_cfd_data(_cfd_problem),
		_cfd_data_neighbor(_cfd_problem),
		_lift_data(_cfd_problem),
		_perturbation(getParam<Real>("perturbation")),
		_penalty(0)
{
}

void CFDFaceKernel::reinit()
{
	_penalty = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume /2.;
	_penalty = (_var_order*_var_order+1)/_penalty;

	_cfd_data.reinit();
	_cfd_data_neighbor.reinit();

	liftOperator();
	fluxRiemann();
}

void CFDFaceKernel::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data.duh[i] = (*_grad_uh[i])[_qp];
		_cfd_data_neighbor.uh[i] = (*_uh_neighbor[i])[_qp];
		_cfd_data_neighbor.duh[i] = (*_grad_uh_neighbor[i])[_qp];
	}

	reinit();
}

void CFDFaceKernel::fluxRiemann()
{
	Real lam = fabs(_cfd_data.vel*_normals[_qp]) + _cfd_data.c;
	lam += fabs(_cfd_data_neighbor.vel*_normals[_qp]) + _cfd_data_neighbor.c;
	lam /= 2.;
	for (int p = 0; p < _n_equation; ++p)
	{
		_flux[p] = 0.5*(_cfd_data.invis_flux[p] + _cfd_data_neighbor.invis_flux[p])*_normals[_qp] +
				    lam*(_cfd_data.uh[p] - _cfd_data_neighbor.uh[p]);
	}
	for (int p = 0; p < _n_equation; ++p)
		_flux[p] -= 0.5*((_cfd_data.vis_flux[p] + _cfd_data_neighbor.vis_flux[p])-_penalty*_lift[p])*_normals[_qp];

}

void CFDFaceKernel::liftOperator()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_lift_data.uh[i] = (_cfd_data.uh[i] + _cfd_data_neighbor.uh[i])/2.;
		_lift_data.duh[i] = (_cfd_data.uh[i] - _cfd_data_neighbor.uh[i])*_normals[_qp];
	}

	_lift_data.reinit();
	for (int p = 0; p < _n_equation; ++p)
		_lift[p] = _lift_data.vis_flux[p];
}

Real CFDFaceKernel::computeQpResidual(Moose::DGResidualType type, unsigned int p)
{
	if(type == Moose::Element)
	{
		return  _flux[p] * _test[_i][_qp] + _lift[p]*_grad_test[_i][_qp];
	}
	if(type == Moose::Neighbor)
	{
		return -_flux[p] * _test_neighbor[_i][_qp] + _lift[p]*_grad_test_neighbor[_i][_qp];
	}
	mooseError("face flux error.");
	return 0.;
}


void CFDFaceKernel::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
	{
		_flux_old[q] = _flux[q];
		_lift_old[q] = _lift[q];
	}

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux_jacobi_variable_ee[p][q] = (_flux[p] - _flux_old[p])/_perturbation;
			_lift_jacobi_variable_ee[p][q] = (_lift[p] - _lift_old[p])/_perturbation;
		}
		_cfd_data.uh[q] -= _perturbation;
	}
	for (int q = 0; q < _n_equation; ++q)
	for (int beta = 0; beta < 3; ++beta)
	{
		_cfd_data.duh[q](beta) += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux_jacobi_grad_variable_ee[p][q](beta) = (_flux[p] - _flux_old[p])/_perturbation;
		}
		_cfd_data.duh[q](beta) -= _perturbation;
	}

	_cfd_data.reinit();
	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data_neighbor.uh[q] += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux_jacobi_variable_en[p][q] = (_flux[p] - _flux_old[p])/_perturbation;
			_lift_jacobi_variable_en[p][q] = (_lift[p] - _lift_old[p])/_perturbation;
		}
		_cfd_data_neighbor.uh[q] -= _perturbation;
	}
	for (int q = 0; q < _n_equation; ++q)
	for (int beta = 0; beta < 3; ++beta)
	{
		_cfd_data_neighbor.duh[q](beta) += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux_jacobi_grad_variable_en[p][q](beta) = (_flux[p] - _flux_old[p])/_perturbation;
		}
		_cfd_data.duh[q](beta) -= _perturbation;
	}
}

Real CFDFaceKernel::computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q)
{
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _flux_jacobi_variable_ee[p][q]*_phi[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_ee[p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
		r += _lift_jacobi_variable_ee[p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _flux_jacobi_variable_en[p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _flux_jacobi_grad_variable_en[p][q]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		r += _lift_jacobi_variable_en[p][q]*_grad_test[_i][_qp]*_phi_neighbor[_j][_qp];
		break;

	case Moose::NeighborElement:
		r = -_flux_jacobi_variable_ee[p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_ee[p][q]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		r += _lift_jacobi_variable_ee[p][q]*_grad_test_neighbor[_i][_qp]*_phi[_j][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_flux_jacobi_variable_en[p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += -_flux_jacobi_grad_variable_en[p][q]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		r += _lift_jacobi_variable_en[p][q]*_grad_test_neighbor[_i][_qp]*_phi_neighbor[_j][_qp];
		break;
	}

//	std::cout << _flux_jacobi_grad_variable_en[p][q] <<std::endl;
	return r;
}

