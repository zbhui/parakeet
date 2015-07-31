
#include "EulerFaceKernel.h"

template<>
InputParameters validParams<EulerFaceKernel>()
{
	  InputParameters params = validParams<MultiDGKernel>();
	  params += validParams<CFDBase>();

	  return params;
}
EulerFaceKernel::EulerFaceKernel(const InputParameters & parameters):
		MultiDGKernel(parameters),
		CFDBase(parameters),
		_cfd_data(_mach, _reynolds, _gamma, _prandtl),
		_cfd_data_neighbor(_mach, _reynolds, _gamma, _prandtl)
{
}

void EulerFaceKernel::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data_neighbor.uh[i] = (*_uh_neighbor[i])[_qp];
	}

	_cfd_data.reinit();
	_cfd_data_neighbor.reinit();

	fluxRiemann();
}

void EulerFaceKernel::fluxRiemann()
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

Real EulerFaceKernel::computeQpResidual(Moose::DGResidualType type, unsigned int p)
{
	if(type == Moose::Element)
	{
		return _flux[p] * _test[_i][_qp];
	}
	if(type == Moose::Neighbor)
	{
		return -_flux[p] * _test_neighbor[_i][_qp];
	}
	mooseError("face flux error.");
	return 0.;
}


void EulerFaceKernel::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
		_flux_old[q] = _flux[q];

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _ds;
		_cfd_data.reinit();
		fluxRiemann();
		for (int p = 0; p < _n_equation; ++p)
		{
			_jacobi_variable_ee[p][q] = (_flux[p] - _flux_old[p])/_ds;
		}
		_cfd_data.uh[q] -= _ds;
	}
	_cfd_data.reinit();

	for (int q = 0; q < _n_equation; ++q)
	{

		_cfd_data_neighbor.uh[q] += _ds;
		_cfd_data_neighbor.reinit();
		fluxRiemann();
		for (int p = 0; p < _n_equation; ++p)
		{
			_jacobi_variable_en[p][q] = (_flux[p] - _flux_old[p])/_ds;
		}
		_cfd_data_neighbor.uh[q] -= _ds;
	}
	_cfd_data_neighbor.reinit();
}

Real EulerFaceKernel::computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q)
{
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _jacobi_variable_ee[p][q]*_phi[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _jacobi_variable_en[p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::NeighborElement:
		r = -_jacobi_variable_ee[p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_jacobi_variable_en[p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		break;
	}

	return r;
}

