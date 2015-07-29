
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
		CFDBase(parameters)
{
}

void EulerFaceKernel::precalculateResidual()
{
	Real ul[10], ur[10];
	valueAtLeftFace(ul);
	valueAtRightFace(ur);
	Point normal = _normals[_qp];
	fluxRiemann(_flux, ul, ur, normal);
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
	Real flux_new[10], flux[10];
	Real ul[10], ur[10];

	valueAtLeftFace(ul);
	valueAtRightFace(ur);
	Point normal = _normals[_qp];
	fluxRiemann(flux, ul, ur, normal);
	for (int q = 0; q < _n_equation; ++q)
	{
		ul[q] += _ds;
		fluxRiemann(flux_new, ul, ur, normal);
		for (int p = 0; p < _n_equation; ++p)
		{
			Real tmp = (flux_new[p] - flux[p])/_ds;
			_jacobi_variable_ee[_qp][p][q] = tmp;
		}
		ul[q] -= _ds;

		ur[q] += _ds;
		fluxRiemann(flux_new, ul, ur, normal);
		for (int p = 0; p < _n_equation; ++p)
		{
			Real tmp = (flux_new[p] - flux[p])/_ds;
			_jacobi_variable_en[_qp][p][q] = tmp;
		}
		ur[q] -= _ds;
	}
}

void EulerFaceKernel::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
{
	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real rho, u, v, w, pre;
	rho = (ul[0] + ur[0])/2.;
	u = (ul[1] + ur [1])/rho/2;
	v = (ul[2] + ur [2])/rho/2;
	w = (ul[3] + ur [3])/rho/2;
	pre = (pressure(ul) + pressure(ur))/2.;
	Real lam = fabs(u*normal(0) + v * normal(1) + w * normal(2)) + sqrt(_gamma*pre/rho);
	for (int eq = 0; eq < _n_equation; ++eq)
	{
		flux[eq] = 0.5*(fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]);
	}

}

Real EulerFaceKernel::computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q)
{
	Real r = 0;
	switch (type)
	{
	case Moose::ElementElement:
		r = _jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::ElementNeighbor:
		r = _jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		break;

	case Moose::NeighborElement:
		r = -_jacobi_variable_ee[_qp][p][q]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		break;

	case Moose::NeighborNeighbor:
		r = -_jacobi_variable_en[_qp][p][q]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		break;
	}

	return r;
}

