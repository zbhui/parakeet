
#include "KOmegaFaceKernel.h"


template<>
InputParameters validParams<KOmegaFaceKernel>()
{
  InputParameters params = validParams<MultiDGKernel>();
  params += validParams<KOmegaModelBase>();

  return params;
}

KOmegaFaceKernel::KOmegaFaceKernel(const InputParameters & parameters):
		MultiDGKernel(parameters),
		KOmegaModelBase(parameters)
{
}

Real KOmegaFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	if(type == Moose::Element)
	{
		return _flux[_qp][_eq] * _test[_i][_qp] + 0.5*_epsilon * _penalty[_qp][_eq]*_grad_test[_i][_qp];
//		return 0.;
	}
	if(type == Moose::Neighbor)
	{
		return -_flux[_qp][_eq] * _test_neighbor[_i][_qp] + 0.5*_epsilon*_penalty_neighbor[_qp][_eq]*_grad_test_neighbor[_i][_qp];
//		return 0.;
	}

	mooseError("face flux error.");
	return 0.;
}

void KOmegaFaceKernel::precalculateResidual()
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;

	Real ul[10], ur[10], uh_bar[10];
	RealGradient dul[10], dur[10], duh[10];
	RealVectorValue vfl[10], vfr[10];

	mooseAssert(_n_equation < 10, "multiFaceKernel方程个数应<10");
	mooseAssert(_qrule->n_points() < 40, "mulitFaceKernel积分点个数应<40");

	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		valueAtLeftFace(ul);
		valueAtRightFace(ur);
		valueGradAtLeftFace(dul);
		valueGradAtRightFace(dur);

		Point normal = _normals[_qp];
		viscousTerm(vfl, ul, dul);
		viscousTerm(vfr, ur, dur);
		fluxRiemann(_flux[_qp], ul, ur, normal);

		for (_eq = 0; _eq < _n_equation; ++_eq)
		{
			uh_bar[_eq] = (ul[_eq]+ur[_eq])/2.;
			duh[_eq] = (ul[_eq]-ur[_eq])*normal;
		}

		viscousTerm(_penalty[_qp], ul, duh);
		viscousTerm(_penalty_neighbor[_qp], ur, duh);

		for (_eq = 0; _eq < _n_equation; ++_eq)
		{
			_flux[_qp][_eq] -= ((vfl[_eq]+vfr[_eq]) - _sigma/h_elem*(_penalty[_qp][_eq] + _penalty_neighbor[_qp][_eq]))/2.*normal;
		}
	}
}

void KOmegaFaceKernel::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
{
	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real rho, u, v, w, pre;
	rho = (ul[0] + ur[0])/2.;
	u = (ul[1] + ur [1])/rho/2.;
	v = (ul[2] + ur [2])/rho/2.;
	w = (ul[3] + ur [3])/rho/2.;
	pre = (pressure(ul) + pressure(ur))/2.;
	Real lam = fabs(u*normal(0) + v * normal(1) + w * normal(2)) + sqrt(_gamma*pre/rho);
	for (int eq = 0; eq < _n_equation; ++eq)
	{
		flux[eq] = 0.5*((fl[eq] + fr[eq])*normal) + lam*(ul[eq] - ur[eq]);
	}
}

Real KOmegaFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	return 0.;
}

Real KOmegaFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	return 0.;
}
