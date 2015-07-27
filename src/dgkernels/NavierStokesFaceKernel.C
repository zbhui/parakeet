
#include "NavierStokesFaceKernel.h"


template<>
InputParameters validParams<NavierStokesFaceKernel>()
{
  InputParameters params = validParams<EulerFaceKernel>();

  return params;
}

NavierStokesFaceKernel::NavierStokesFaceKernel(const std::string & name, InputParameters parameters):
		EulerFaceKernel(name, parameters)
{
}

Real NavierStokesFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	if(type == Moose::Element)
		return _flux[_qp][_eq] * _test[_i][_qp] + 0.5*_epsilon * _penalty[_qp][_eq]*_grad_test[_i][_qp];

	if(type == Moose::Neighbor)
		return -_flux[_qp][_eq] * _test_neighbor[_i][_qp] + 0.5*_epsilon*_penalty_neighbor[_qp][_eq]*_grad_test_neighbor[_i][_qp];

	mooseError("face flux error.");
	return 0.;
}

void NavierStokesFaceKernel::precalculateResidual()
{
	mooseAssert(_n_equation < 10, "multiFaceKernel方程个数应<10");
	mooseAssert(_qrule->n_points() < 40, "mulitFaceKernel积分点个数应<40");

	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;

	Real ul[10], ur[10];
	RealGradient dul[10], dur[10];
	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		valueAtLeftFace(ul);
		valueAtRightFace(ur);
		valueGradAtLeftFace(dul);
		valueGradAtRightFace(dur);

		Point normal = _normals[_qp];
		fluxTerm(_flux[_qp], ul, ur, dul, dur);
		penaltyTerm(_penalty[_qp], ul, ur);
		penaltyTermNeighbor(_penalty_neighbor[_qp], ul, ur);
		for (int p = 0; p < _n_equation; ++p)
			_flux[_qp][p] += _sigma/h_elem*(_penalty[_qp][p] + _penalty_neighbor[_qp][p])/2.*normal;

		computeJacobianVariable(ul, ur, dul, dur);
		computeJacobianGradient(ul, ur, dul, dur);
	}
}

Real NavierStokesFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	  Real r = 0;

	  switch (type)
	  {
	  case Moose::ElementElement:
		  r = _jacobi_variable_ee[_qp][_ep][_eq]*_phi[_j][_qp]*_test[_i][_qp];
		  r -= _jacobi_gradient_ee[_qp][_ep][_eq]*_grad_phi[_j][_qp]*_test[_i][_qp];
		  r += 0.5*_epsilon*_jacobi_penalty_variable_ee[_qp][_ep][_eq]*_phi[_j][_qp]*_grad_test[_i][_qp];
	    break;

	  case Moose::ElementNeighbor:
		  r = _jacobi_variable_en[_qp][_ep][_eq]*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		  r -= _jacobi_gradient_en[_qp][_ep][_eq]*_grad_phi_neighbor[_j][_qp]*_test[_i][_qp];
		  r += 0.5*_epsilon*_jacobi_penalty_variable_en[_qp][_ep][_eq]*_phi_neighbor[_j][_qp]*_grad_test[_i][_qp];
	    break;

	  case Moose::NeighborElement:
		  r = _jacobi_variable_ne[_qp][_ep][_eq]*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		  r -= _jacobi_gradient_ne[_qp][_ep][_eq]*_grad_phi[_j][_qp]*_test_neighbor[_i][_qp];
		  r += 0.5*_epsilon*_jacobi_penalty_variable_ne[_qp][_ep][_eq]*_phi[_j][_qp]*_grad_test_neighbor[_i][_qp];
	    break;

	  case Moose::NeighborNeighbor:
		  r = _jacobi_variable_nn[_qp][_ep][_eq]*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		  r -= _jacobi_gradient_nn[_qp][_ep][_eq]*_grad_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		  r += 0.5*_epsilon*_jacobi_penalty_variable_nn[_qp][_ep][_eq]*_phi_neighbor[_j][_qp]*_grad_test_neighbor[_i][_qp];
	    break;
	  }

//	  return 0.;
	  return r;
}

Real NavierStokesFaceKernel::computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar)
{
	return 0.;
}

void NavierStokesFaceKernel::fluxTerm(Real* flux, Real* ul, Real* ur, RealGradient* dul, RealGradient* dur)
{
	RealVectorValue vfl[10], vfr[10];
	Point normal = _normals[_qp];

	viscousTerm(vfl, ul, dul);
	viscousTerm(vfr, ur, dur);
	fluxRiemann(flux, ul, ur, normal);
	for (int q = 0; q < _n_equation; ++q)
		flux[q] -= (vfl[q]+vfr[q])*normal/2.;
}

void NavierStokesFaceKernel::penaltyTerm(RealVectorValue* penalty, Real* ul, Real* ur)
{
	RealGradient duh[10];
	Point normal = _normals[_qp];
	for (int q = 0; q < _n_equation; ++q)
		duh[q] = (ul[q]-ur[q])*normal;

	viscousTerm(penalty, ul, duh);
}

void NavierStokesFaceKernel::penaltyTermNeighbor(RealVectorValue* penalty_neighbor, Real* ul, Real* ur)
{
	RealGradient duh[10];
	Point normal = _normals[_qp];
	for (int q = 0; q < _n_equation; ++q)
		duh[q] = (ul[q]-ur[q])*normal;

	viscousTerm(penalty_neighbor, ur, duh);
}

void NavierStokesFaceKernel::computeJacobianVariable(Real* ul, Real* ur, RealGradient* dul, RealGradient* dur)
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;

	Real flux_new[10];
	RealVectorValue penalty_new[10], penalty_neighbor_new[10];
	Point normal = _normals[_qp];
	for (int q = 0; q < _n_equation; ++q)
	{
		ul[q] += _ds;
		fluxTerm(flux_new, ul, ur, dul, dur);
		penaltyTerm(penalty_new, ul, ur);
		penaltyTermNeighbor(penalty_neighbor_new, ul, ur);
		for (int p = 0; p < _n_equation; ++p)
			flux_new[p] += _sigma/h_elem*(penalty_new[p] + penalty_neighbor_new[p])/2.*normal;

		for (int p = 0; p < _n_equation; ++p)
		{
			Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
			_jacobi_variable_ee[_qp][p][q] = tmp;
			_jacobi_variable_ne[_qp][p][q] = -tmp;

			_jacobi_penalty_variable_ee[_qp][p][q] = (penalty_new[p] -  _penalty[_qp][p])/_ds;
			_jacobi_penalty_variable_ne[_qp][p][q] = (penalty_neighbor_new[p] -  _penalty_neighbor[_qp][p])/_ds;
		}
		ul[q] -= _ds;

		ur[q] += _ds;
		fluxTerm(flux_new, ul, ur, dul, dur);
		penaltyTerm(penalty_new, ul, ur);
		penaltyTermNeighbor(penalty_neighbor_new, ul, ur);
		for (int p = 0; p < _n_equation; ++p)
			flux_new[p] += _sigma/h_elem*(penalty_new[p] + penalty_neighbor_new[p])/2.*normal;

		for (int p = 0; p < _n_equation; ++p)
		{
			Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
			_jacobi_variable_en[_qp][p][q] = tmp;
			_jacobi_variable_nn[_qp][p][q] = -tmp;

			_jacobi_penalty_variable_en[_qp][p][q] = (penalty_new[p] -  _penalty[_qp][p])/_ds;
			_jacobi_penalty_variable_nn[_qp][p][q] = (penalty_neighbor_new[p] -  _penalty_neighbor[_qp][p])/_ds;
		}
		ur[q] -= _ds;
	}
}

void NavierStokesFaceKernel::computeJacobianGradient(Real* ul, Real* ur, RealGradient* dul, RealGradient* dur)
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_neighbor_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;

	Real flux_new[10];
	RealVectorValue penalty_new[10], penalty_neighbor_new[10];
	RealVectorValue vis_term_left[10], vis_term_right[10], vis_term_new[10];
	Point normal = _normals[_qp];

	viscousTerm(vis_term_left, ul, dul);
	viscousTerm(vis_term_right, ur, dur);
	for (int beta = 0; beta < 3; ++beta)
	for (int q = 0; q < _n_equation; ++q)
	{
		dul[q](beta) += _ds;
		viscousTerm(vis_term_new, ul, dul);
		for (int p = 0; p < _n_equation; ++p)
		{
			Real tmp = (vis_term_new[p] - vis_term_left[p])/_ds * normal;
			_jacobi_gradient_ee[_qp][p][q](beta) = 0.5* tmp;
			_jacobi_gradient_ne[_qp][p][q](beta) = -0.5*tmp;
		}
		dul[q](beta) -= _ds;

		dur[q](beta) += _ds;
		viscousTerm(vis_term_new, ur, dur);
		for (int p = 0; p < _n_equation; ++p)
		{
			Real tmp = (vis_term_new[p] - vis_term_right[p])/_ds * normal;
			_jacobi_gradient_en[_qp][p][q](beta) = 0.5*tmp;
			_jacobi_gradient_nn[_qp][p][q](beta) = -0.5*tmp;
		}
		dur[q](beta) -= _ds;
	}
}
