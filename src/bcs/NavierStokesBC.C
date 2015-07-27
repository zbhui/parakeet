#include "NavierStokesBC.h"

// N-S方程边界条件
template<>
InputParameters validParams<NavierStokesBC>()
{
	InputParameters params = validParams<EulerBC>();
	return params;
}
NavierStokesBC::NavierStokesBC(const InputParameters & parameters):
		EulerBC(parameters)
{
}
void NavierStokesBC::precalculateResidual()
{
	mooseAssert(_n_equation < 10, "multiBC方程个数应<10");
	mooseAssert(_qrule->n_points() < 40, "mulitBC积分点个数应<40");

	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_current_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;

	Real ul[10], ur[10];
	RealGradient dul[10], dur[10], duh[10];
	RealVectorValue vfl[10], vfr[10];
	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		valueAtLeftFace(_ul[_qp]);
		valueAtRightFace(_ur[_qp]);
		valueAtLeftFace(ul);
		valueAtRightFace(ur);
		valueGradAtLeftFace(dul);
		valueGradAtRightFace(dur);

		Point normal = _normals[_qp];

		viscousTerm(vfl, ul, dul);
		viscousTerm(vfr, ur, dur);
		fluxRiemann(_flux[_qp], ul, ur, normal);

		for (int eq = 0; eq < _n_equation; ++eq)
			duh[_eq] = (ul[_eq]-ur[_eq])*normal;

		viscousTerm(_penalty[_qp], ul, duh);
		viscousTerm(_penalty_neighbor[_qp], ur, duh);

		for (_eq = 0; _eq < _n_equation; ++_eq)
			_flux[_qp][_eq] -= ((vfl[_eq]+vfr[_eq]) - _sigma/h_elem*(_penalty[_qp][_eq] + _penalty_neighbor[_qp][_eq]))/2.*normal;
	}
}

void NavierStokesBC::computeJacobianVariable(Real* ul)
{
//	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
//	const double h_elem = (_current_elem_volume+_current_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;
//
//	Real flux_new[10];
//	RealVectorValue penalty_new[10], penalty_neighbor_new[10];
//	Point normal = _normals[_qp];
//	for (int q = 0; q < _n_equation; ++q)
//	{
//		ul[q] += _ds;
//		fluxTerm(flux_new, ul, ur, dul, dur);
//		penaltyTerm(penalty_new, ul, ur);
//		penaltyTermNeighbor(penalty_neighbor_new, ul, ur);
//		for (int p = 0; p < _n_equation; ++p)
//			flux_new[p] += _sigma/h_elem*(penalty_new[p] + penalty_neighbor_new[p])/2.*normal;
//
//		for (int p = 0; p < _n_equation; ++p)
//		{
//			Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
//			_jacobi_variable_ee[_qp][p][q] = tmp;
//			_jacobi_variable_ne[_qp][p][q] = -tmp;
//
//			_jacobi_penalty_variable_ee[_qp][p][q] = (penalty_new[p] -  _penalty[_qp][p])/_ds;
//			_jacobi_penalty_variable_ne[_qp][p][q] = (penalty_neighbor_new[p] -  _penalty_neighbor[_qp][p])/_ds;
//		}
//		ul[q] -= _ds;
//
//		ur[q] += _ds;
//		fluxTerm(flux_new, ul, ur, dul, dur);
//		penaltyTerm(penalty_new, ul, ur);
//		penaltyTermNeighbor(penalty_neighbor_new, ul, ur);
//		for (int p = 0; p < _n_equation; ++p)
//			flux_new[p] += _sigma/h_elem*(penalty_new[p] + penalty_neighbor_new[p])/2.*normal;
//
//		for (int p = 0; p < _n_equation; ++p)
//		{
//			Real tmp = (flux_new[p] - _flux[_qp][p])/_ds;
//			_jacobi_variable_en[_qp][p][q] = tmp;
//			_jacobi_variable_nn[_qp][p][q] = -tmp;
//
//			_jacobi_penalty_variable_en[_qp][p][q] = (penalty_new[p] -  _penalty[_qp][p])/_ds;
//			_jacobi_penalty_variable_nn[_qp][p][q] = (penalty_neighbor_new[p] -  _penalty_neighbor[_qp][p])/_ds;
//		}
//		ur[q] -= _ds;
//	}
}

Real NavierStokesBC::computeQpResidual()
{
	return _flux[_qp][_eq] * _test[_i][_qp] + 0.5*_epsilon * _penalty[_qp][_eq]* _grad_test[_i][_qp];
}

Real NavierStokesBC::computeQpJacobian()
{
	return 0.;
}

Real NavierStokesBC::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.;
}

void NavierStokesBC::wallBC(Real* ur)
{
	Real ui[5];
	valueAtLeftFace(ui);
    Real pre = pressure(ui);
    Real twall = 1.;
//	Real pre = ui[0]*twall/_gamma/_mach/_mach;

//    ur[ 0 ] = ui[0];
    ur[0] = _gamma*_mach*_mach*pre/twall;
    ur[1] = 0.;
    ur[2] = 0.;
    ur[3] = 0.;
    ur[4] = pre/(_gamma-1) + 0.5*( ur[1]*ur[1] + ur[2]*ur[2] + ur[3]*ur[3] )/ur[0];
}

void NavierStokesBC::farFieldBC(Real* ur)
{
	EulerBC::farFieldBC(ur);
}

void NavierStokesBC::symmetricBC(Real* ur)
{
	EulerBC::wallBC(ur);
}

void NavierStokesBC::pressureOutBC(Real* ur)
{
	Real pR = 1 /_gamma/_mach/_mach;

	Real ui[5];
	valueAtLeftFace(ui);
	ur[0] = ui[0];
	ur[1] = ui[1];
	ur[2] = ui[2];
	ur[3] = ui[3];
	ur[4] = pR/(_gamma-1) + 0.5*( ur[1]*ur[1] + ur[2]*ur[2] + ur[3]*ur[3] )/ur[0];;
}
