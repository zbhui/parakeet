//#include "NavierStokesCellKernel.h"
//
//template<>
//InputParameters validParams<NavierStokesCellKernel>()
//{
//  InputParameters params = validParams<EulerCellKernel>();
//
//  return params;
//}
//
//NavierStokesCellKernel::NavierStokesCellKernel(const InputParameters & parameters):
//		EulerCellKernel(parameters)
//{
//}
//
//void NavierStokesCellKernel::precalculateResidual()
//{
//	Real uh[10];
//	RealGradient duh[10];
//	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//	{
//		valueAtCellPoint(uh);
//		valueGradAtCellPoint(duh);
//		fluxTerm(_flux[_qp], uh, duh);
//
//		computeJacobianVariable(uh, duh);
//		computeJacobianGradient(uh, duh);
//	}
//}
//
//Real NavierStokesCellKernel::computeQpResidual(unsigned int p)
//{
//	return 0;
//}
//
//Real NavierStokesCellKernel::computeQpJacobian(unsigned int p, unsigned int q)
//{
//	return -(_jacobi_variable[_qp][p][q]*_phi[_j][_qp] - _jacobi_gradient[_qp][p][q]*_grad_phi[_j][_qp])*_grad_test[_i][_qp];
//}
//
//void NavierStokesCellKernel::fluxTerm(RealVectorValue* flux_vector, Real* uh, RealGradient* duh)
//{
//	RealVectorValue vis_term[10];
//	inviscousTerm(flux_vector, uh);
//	viscousTerm(vis_term, uh, duh);
//	for (int eq = 0; eq < _n_equation; ++eq)
//		flux_vector[eq] -= vis_term[eq];
//}
//
//void NavierStokesCellKernel::computeJacobianVariable(Real* uh,	RealGradient* duh)
//{
//	RealVectorValue flux_vector_new[10];
//	for (int q = 0; q < _n_equation; ++q)
//	{
//		uh[q] += _ds;
//		fluxTerm(flux_vector_new, uh, duh);
//		for (int p = 0; p < _n_equation; ++p)
//			_jacobi_variable[_qp][p][q] = (flux_vector_new[p] - _flux[_qp][p])/_ds;
//
//		uh[q] -= _ds;
//	}
//}
//
//void NavierStokesCellKernel::computeJacobianGradient(Real* uh, RealGradient* duh)
//{
//	RealVectorValue vis_term[10], vis_term_new[10];
//	viscousTerm(vis_term, uh, duh);
//	for (int q = 0; q < _n_equation; ++q)
//		for (int beta = 0; beta < 3; ++beta)
//		{
//			duh[q](beta) += _ds;
//			viscousTerm(vis_term_new, uh, duh);
//			for (int p = 0; p < _n_equation; ++p)
//				for (int alpha = 0; alpha < 3; ++alpha)
//					_jacobi_gradient[_qp][p][q](alpha, beta) = (vis_term_new[p](alpha) - vis_term[p](alpha))/_ds;
//
//			duh[q](beta) -= _ds;
//		}
//}
