//
//#include "KOmegaCellKernel.h"
//
//template<>
//InputParameters validParams<KOmegaCellKernel>()
//{
//	  InputParameters params = validParams<MultiKernel>();
//	  params += validParams<CFDBase>();
//
//	  return params;
//}
//KOmegaCellKernel::KOmegaCellKernel(const InputParameters & parameters):
//		MultiKernel(parameters),
//		KOmegaModelBase(parameters)
//{
//}
//
//Real KOmegaCellKernel::computeQpResidual(unsigned int p)
//{
//	return 0;
//}
//
//Real KOmegaCellKernel::computeQpJacobian(unsigned int p, unsigned int q)
//{
//	return 0;
//}
//
//void KOmegaCellKernel::precalculateResidual()
//{
//	Real uh[10];
//	RealGradient duh[10];
//	RealVectorValue vis_term[10];
//
//	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//	{
//		valueAtCellPoint(uh);
//		valueGradAtCellPoint(duh);
//
//		inviscousTerm(_flux_vector[_qp], uh);
//		viscousTerm(vis_term, uh, duh);
//		for (unsigned int eq = 0; eq < _n_equation; ++eq)
//		{
//			_flux_vector[_qp][eq] -= vis_term[eq];
//		}
//		sourceTerm(_source_term[_qp], uh, duh);
//	}
//}
