//
//#pragma once
//
//#include "KOmegaModelBase.h"
//#include "CFDBC.h"
//
//class KOmegaBC;
//
//template<>
//InputParameters validParams<KOmegaBC>();
//
//class KOmegaBC :
//public CFDBC,
//public KOmegaModelBase
//{
//public:
//	KOmegaBC(const InputParameters & params);
//	virtual ~KOmegaBC(){}
//
//protected:
////	Real _tu_infty;
////	Real _r_mu;
//	Real _flux[40][10];
//	RealVectorValue _penalty[40][10];
//	RealVectorValue _penalty_neighbor[40][10];
//
//	void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
//	virtual Real computeQpResidual();
//	virtual Real computeQpJacobian();
//	virtual Real computeQpOffDiagJacobian(unsigned int jvar);
//
//	virtual void precalculateResidual();
//
//	virtual void wallBC(Real *ur);
//	virtual void farFieldBC(Real *ur);
//	virtual void symmetricBC(Real *ur);
//};
