//
//#pragma once
//
//#include "NavierStokesBC.h"
//#include "CouetteFlowIC.h"
//
//class CouetteFlowBC;
//
//template<>
//InputParameters validParams<CouetteFlowBC>();
//
//class CouetteFlowBC :
//public NavierStokesBC
//{
//public:
//	CouetteFlowBC(const InputParameters & params);
//	virtual ~CouetteFlowBC(){}
//
//
//protected:
//	virtual void valueAtRightFace(Real *ur);
//
//private:
//
//  Real density(const Point &p);
//  Real x_momentum(const Point &p);
//  Real y_momentum(const Point &p);
//  Real z_momentum(const Point &p);
//  Real total_energy(const Point &p);
//  Real temperature(const Point &p);
//};
