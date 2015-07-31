//
//#pragma once
//
//#include "EulerBC.h"
//
//
//class NavierStokesBC :
//public EulerBC
//{
//public:
//	NavierStokesBC(const InputParameters & parameters);
//	virtual ~NavierStokesBC(){}
//
//
//protected:
//	  RealVectorValue _penalty[40][10];
//	  RealVectorValue _penalty_neighbor[40][10];
//
//	  virtual Real computeQpResidual(unsigned int p);
//	  virtual Real computeQpJacobian(unsigned int p, unsigned int q);
//
//	  virtual void precalculateResidual();
//	  virtual void precalculateJacobian();
//	  virtual void computeJacobianVariable(Real *ul);
//
//	  virtual void wall(Real *ur);
//	  virtual void farFieldBC(Real *ur);
//	  virtual void symmetricBC(Real *ur);
//	  virtual void pressureOutBC(Real *ur);
//};
//
//template<>
//InputParameters validParams<NavierStokesBC>();
