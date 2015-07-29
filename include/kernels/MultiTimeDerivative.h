///****************************************************************/
//
//#pragma once
//
//#include "MultiTimeKernel.h"
//
//// Forward Declaration
//class MultiTimeDerivative;
//
//template<>
//InputParameters validParams<MultiTimeDerivative>();
//
//class MultiTimeDerivative : public MultiTimeKernel
//{
//public:
//  MultiTimeDerivative(const InputParameters & parameters);
//
//  virtual void computeJacobian();
//
//protected:
//  virtual Real computeQpResidual();
//  virtual Real computeQpJacobian();
//
//  virtual void precalculateResidual();
//
//  Real _duh_dt[20][20];
//  bool _lumping;
//};
//
