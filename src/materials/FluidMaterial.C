//#include "FluidMaterial.h"
//
//template<>
//InputParameters validParams<FluidMaterial>()
//{
//  InputParameters params = validParams<Material>();
//  params += validParams<CFDBase>();
//
//  params.addRequiredCoupledVar("variables", "守恒变量");
//
//  return params;
//}
//
//
//FluidMaterial::FluidMaterial(const InputParameters & parameters):
//		Material(parameters),
//		CFDBase(parameters),
//
//		_viscous_stress_tensor(declareProperty<RealTensorValue>("_viscous_stress_tensor")),
//		_test_viscosity(declareProperty<Real>("viscosity")),
//
//		_rho(coupledValue("variables", 0)),
//		_rhou(coupledValue("variables", 1)),
//		_rhov(coupledValue("variables", 2)),
//		_rhow(coupledValue("variables", 3)),
//		_rhoe(coupledValue("variables", 4)),
//
//		_grad_rho(coupledGradient("variables", 0)),
//		_grad_rhou(coupledGradient("variables", 1)),
//		_grad_rhov(coupledGradient("variables", 2)),
//		_grad_rhow(coupledGradient("variables", 3)),
//		_grad_rhoe(coupledGradient("variables", 4))
//
//
//{
//}
//
//
///**
// * Must be called _after_ the child class computes dynamic_viscosity.
// */
//void
//FluidMaterial::computeQpProperties()
//{
//	RealTensorValue grad_outter_u;
//	_viscous_stress_tensor[_qp] = grad_outter_u;
//
//	_test_viscosity[_qp] = 1.0;
//}
//
//
//
