//#pragma once
//
//#include "Material.h"
//#include "CFDBase.h"
//
//class FluidMaterial;
//
//template<>
//InputParameters validParams<FluidMaterial>();
//
///**
// * @warning Material 为单元积分点上的材料属性提供接口，但不能为面积分提供帮助
// */
//class FluidMaterial :
//public Material,
//public CFDBase
//{
//public:
//	FluidMaterial(const InputParameters & parameters);
//
//protected:
//  virtual void computeQpProperties();
//
//  MaterialProperty<RealTensorValue>& _viscous_stress_tensor;
//  MaterialProperty<Real>&            _test_viscosity;
////  MaterialProperty<Real>&            _dynamic_viscosity;
//
//protected:
//	VariableValue& _rho;
//	VariableValue& _rhou;
//	VariableValue& _rhov;
//	VariableValue& _rhow;
//	VariableValue& _rhoe;
//
//	// Gradients
//	VariableGradient& _grad_rho;
//	VariableGradient& _grad_rhou;
//	VariableGradient& _grad_rhov;
//	VariableGradient& _grad_rhow;
//	VariableGradient& _grad_rhoe;
//
//};
//
