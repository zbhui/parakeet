#pragma once

#include "MultiKernel.h"
#include "CFDBase.h"

// 前置声明
class EulerCellKernel;

template<>
InputParameters validParams<EulerCellKernel>();

class EulerCellKernel :
public MultiKernel,
public CFDBase
{
public:
	EulerCellKernel(const InputParameters & parameters);
	virtual ~EulerCellKernel(){}

protected:
	RealVectorValue _flux[40][10];
	RealVectorValue _jacobi_variable[40][10][10];

	virtual void precalculateResidual();
	virtual void precalculateJacobian();

	virtual Real computeQpResidual(unsigned int p);
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);
};
