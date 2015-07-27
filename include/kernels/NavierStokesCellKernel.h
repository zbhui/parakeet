#pragma once

#include "EulerCellKernel.h"

// 前置声明
class NavierStokesCellKernel;

template<>
InputParameters validParams<NavierStokesCellKernel>();

class NavierStokesCellKernel :
public EulerCellKernel
{
public:
	NavierStokesCellKernel(const InputParameters & parameters);
	virtual ~NavierStokesCellKernel(){}

protected:
	RealTensorValue _jacobi_gradient[40][10][10];

	virtual Real computeQpResidual(unsigned int p);
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);

	virtual void precalculateResidual();

	void fluxTerm(RealVectorValue *flux_vector, Real *uh, RealGradient *duh);
	void computeJacobianVariable(Real *uh, RealGradient *duh);
	void computeJacobianGradient(Real *uh, RealGradient *duh);
};
