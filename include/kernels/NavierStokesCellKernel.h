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
	NavierStokesCellKernel(const std::string & name, InputParameters parameters);
	virtual ~NavierStokesCellKernel(){}

protected:
	RealTensorValue _jacobi_gradient[40][10][10];

	virtual Real computeQpJacobian();
	virtual Real computeQpOffDiagJacobian(unsigned int jvar);

	virtual void precalculateResidual();

	void fluxTerm(RealVectorValue *flux_vector, Real *uh, RealGradient *duh);
	void computeJacobianVariable(Real *uh, RealGradient *duh);
	void computeJacobianGradient(Real *uh, RealGradient *duh);
};
