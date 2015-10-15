
#pragma once

#include "DGKernel.h"
#include "MultiVariableInterface.h"

class MultiDGKernel :
  public DGKernel,
  public MultiVariableInterface
{
public:
	MultiDGKernel(const InputParameters & parameters);

	virtual ~MultiDGKernel(){};

	virtual void computeResidual();
	virtual void computeElemNeighResidual(Moose::DGResidualType type);
	virtual void  precalculateResidual() = 0;
	virtual Real computeQpResidual(Moose::DGResidualType type, unsigned int p) = 0;

	virtual void computeOffDiagJacobian(unsigned int jvar);
	virtual void computeJacobian();
	virtual void computeElemNeighJacobian(Moose::DGJacobianType type);
	virtual void precalculateJacobian() = 0;
	virtual Real computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q) = 0;

protected:
	FEProblem & _fe_problem;
	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
};


template<>
InputParameters validParams<MultiDGKernel>();

