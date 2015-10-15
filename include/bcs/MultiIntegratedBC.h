
#pragma once

#include "IntegratedBC.h"
#include "MultiVariableInterface.h"

class MultiIntegratedBC :
public IntegratedBC,
public MultiVariableInterface
{
public:
	  MultiIntegratedBC(const InputParameters & params);
	  virtual ~MultiIntegratedBC(){}

	  virtual void precalculateResidual() = 0;
	  virtual Real computeQpResidual(unsigned int p) = 0;
	  virtual void computeResidual();

	  void computeJacobianBlock(unsigned int jvar);
	  virtual void computeJacobian();
	  virtual void precalculateJacobian() = 0;
	  virtual Real computeQpJacobian(unsigned int p, unsigned int q) = 0;

protected:
	  virtual Real computeQpResidual();
	  virtual Real computeQpJacobian();
};

template<>
InputParameters validParams<MultiIntegratedBC>();
