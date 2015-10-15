
#pragma once

#include "Kernel.h"
#include "MultiVariableInterface.h"

class MultiKernel :
  public Kernel,
  public MultiVariableInterface
{
public:
  MultiKernel(const InputParameters & parameters);

  virtual ~MultiKernel(){};

  virtual void computeResidual();
  virtual void precalculateResidual() = 0;
  virtual Real computeQpResidual(unsigned int p) = 0;

  virtual void computeOffDiagJacobian(unsigned int jvar);
  virtual void computeJacobian();
  virtual void precalculateJacobian() = 0;
  virtual Real computeQpJacobian(unsigned int p, unsigned int q) = 0;

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);
};


template<>
InputParameters validParams<MultiKernel>();
