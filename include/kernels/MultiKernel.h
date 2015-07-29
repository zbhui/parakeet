
#pragma once

#include "Kernel.h"

class MultiKernel :
  public Kernel
{
public:
	MultiKernel(const InputParameters & parameters);

  virtual ~MultiKernel(){};

  virtual void computeResidual();
  virtual void precalculateResidual() = 0;
  virtual Real computeQpResidual(unsigned int p) = 0;

  virtual void computeJacobian();
  virtual void precalculateJacobian() = 0;
  virtual Real computeQpJacobian(unsigned int p, unsigned int q) = 0;


  virtual Real computeQpResidual();
  virtual void computeOffDiagJacobian(unsigned int jvar);
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);

protected:

  virtual void valueAtCellPoint(Real *uh);
  virtual void valueGradAtCellPoint(RealGradient *duh);
  unsigned int _n_equation;
  std::vector<NonlinearVariableName> _variables;
  std::vector<VariableValue*> _uh;
  std::vector<VariableGradient*> _grad_uh;
  std::vector<VariableValue *> _uh_dot;
  std::vector<VariableValue *> _duh_dot_du;
};


template<>
InputParameters validParams<MultiKernel>();
