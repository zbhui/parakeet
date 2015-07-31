
#pragma once

#include "IntegratedBC.h"

class MultiIntegratedBC;

template<>
InputParameters validParams<MultiIntegratedBC>();

class MultiIntegratedBC :
public IntegratedBC
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

	  std::vector<NonlinearVariableName> _variables;
	  unsigned int _n_equation;
	  std::vector<VariableValue*> _uh;
	  std::vector<VariableGradient*> _grad_uh;
};

