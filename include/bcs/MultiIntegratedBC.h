
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

	  virtual void computeResidual();
	  virtual void precalculateResidual() = 0;
	  virtual Real computeQpResidual(unsigned int p) = 0;

	  virtual void computeJacobian();
	  virtual void precalculateJacobian() = 0;
	  virtual Real computeQpJacobian(unsigned int p, unsigned int q) = 0;

	  virtual Real computeQpResidual();
	  void computeJacobianBlock(unsigned int jvar);
//	  void computeJacobianBlockScalar(unsigned int jvar);

protected:
	  virtual void valueAtLeftFace(Real *ul);
	  virtual void valueAtRightFace(Real *ur);
	  virtual void valueGradAtLeftFace(RealGradient *dul);
	  virtual void valueGradAtRightFace(RealGradient *dur);

	  std::vector<NonlinearVariableName> _variables;
	  unsigned int _ep;
	  unsigned int _eq;
	  unsigned int _n_equation;
	  std::vector<VariableValue*> _uh;
	  std::vector<VariableGradient*> _grad_uh;
};

