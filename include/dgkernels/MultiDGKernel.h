
#pragma once

#include "DGKernel.h"

class MultiDGKernel;


template<>
InputParameters validParams<MultiDGKernel>();

/**
 * \p MultiDGKernel 计算多个变量的单元积分
 */
class MultiDGKernel :
  public DGKernel
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

	unsigned int _n_equation;
	std::vector<NonlinearVariableName> _variables;
	std::vector<VariableValue*> _uh;
	std::vector<VariableValue*> _uh_neighbor;
	std::vector<VariableGradient*> _grad_uh;
	std::vector<VariableGradient*> _grad_uh_neighbor;
};

