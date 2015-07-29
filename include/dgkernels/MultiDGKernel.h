
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
  virtual void  precalculateResidual() = 0;
  virtual Real computeQpResidual(Moose::DGResidualType type, unsigned int p) = 0;
  virtual Real computeQpResidual(Moose::DGResidualType type);

  virtual void precalculateJacobian() = 0;
  virtual void computeJacobian();
  virtual Real computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q) = 0;
  virtual void computeOffDiagJacobian(unsigned int jvar);

  virtual void computeElemNeighResidual(Moose::DGResidualType type);
  virtual void computeElemNeighJacobian(Moose::DGJacobianType type);

protected:

  virtual void valueAtLeftFace(Real *ul);
  virtual void valueAtRightFace(Real *ur);
  virtual void valueGradAtLeftFace(RealGradient *dul);
  virtual void valueGradAtRightFace(RealGradient *dur);

  virtual Real computeQpJacobian(Moose::DGJacobianType type);

  unsigned int _n_equation;

  /// 当前kernel的变量名
  std::vector<NonlinearVariableName> _variables;

  /// 积分点上的左变量值
  std::vector<VariableValue*> _uh;

  /// 积分点上的右变量值
  std::vector<VariableValue*> _uh_neighbor;

  /// 积分点上的左变量梯度
  std::vector<VariableGradient*> _grad_uh;

  /// 积分点上的右变量梯度
  std::vector<VariableGradient*> _grad_uh_neighbor;
};

