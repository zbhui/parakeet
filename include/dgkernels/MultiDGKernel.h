
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
	MultiDGKernel(const std::string & name, InputParameters parameters);

  virtual ~MultiDGKernel(){};

  /// 当前side对residual的贡献
  virtual void computeResidual();
  virtual void computeJacobian();
  virtual void computeOffDiagJacobian(unsigned int jvar);

   ///Computes the element/neighbor-element/neighbor Jacobian
  virtual void computeElemNeighResidual(Moose::DGResidualType type);
  virtual void computeElemNeighJacobian(Moose::DGJacobianType type);
  virtual void computeOffDiagElemNeighJacobian(Moose::DGJacobianType type,unsigned int jvar);

protected:

  virtual void valueAtLeftFace(Real *ul);
  virtual void valueAtRightFace(Real *ur);
  virtual void valueGradAtLeftFace(RealGradient *dul);
  virtual void valueGradAtRightFace(RealGradient *dur);

  /// Compute this Kernel's contribution to the Jacobian at the current quadrature point
  virtual Real computeQpJacobian(Moose::DGJacobianType type);
//  virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

  virtual void  precalculateResidual() = 0;
  virtual void precalculateJacobian(){};

  /// 方程索引
  unsigned int _ep;
  unsigned int _eq;

  /// 方程个数
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

