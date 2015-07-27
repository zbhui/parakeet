
#pragma once

#include "Kernel.h"

class MultiKernel;


template<>
InputParameters validParams<MultiKernel>();

/**
 * \p MultiKernel 计算多个变量的单元积分
 */
class MultiKernel :
  public Kernel
{
public:
	MultiKernel(const InputParameters & parameters);

  virtual ~MultiKernel(){};

  /// 当前kernel对residual的贡献
  virtual void computeResidual();
  virtual Real computeQpResidual();
  virtual Real computeQpResidual(unsigned int p) = 0;
  virtual void computeJacobian();

  virtual Real computeQpJacobian(unsigned int p, unsigned int q) = 0;
  virtual void computeOffDiagJacobian(unsigned int jvar);
  virtual void computeOffDiagJacobianScalar(unsigned int jvar);

protected:

  /// 组合积分点上的变量
  virtual void valueAtCellPoint(Real *uh);

  virtual void valueGradAtCellPoint(RealGradient *duh);

  /// 在\p computeReisudal 之前的数据准备
  virtual void precalculateResidual() = 0;
  /// 在\p computeJacobian 之前的数据准备
  virtual void precalculateJacobian();

  /// Compute this Kernel's contribution to the Jacobian at the current quadrature point
//  virtual Real computeQpJacobian();

  /// This is the virtual that derived classes should override for computing an off-diagonal Jacobian component.
//  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// This callback is used for Kernels that need to perturb residual calculations

  /// 方程索引
  unsigned int _eq;
  unsigned int _ep;

  /// 方程个数
  unsigned int _n_equation;

  /// 当前kernel的变量名
  std::vector<NonlinearVariableName> _variables;

  /// 积分点上的变量值
  std::vector<VariableValue*> _uh;

  /// 积分点上的变量梯度
  std::vector<VariableGradient*> _grad_uh;

  /// 变量的时间导数
  std::vector<VariableValue *> _uh_dot;

  /// Derivative of u_dot with respect to u
  std::vector<VariableValue *> _duh_dot_du;
};

