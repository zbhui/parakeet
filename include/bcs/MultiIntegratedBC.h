
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
	  virtual void computeJacobian();
	  /**
	   * Computes d-ivar-residual / d-jvar...
	   */
	  void computeJacobianBlock(unsigned int jvar);
	  /**
	   * Computes jacobian block with respect to a scalar variable
	   * @param jvar The number of the scalar variable
	   */
//	  void computeJacobianBlockScalar(unsigned int jvar);

protected:
	  virtual void valueAtLeftFace(Real *ul);
	  virtual void valueAtRightFace(Real *ur);
	  virtual void valueGradAtLeftFace(RealGradient *dul);
	  virtual void valueGradAtRightFace(RealGradient *dur);

	  virtual void precalculateResidual() = 0;
	  virtual void precalculateJacobian(){}

	  /// 当前kernel的变量名
	  std::vector<NonlinearVariableName> _variables;
	  /// 方程索引
	  unsigned int _ep;
	  unsigned int _eq;

	  /// 方程个数
	  unsigned int _n_equation;

	  /// 积分点上的左变量值
	  std::vector<VariableValue*> _uh;

	  /// 积分点上的左变量梯度
	  std::vector<VariableGradient*> _grad_uh;
};

