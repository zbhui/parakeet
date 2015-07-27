
#pragma once

#include "MooseVariableBase.h"
#include "MooseObject.h"
#include "Eigen/Dense"
typedef Eigen::Matrix<Real, 5, 5> Matrix5x5;

class CFDBase;

template<>
InputParameters validParams<CFDBase>();

/**
 * CFD base class.
 */
class CFDBase
{
public:
	CFDBase(const std::string & name, InputParameters parameters);
  virtual ~CFDBase(){};

protected:
//  virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
//  virtual void fluxViscous(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur, Point &normal);

  Real pressure(Real *uh);
  Real enthalpy(Real *uh);
  Real temperature(Real  *uh);
  Real mach_local(Real *uh);
  Real physicalViscosity(Real *uh);


  virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
  virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);

  virtual void inviscousJacobian(Matrix5x5 *inviscous_jac, Real *uh);
  virtual void inviscousJacobianFD(Matrix5x5 *inviscous_jac, Real *uh);
protected:

//  int _n_equation;

/// Required parameters
  Real _gamma;
  Real _prandtl;
  Real _reynolds;
  Real _mach;

  Real _attack;
  Real _slide;

  Real _epsilon;
  Real _sigma;

  /// 有限差分求Jacobian矩阵的变量增量
  Real _ds;
};

