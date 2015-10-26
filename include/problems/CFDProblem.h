
#pragma once

#include "FEProblem.h"
#include "Attitude.h"
#include "Eigen/Geometry"
using Eigen::Quaterniond;

class CFDProblem : public FEProblem
{
public:
	CFDProblem(const InputParameters & params);
	~CFDProblem(){}

	virtual Real initialCondition(const Point & point, int eq);
	virtual Real boundaryCondition(Real t, const Point & point, int eq);
	virtual Real valueExact(Real t, const Point &p, int eq);
	static MooseEnum viscousType();
	static MooseEnum fluxRiemannType();

private:
	typedef FEProblem Parent;
	virtual Real density(Real t, const Point &p);
	virtual Real momentumX(Real t, const Point &p);
	virtual Real momentumY(Real t, const Point &p);
	virtual Real momentumZ(Real t, const Point &p);
	virtual Real energyTotal(Real t, const Point &p);
	virtual void computeJacobian(NonlinearImplicitSystem & sys, const NumericVector<Number> & soln, SparseMatrix<Number> &  jacobian);
	virtual void initialSetup();
public:
	int _var_order;
	MooseEnum _vis_type;
	MooseEnum _flux_type;
	Real _mach;
	Real _gamma;
	Real _reynolds;
	Real _prandtl;

	Real _attack;			/// 攻角
	Real _sideslip;		///侧滑角
	Real _pitch;			///俯仰角
	Real _yaw;			///偏航角
	Real _roll;			///滚转角

	Attitude _attitude;
	Real _velocity;

	unsigned int _jacobian_delay;
};

template<>
InputParameters validParams<CFDProblem>();
