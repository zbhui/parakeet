
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
	virtual Real boundaryCondition(Real t, const Point & point, int eq){ return 0;}
	static MooseEnum getViscousType();
private:
	virtual Real density(const Point &p);
	virtual Real momentumX(const Point &p);
	virtual Real momentumY(const Point &p);
	virtual Real momentumZ(const Point &p);
	virtual Real energyTotal(const Point &p);
	virtual void computeJacobian(NonlinearImplicitSystem & sys, const NumericVector<Number> & soln, SparseMatrix<Number> &  jacobian);

public:
	int _var_order;
	MooseEnum _vis_type;
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
};

template<>
InputParameters validParams<CFDProblem>();
