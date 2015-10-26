
#include "CFDProblem.h"

using namespace Eigen;

template<>
InputParameters validParams<CFDProblem>()
{
  InputParameters params = validParams<FEProblem>();
  MooseEnum vis_type(CFDProblem::viscousType());
  params.addParam<MooseEnum>("vis_type", vis_type, "粘性计算方法");
  MooseEnum flux_type(CFDProblem::fluxRiemannType());
  params.addParam<MooseEnum>("flux_type", flux_type, "Riemann通量类型");
  params.addParam<Real>("mach",  0.1, "马赫数");
  params.addParam<Real>("gamma", 1.4, "比热比");
  params.addParam<Real>("reynolds", 1, "雷诺数");
  params.addParam<Real>("prandtl", 0.72, "prandtl数");
  params.addParam<Real>("attack", 0., "攻角");
  params.addParam<Real>("sideslip", 0., "侧滑角");
  params.addParam<Real>("pitch", 0., "俯仰角");
  params.addParam<Real>("yaw", 180., "偏航角");
  params.addParam<Real>("roll", -90., "滚转角");
  params.addParam<Real>("init_vel", 1.0, "初始流体速度");
  params.addParam<unsigned int>("jacobian_delay", 1, "jacobian矩阵更新频率");
  return params;
}

CFDProblem::CFDProblem(const InputParameters &params) :
	FEProblem(params),
	_var_order(1),
	_vis_type(getParam<MooseEnum>("vis_type")),
	_flux_type(getParam<MooseEnum>("flux_type")),
	_mach(getParam<Real>("mach")),
	_gamma(getParam<Real>("gamma")),
	_reynolds(getParam<Real>("reynolds")),
	_prandtl(getParam<Real>("prandtl")),

	_attack(getParam<Real>("attack")*libMesh::pi/180),
	_sideslip(getParam<Real>("sideslip")*libMesh::pi/180),
	_pitch(getParam<Real>("pitch")*libMesh::pi/180),
	_yaw(getParam<Real>("yaw")*libMesh::pi/180),
	_roll(getParam<Real>("roll")*libMesh::pi/180),

	_attitude(_attack, _sideslip, _pitch, _yaw, _roll),
    _velocity(getParam<Real>("init_vel")),
	_jacobian_delay(getParam<unsigned int>("jacobian_delay"))
{
}

void CFDProblem::initialSetup()
{
	FEProblem::initialSetup();
	computeIndicatorsAndMarkers();
}

void CFDProblem::computeJacobian(NonlinearImplicitSystem & sys, const NumericVector<Number> & soln, SparseMatrix<Number> &  jacobian)
{
	if((_t_step-1) % _jacobian_delay == 0 || _app.isRecovering() || _app.isRestarting() || !converged())
		Parent::computeJacobian(sys, soln, jacobian);
//	 if (!_has_jacobian || !_const_jacobian)
//	  {
//	    _nl.setSolution(soln);
//
//	    _nl.zeroVariablesForJacobian();
//	    _aux.zeroVariablesForJacobian();
//
//	    unsigned int n_threads = libMesh::n_threads();
//
//	    // Random interface objects
////	    for (std::map<std::string, RandomData *>::iterator it = _random_data_objects.begin();
////	         it != _random_data_objects.end();
////	         ++it)
////	      it->second->updateSeeds(EXEC_NONLINEAR);
//
//	    execTransfers(EXEC_NONLINEAR);
//	    execMultiApps(EXEC_NONLINEAR);
//
//	    computeUserObjects(EXEC_NONLINEAR, UserObjectWarehouse::PRE_AUX);
//
//	    if (_displaced_problem != NULL)
//	      _displaced_problem->updateMesh(soln, *_aux.currentSolution());
//
//	    for (unsigned int i=0; i<n_threads; i++)
//	    {
//	      _materials[i].jacobianSetup();
//
//	      for (std::map<std::string, MooseSharedPointer<Function> >::iterator vit = _functions[i].begin();
//	          vit != _functions[i].end();
//	          ++vit)
//	        vit->second->jacobianSetup();
//	    }
//
//	    _aux.jacobianSetup();
//
//	    _aux.compute(EXEC_NONLINEAR);
//
//	    computeUserObjects(EXEC_NONLINEAR, UserObjectWarehouse::POST_AUX);
//
////	    _app.getOutputWarehouse().jacobianSetup();
//
////	    if(_t_step % 1 == 1)
//	    {
//	    std::cout << _t_step <<std::endl;
//	    _nl.computeJacobian(jacobian);
//
//	    _has_jacobian = true;
//	    }
//	  }
//
//	  if (_solver_params._type == Moose::ST_JFNK || _solver_params._type == Moose::ST_PJFNK)
//	  {
//	    // This call is here to make sure the residual vector is up to date with any decisions that have been made in
//	    // the Jacobian evaluation.  That is important in JFNK because that residual is used for finite differencing
//	    computeResidual(sys, soln, *sys.rhs);
//	    sys.rhs->close();
//	  }
}


MooseEnum CFDProblem::viscousType()
{
  return MooseEnum("INVISCOUS CONSTANT  SUTHERLAND", "CONSTANT");
}

MooseEnum CFDProblem::fluxRiemannType()
{
  return MooseEnum("Lax-F HLL  HLLC-PV HLLC-Roe", "Lax-F");
}

Real CFDProblem::valueExact(Real t, const Point &p, int eq)
{
	switch (eq) {
		case 0:
			return density(t, p);
			break;
		case 1:
			return momentumX(t, p);
			break;
		case 2:
			return momentumY(t, p);
			break;
		case 3:
			return momentumZ(t, p);
			break;
		case 4:
			return energyTotal(t, p);
			break;
		default:
			return 0.0;
			break;
	}
}

Real CFDProblem::initialCondition(const Point& p, int eq)
{
	Real t = 0;
	switch (eq) {
			case 0:
				return density(t, p);
				break;
			case 1:
				return momentumX(t, p);
				break;
			case 2:
				return momentumY(t, p);
				break;
			case 3:
				return momentumZ(t, p);
				break;
			case 4:
				return energyTotal(t, p);
				break;
			default:
				return 0.0;
				break;
		}
}

Real CFDProblem::boundaryCondition(Real t, const Point& p, int eq)
{
	switch (eq) {
			case 0:
				return density(t, p);
				break;
			case 1:
				return momentumX(t, p);
				break;
			case 2:
				return momentumY(t, p);
				break;
			case 3:
				return momentumZ(t, p);
				break;
			case 4:
				return energyTotal(t, p);
				break;
			default:
				return 0.0;
				break;
		}
}

Real CFDProblem::density(Real t, const Point &p)
{
	return 1.0;
}

Real CFDProblem::momentumX(Real t, const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	return density(t, p)*vel(0);
}

Real CFDProblem::momentumY(Real t, const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	return density(t, p)*vel(1);
}

Real CFDProblem::momentumZ(Real t, const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	if(_mesh.dimension() == 2)
		return 0.;
	else if(_mesh.dimension() == 3)
		return density(t, p)*vel(2);
	else
	{
		mooseError("一维问题此处需要调试");
		return 0.;
	}
}

Real CFDProblem::energyTotal(Real t, const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(t, p)*(_velocity*_velocity);
}
