
#include "CFDProblem.h"

#include "MooseApp.h"
using namespace Eigen;

template<>
InputParameters validParams<CFDProblem>()
{
  InputParameters params = validParams<FEProblem>();
  MooseEnum vis_type(CFDProblem::getViscousType());
  params.addParam<MooseEnum>("vis_type", vis_type, "粘性计算方法");
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
  return params;
}

CFDProblem::CFDProblem(const InputParameters &params) :
	FEProblem(params),
	_vis_type(getParam<MooseEnum>("vis_type")),
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
    _velocity(getParam<Real>("init_vel"))
{
}

MooseEnum CFDProblem::getViscousType()
{
  return MooseEnum("INVISCOUS CONSTANT  SUTHERLAND", "CONSTANT");
}
Real CFDProblem::initialCondition(const Point& p, int eq)
{
	switch (eq) {
		case 0:
			return density(p);
			break;
		case 1:
			return momentumX(p);
			break;
		case 2:
			return momentumY(p);
			break;
		case 3:
			return momentumZ(p);
			break;
		case 4:
			return energyTotal(p);
			break;
		default:
			return 0.0;
			break;
	}
}

Real CFDProblem::density(const Point &p)
{
	return 1.0;
}

Real CFDProblem::momentumX(const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(0);
}

Real CFDProblem::momentumY(const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	return density(p)*vel(1);
}

Real CFDProblem::momentumZ(const Point &p)
{
	Vector3d vel = _velocity*(_attitude.earthFromWind()*Vector3d::UnitX());
	if(_mesh.dimension() == 2)
		return 0.;
	else if(_mesh.dimension() == 3)
		return density(p)*vel(2);
	else
	{
		mooseError("一维问题此处需要调试");
		return 0.;
	}
}

Real CFDProblem::energyTotal(const Point &p)
{
	Real pre = 1./_gamma/_mach/_mach;
	return pre/(_gamma-1) + 0.5*density(p)*(_velocity*_velocity);
}
