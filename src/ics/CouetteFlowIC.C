
#include "CouetteFlowIC.h"

template<>
InputParameters validParams<CouetteFlowIC>()
{
  InputParameters params = validParams<CFDInitialCondition>();
  return params;
}

CouetteFlowIC::CouetteFlowIC(const std::string & name, InputParameters parameters) :
    CFDInitialCondition(name, parameters)
{}

Real CouetteFlowIC::density(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real tem = temperature(p);
	Real pre = 1./(_gamma * _mach * _mach);
	return pre * _gamma * _mach * _mach/tem;

}

Real CouetteFlowIC::x_momentum(const Point &p)
{
	Real rho = density(p);
	Real u = p(1)/2.0;
	return rho * u;
}

Real CouetteFlowIC::y_momentum(const Point &p)
{
	Real rho = density(p);
	Real v = 0.;
	return rho * v;
}

Real CouetteFlowIC::z_momentum(const Point &p)
{
	return 0.0;
}

Real CouetteFlowIC::total_energy(const Point &p)
{
	Real rho = density(p);
	Real tem = temperature(p);
	Real pre = 1./(_gamma * _mach * _mach);
	Real u = p(1)/2.0;
	Real v = 0;
	Real w = 0;
	return pre/(_gamma-1)+0.5*rho * (u*u + v*v + w*w);
}

Real CouetteFlowIC::temperature(const Point& p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
	Real T2=0.85,T1=0.8;
	return T1 + ( T2 - T1 ) * y / 2 + 0.5 * _prandtl * (_gamma - 1) * _mach * _mach * y / 2 * ( 1 - y / 2 );
}
