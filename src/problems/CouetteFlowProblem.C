
#include "CouetteFlowProblem.h"

template<>
InputParameters validParams<CouetteFlowProblem>()
{
  InputParameters params = validParams<NavierStokesProblem>();

  return params;
}

CouetteFlowProblem::CouetteFlowProblem(const InputParameters &params) :
	NavierStokesProblem(params)
{
}

Real CouetteFlowProblem::density(Real t, const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real tem = temperature(t, p);
	Real pre = 1./(_gamma * _mach * _mach);
	return pre * _gamma * _mach * _mach/tem;

}

Real CouetteFlowProblem::momentumX(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real u = p(1)/2.0;
	return rho * u;
}

Real CouetteFlowProblem::momentumY(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real v = 0.;
	return rho * v;
}

Real CouetteFlowProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real CouetteFlowProblem::energyTotal(Real t, const Point &p)
{
	Real rho = density(t, p);
	Real tem = temperature(t, p);
	Real pre = 1./(_gamma * _mach * _mach);
	Real u = p(1)/2.0;
	Real v = 0;
	Real w = 0;
	return pre/(_gamma-1)+0.5*rho * (u*u + v*v + w*w);
}

Real CouetteFlowProblem::temperature(Real t, const Point& p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
	Real T2=0.85,T1=0.8;
	return T1 + ( T2 - T1 ) * y / 2 + 0.5 * _prandtl * (_gamma - 1) * _mach * _mach * y / 2 * ( 1 - y / 2 );
}

