
#include "ShockVortexProblem.h"
#include "boost/assign.hpp"
using namespace boost::assign;

template<>
InputParameters validParams<ShockVortexProblem>()
{
  InputParameters params = validParams<Riemann2DProblem>();

  return params;
}

ShockVortexProblem::ShockVortexProblem(const InputParameters &params) :
	Riemann2DProblem(params),
	_p0(0.25, 0.5, 0),
	_e(0.3),
	_a(0.204),
	_rc(0.05)
{
	_mach = 1.1;
	_shock_depart[0] += 1,1.1*sqrt(_gamma), 0, 0, 1;
	Real mach2 = _mach*_mach;
	Real rho = (_gamma+1)*mach2/((_gamma-1)*mach2+2);
	Real p = 2*_gamma/(_gamma+1)*mach2 - (_gamma-1)/(_gamma+1);
	Real u = 1.1*sqrt(_gamma)/rho;
	_shock_depart[1] += rho, u, 0, 0, p;
}

Real ShockVortexProblem::density(Real t, const Point &p)
{
	Point normal(p - _p0);
	Real r = (p-_p0).size();
	Real tau = r/_rc;
	Real du =  _e/_rc*exp(_a*(1-tau*tau))*normal(1);
	Real dv = -_e/_rc*exp(_a*(1-tau*tau))*normal(0);
	Real dT = -(_gamma-1)/4./_a/_gamma*_e*_e*exp(2*_a*(1-tau*tau));

	Real u = _shock_depart[pointLocator(p)][1]+du;
	Real T = pressure(t, p )/_shock_depart[pointLocator(p)][0] + dT;
	Real rho = pow( T, 1/(_gamma-1));

	return rho;
}

Real ShockVortexProblem::momentumX(Real t, const Point &p)
{
	Point normal(p - _p0);
	Real r = (p-_p0).size();
	Real tau = r/_rc;
	Real du = _e/_rc*exp(_a*(1-tau*tau))*normal(1);
	Real dv = -_e/_rc*exp(_a*(1-tau*tau))*normal(0);
	Real dT = -(_gamma-1)/4./_a/_gamma*_e*_e*exp(2*_a*(1-tau*tau));

	Real u = _shock_depart[pointLocator(p)][1]+du;
	Real T = pressure(t, p )/_shock_depart[pointLocator(p)][0] + dT;
	Real rho = pow( T, 1/(_gamma-1));

	return rho*u;
}

Real ShockVortexProblem::momentumY(Real t, const Point &p)
{
	Point normal(p - _p0);
	Real r = (p-_p0).size();
	Real tau = r/_rc;
	Real du =  _e/_rc*exp(_a*(1-tau*tau))*normal(1);
	Real dv = -_e/_rc*exp(_a*(1-tau*tau))*normal(0);
	Real dT = -(_gamma-1)/4./_a/_gamma*_e*_e*exp(2*_a*(1-tau*tau));

	Real u = _shock_depart[pointLocator(p)][2]+dv;
	Real T = pressure(t, p )/_shock_depart[pointLocator(p)][0] + dT;
	Real rho = pow( T, 1/(_gamma-1));

	return rho*u;
}

Real ShockVortexProblem::momentumZ(Real t, const Point &p)
{
	Point normal(p - _p0);
	Real r = (p-_p0).size();
	Real tau = r/_rc;
	Real du =  _e/_rc*exp(_a*(1-tau*tau))*normal(1);
	Real dv = -_e/_rc*exp(_a*(1-tau*tau))*normal(0);
	Real dT = -(_gamma-1)/4./_a/_gamma*_e*_e*exp(2*_a*(1-tau*tau));

	Real u = _shock_depart[pointLocator(p)][3];
	Real T = pressure(t, p )/_shock_depart[pointLocator(p)][0] + dT;
	Real rho = pow( T, 1/(_gamma-1));

	return rho*u;
}

Real ShockVortexProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));

	Point normal(p - _p0);
	Real r = (p-_p0).size();
	Real tau = r/_rc;
	Real du =  _e/_rc*exp(_a*(1-tau*tau))*normal(1);
	Real dv = -_e/_rc*exp(_a*(1-tau*tau))*normal(0);
	Real dT = -(_gamma-1)/4./_a/_gamma*_e*_e*exp(2*_a*(1-tau*tau));

	Real u = _shock_depart[pointLocator(p)][3];
	Real T = pressure(t, p )/_shock_depart[pointLocator(p)][0] + dT;
	Real rho = pow( T, 1/(_gamma-1));

	return rho*T/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real ShockVortexProblem::pressure(Real t, const Point &p)
{
	return _shock_depart[pointLocator(p)][4];
}

int ShockVortexProblem::pointLocator(const Point& p)
{
	Real x = p(0);
	Real y = p(1);
    if (x > 0.5 )
		return 0;
    else
    	return 1;
}
