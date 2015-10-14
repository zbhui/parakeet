
#include "ShockVortexProblem.h"
#include "boost/assign.hpp"
using namespace boost::assign;

template<>
InputParameters validParams<ShockVortexProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

ShockVortexProblem::ShockVortexProblem(const InputParameters &params) :
	EulerProblem(params)
{
	_mach = 1.1;
	_initial_condition[0] += 1,1.1*sqrt(_gamma), 0, 0, 1;
	Real mach2 = _mach*_mach;
	Real rho = (_gamma+1)*mach2/((_gamma-1)*mach2+2);
	Real p = 2*_gamma/(_gamma+1)*mach2 - (_gamma-1)/(_gamma+1);
	Real u = 1.1*sqrt(_gamma)/rho;
	_initial_condition[1] += rho, u, 0, 0, p;
}

Real ShockVortexProblem::density(Real t, const Point &p)
{
	Point center(0.25, 0.5);

	return _initial_condition[pointLocator(p)][0];
}

Real ShockVortexProblem::momentumX(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][1];
	Real rho = density(t, p);

	return rho*u;
}

Real ShockVortexProblem::momentumY(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][2];
	Real rho = density(t, p);

	return rho*u;
}

Real ShockVortexProblem::momentumZ(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][3];
	Real rho = density(t, p);

	return rho*u;
}

Real ShockVortexProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real ShockVortexProblem::pressure(Real t, const Point &p)
{
	return _initial_condition[pointLocator(p)][4];
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
