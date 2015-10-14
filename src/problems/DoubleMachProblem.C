
#include "DoubleMachProblem.h"
#include "boost/assign.hpp"
using namespace boost::assign;

template<>
InputParameters validParams<DoubleMachProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

DoubleMachProblem::DoubleMachProblem(const InputParameters &params) :
	EulerProblem(params)
{
	_mach = 10;
	Real mach2 = _mach*_mach;
	Real r1 =1.4;
	Real p1 = 1;
	Real c = sqrt(_gamma*p1/r1);

	_initial_condition[0] += r1,0, 0, 0, p1;
	Real r2 = (_gamma+1)*mach2/((_gamma-1)*mach2+2)*r1;
	Real p2 = (2*_gamma/(_gamma+1)*mach2 - (_gamma-1)/(_gamma+1))*p1;
	Real u = 2*c/(_gamma+1)*(_mach-1./_mach);
	_initial_condition[1] += r2, u*cos(_attack), u*sin(_attack), 0, p2;
}

Real DoubleMachProblem::density(Real t, const Point &p)
{
	Point center(0.25, 0.5);

	return _initial_condition[pointLocator(p)][0];
}

Real DoubleMachProblem::momentumX(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][1];
	Real rho = density(t, p);

	return rho*u;
}

Real DoubleMachProblem::momentumY(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][2];
	Real rho = density(t, p);

	return rho*u;
}

Real DoubleMachProblem::momentumZ(Real t, const Point &p)
{
	Real u = _initial_condition[pointLocator(p)][3];
	Real rho = density(t, p);

	return rho*u;
}

Real DoubleMachProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real DoubleMachProblem::pressure(Real t, const Point &p)
{
	return _initial_condition[pointLocator(p)][4];
}

int DoubleMachProblem::pointLocator(const Point& p)
{

	Real r1, p1, r2, p2;
	r1 = _initial_condition[0][0];
	r2 = _initial_condition[1][0];

	p1 = _initial_condition[0][4];
	p2 = _initial_condition[1][4];

	Real vs = 10;

	Point normal(cos(_attack), sin(_attack));
	Point p0(1./6, 0);
	p0 += vs*_time*normal;
	bool outer = (p-p0)*normal > 0;

    if ((p-p0)*normal > 0)
		return 0;
    else
    	return 1;

}
