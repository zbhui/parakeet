
#include "SodProblem.h"

template<>
InputParameters validParams<SodProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

SodProblem::SodProblem(const InputParameters &params) :
	EulerProblem(params)
{
}

Real SodProblem::density(Real t, const Point &p)
{
	Real x = p(0);
    if (x < 0.5)
		return 1;
    else
    	return 0.125;
}

Real SodProblem::momentumX(Real t, const Point &p)
{
	Real x = p(0);
	Real u = 0;
	Real rho = density(t, p);
    if (x < 0.5)
		u = 0;
    else
    	u = 0;

	return rho*u;
}

Real SodProblem::momentumY(Real t, const Point &p)
{
	return 0.0;
}

Real SodProblem::momentumZ(Real t, const Point &p)
{
	return 0.0;
}

Real SodProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real SodProblem::pressure(Real t, const Point &p)
{
	Real x = p(0);
    if (x < 0.5)
		return 1;
    else
    	return 0.1;
}

Real SodProblem::valueExact(Real t, const Point& p, int eq)
{
	mooseError("Sod问题的精确解是什么？？" << eq);
}
