
#include "Riemann2DProblem.h"
#include "boost/assign.hpp"
using namespace boost::assign;

template<>
InputParameters validParams<Riemann2DProblem>()
{
  InputParameters params = validParams<EulerProblem>();

  return params;
}

Riemann2DProblem::Riemann2DProblem(const InputParameters &params) :
	EulerProblem(params)
{
//	setShockDepart();
//	_shock_depart[0] += 0.5313,0,0,0,0.4;
//	_shock_depart[1] += 1,0.7276,0,0,1;
//	_shock_depart[2] += 0.8,0,0,0,1;
//	_shock_depart[3] += 1,0,0.7276,0,1;
}

Real Riemann2DProblem::density(Real t, const Point &p)
{
	return _shock_depart[pointLocator(p)][0];
}

Real Riemann2DProblem::momentumX(Real t, const Point &p)
{
	Real u = _shock_depart[pointLocator(p)][1];
	Real rho = density(t, p);

	return rho*u;
}

Real Riemann2DProblem::momentumY(Real t, const Point &p)
{
	Real u = _shock_depart[pointLocator(p)][2];
	Real rho = density(t, p);

	return rho*u;
}

Real Riemann2DProblem::momentumZ(Real t, const Point &p)
{
	Real u = _shock_depart[pointLocator(p)][3];
	Real rho = density(t, p);

	return rho*u;
}

Real Riemann2DProblem::energyTotal(Real t, const Point &p)
{
	RealVectorValue momentum(momentumX(t, p), momentumY(t, p), momentumZ(t, p));
	Real rho = density(t, p);
	Real pre = pressure(t, p);

	return pre/(_gamma-1) +0.5*momentum.size_sq()/rho;
}

Real Riemann2DProblem::pressure(Real t, const Point &p)
{
	return _shock_depart[pointLocator(p)][4];
}

//int Riemann2DProblem::pointLocator(const Point& p)
//{
//	Real x = p(0);
//	Real y = p(1);
//    if (x > 0.5 && y > 0.5)
//		return 0;
//    else if (x < 0.5 && y > 0.5)
//    	return 1;
//    else if (x < 0.5 && y < 0.5)
//    	return 2;
//    else
//    	return 3;
//}
