
#pragma once

#include "EulerProblem.h"

class Riemann2DProblem : public EulerProblem
{
public:
	Riemann2DProblem(const InputParameters &params);
private:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real pressure(Real t, const Point &p);

	std::vector<Real> _initial_condition[4];
	int pointLocator( const Point &p);
};

template<>
InputParameters validParams<Riemann2DProblem>();
