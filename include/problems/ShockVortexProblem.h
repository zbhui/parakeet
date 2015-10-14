
#pragma once

#include "Riemann2DProblem.h"

class ShockVortexProblem : public Riemann2DProblem
{
public:
	ShockVortexProblem(const InputParameters &params);

protected:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real pressure(Real t, const Point &p);

	int pointLocator( const Point &p);
	Point _p0;
	Real _e, _a, _rc;
};

template<>
InputParameters validParams<ShockVortexProblem>();
