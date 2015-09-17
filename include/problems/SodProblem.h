
#pragma once

#include "EulerProblem.h"

class SodProblem : public EulerProblem
{
public:
	SodProblem(const InputParameters &params);

	Real valueExact(Real t, const Point& p, int eq);
private:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real pressure(Real t, const Point &p);
};

template<>
InputParameters validParams<SodProblem>();
