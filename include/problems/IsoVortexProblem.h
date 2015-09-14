
#pragma once

#include "EulerProblem.h"

class IsoVortexProblem : public EulerProblem
{
public:
	IsoVortexProblem(const InputParameters &params);

private:
 	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
};

template<>
InputParameters validParams<IsoVortexProblem>();
