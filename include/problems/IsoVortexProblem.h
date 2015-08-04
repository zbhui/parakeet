
#pragma once

#include "EulerProblem.h"

class IsoVortexProblem : public EulerProblem
{
public:
	IsoVortexProblem(const InputParameters &params);
	virtual Real initialCondition(const Point & point, int eq);
	virtual Real boundaryCondition(Real t, const Point & point, int eq);
private:
	Real valueExact(Real t, const Point &p, int eq);

	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
};

template<>
InputParameters validParams<IsoVortexProblem>();
