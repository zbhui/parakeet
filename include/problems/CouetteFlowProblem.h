
#pragma once

#include "NavierStokesProblem.h"

class CouetteFlowProblem : public NavierStokesProblem
{
public:
	CouetteFlowProblem(const InputParameters &params);

	virtual Real initialCondition(const Point & point, int eq);
	Real valueExact(Real t, const Point &p, int eq);
private:
	Real density(Real t, const Point &p);
	Real momentumX(Real t, const Point &p);
	Real momentumY(Real t, const Point &p);
	Real momentumZ(Real t, const Point &p);
	Real energyTotal(Real t, const Point &p);
	Real temperature(Real t, const Point &p);

};

template<>
InputParameters validParams<CouetteFlowProblem>();
