
#pragma once

#include "Riemann2DProblem.h"

class ForwardStepProblem : public Riemann2DProblem
{
public:
	ForwardStepProblem(const InputParameters &params);

protected:
	virtual int pointLocator( const Point &p);
};

template<>
InputParameters validParams<ForwardStepProblem>();
