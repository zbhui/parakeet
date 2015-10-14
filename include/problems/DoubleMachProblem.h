
#pragma once

#include "Riemann2DProblem.h"

class DoubleMachProblem : public Riemann2DProblem
{
public:
	DoubleMachProblem(const InputParameters &params);

protected:
	int pointLocator( const Point &p);
};

template<>
InputParameters validParams<DoubleMachProblem>();
