
#pragma once

#include "InitialCondition.h"

class CFDProblem;

class CFDInitialCondition : public InitialCondition
{
public:
	CFDInitialCondition(const InputParameters & parameters);
	virtual Real value(const Point & p);
	void compute();
protected:
	CFDProblem &_cfd_problem;
	unsigned int _component;
	virtual Real value(int component, const Point & p);

	bool _constant_ic;
};

template<>
InputParameters validParams<CFDInitialCondition>();
