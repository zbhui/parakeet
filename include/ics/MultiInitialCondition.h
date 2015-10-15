
#pragma once

#include "InitialCondition.h"


class MultiInitialCondition : public InitialCondition
{
public:
	MultiInitialCondition(const InputParameters & parameters);
	virtual Real value(const Point & p);
	void compute();
protected:
	virtual Real value(int component, const Point & p);

	bool _constant_ic;
};

template<>
InputParameters validParams<MultiInitialCondition>();
