
#pragma once

#include "InitialCondition.h"
#include "KOmegaModelBase.h"

class KOmegaIC;

template<>
InputParameters validParams<KOmegaIC>();

/**
 * KOmegaIC问题初始条件
 */
class KOmegaIC :
public InitialCondition,
public KOmegaModelBase
{
public:

	KOmegaIC(const InputParameters & parameters);
	~KOmegaIC(){}
protected:
	virtual Real value(const Point & p);

	int _component;
};
