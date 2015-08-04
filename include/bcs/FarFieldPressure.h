
#pragma once

#include "CFDBC.h"

class FarFieldPressure : public CFDBC
{
public:
	FarFieldPressure(const InputParameters & params);
	virtual ~FarFieldPressure(){}

protected:
	virtual void boundaryCondition();

	Real _gamma;
	Real _mach;
	Real _rho_inf;
	Real _vel_inf;
};

template<>
InputParameters validParams<FarFieldPressure>();
