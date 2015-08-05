
#pragma once

#include "CFDBC.h"

class FarFieldRiemann : public CFDBC
{
public:
	FarFieldRiemann(const InputParameters & params);
	virtual ~FarFieldRiemann(){}

protected:
	Real _rho_inf;
	Real _vel_inf;
	Real _tem_inf;

	virtual void boundaryCondition();
};

template<>
InputParameters validParams<FarFieldRiemann>();
