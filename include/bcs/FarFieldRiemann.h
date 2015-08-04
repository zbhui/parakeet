
#pragma once

#include "CFDBC.h"

class FarFieldRiemann : public CFDBC
{
public:
	FarFieldRiemann(const InputParameters & params);
	virtual ~FarFieldRiemann(){}


protected:
	virtual void boundaryCondition();
};

template<>
InputParameters validParams<FarFieldRiemann>();
