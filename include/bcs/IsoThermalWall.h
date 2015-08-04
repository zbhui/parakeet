
#pragma once

#include "CFDBC.h"

class IsoThermalWall : public CFDBC
{
public:
	IsoThermalWall(const InputParameters & params);
	virtual ~IsoThermalWall(){}

protected:
	virtual void boundaryCondition();
};

template<>
InputParameters validParams<IsoThermalWall>();
