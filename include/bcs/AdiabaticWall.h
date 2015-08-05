
#pragma once

#include "CFDBC.h"

class AdiabaticWall : public CFDBC
{
public:
	AdiabaticWall(const InputParameters & params);
	virtual ~AdiabaticWall(){}

protected:
	virtual void boundaryCondition();
};

template<>
InputParameters validParams<AdiabaticWall>();
