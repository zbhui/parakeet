
#pragma once

#include "CFDBC.h"

class SlipWall : public CFDBC
{
public:
	SlipWall(const InputParameters & params);
	virtual ~SlipWall(){}


protected:
	virtual void boundaryCondition();
};

template<>
InputParameters validParams<SlipWall>();
