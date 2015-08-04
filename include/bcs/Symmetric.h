
#pragma once

#include "CFDBC.h"

class Symmetric : public CFDBC
{
public:
	Symmetric(const InputParameters & params);
	virtual ~Symmetric(){}


protected:
	virtual void boundaryCondition();
};

template<>
InputParameters validParams<Symmetric>();
