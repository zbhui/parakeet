
#pragma once
#include "CFDBC.h"
#include "CFDDataPack.h"

class EulerBC :
		public CFDBC
{
public:
	EulerBC(const InputParameters & params);
	virtual ~EulerBC(){}


protected:
	MooseEnum _bc_type;

	virtual void boundaryCondition();

	virtual void wall();
	virtual void farfield();
	virtual void symmetric();
};

template<>
InputParameters validParams<EulerBC>();
