
#pragma once

#include "MooseEnum.h"

#include "MultiIntegratedBC.h"
#include "CFDBase.h"

class CFDBC :
public MultiIntegratedBC
{
public:
	CFDBC(const InputParameters & params);
	virtual ~CFDBC(){}

protected:
	MooseEnum _bc_type;

};

template<>
InputParameters validParams<CFDBC>();
