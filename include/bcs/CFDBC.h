
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

	 virtual void valueAtRightFace(Real *ur);

	 virtual void wallBC(Real *ur) = 0;
	 virtual void farFieldBC(Real *ur) = 0;
	 virtual void symmetricBC(Real *ur) = 0;

private:
	MooseEnum _bc_type;

};

template<>
InputParameters validParams<CFDBC>();
