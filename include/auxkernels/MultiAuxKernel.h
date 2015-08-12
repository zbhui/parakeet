
#pragma once

#include "AuxKernel.h"

class MultiAuxKernel : public AuxKernel
{
public:
	MultiAuxKernel(const InputParameters & parameters);
	virtual ~MultiAuxKernel(){};

	virtual void compute();

protected:
	std::vector<AuxVariableName> _aux_variables;
	int _ivar;
};


template<>
InputParameters validParams<MultiAuxKernel>();
