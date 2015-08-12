
#pragma once

#include "MultiAuxKernel.h"
#include "CFDDataPack.h"

class CFDProblem;

class NSAuxVariable : public MultiAuxKernel
{
public:
  NSAuxVariable(const InputParameters & parameters);

protected:
    virtual Real computeValue();

protected:
	FEProblem & _fe_problem;
	CFDProblem &_cfd_problem;
	CFDDataPack _cfd_data;
	std::vector<VariableValue*> _uh;
};

template<>
InputParameters validParams<NSAuxVariable>();
