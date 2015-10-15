
#pragma once

#include "MooseVariable.h"
#include "MooseVariableScalar.h"
#include "InputParameters.h"

class MultiVariableInterface
{
public:
	MultiVariableInterface(const InputParameters & parameters);

	virtual ~MultiVariableInterface(){};

protected:
	FEProblem & _mv_fe_problem;
	SystemBase & _sys;
	THREAD_ID _tid;

	unsigned int _n_equation;
	std::vector<NonlinearVariableName> _variables;
	std::vector<VariableValue*> _uh;
	std::vector<VariableValue*> _uh_neighbor;
	std::vector<VariableGradient*> _grad_uh;
	std::vector<VariableGradient*> _grad_uh_neighbor;
	int _var_order;
};

template<>
InputParameters validParams<MultiVariableInterface>();
