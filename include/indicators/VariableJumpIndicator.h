
#pragma once

#include "InternalSideIndicator.h"
#include <vector>
using std::vector;

class VariableJumpIndicator : public InternalSideIndicator
{
public:
	VariableJumpIndicator(const InputParameters &parameters);
	virtual ~VariableJumpIndicator(){};

protected:
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	vector<VariableName> _aux_variables;
	int _n_variables;
	int _var_order;

	vector<VariableValue*> _uh;
	vector<VariableValue*> _uh_neighbor;
	vector<VariableGradient*> _grad_uh;
	vector<VariableGradient*> _grad_uh_neighbor;

	bool _is_implicit;
	virtual Real computeQpIntegral();
	void computeIndicator();
	void finalize();
};

template<>
InputParameters validParams<VariableJumpIndicator>();
