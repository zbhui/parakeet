
#pragma once

#include "TransientInterface.h"
#include "InternalSideIndicator.h"
#include "CFDDataPack.h"
#include <vector>
using std::vector;

class CFDProblem;

class FluxJumpIndicator :
public InternalSideIndicator,
public TransientInterface
{
public:
	FluxJumpIndicator(const InputParameters &parameters);
	virtual ~FluxJumpIndicator(){};

protected:
	CFDProblem &_cfd_problem;
	CFDDataPack _cfd_data, _cfd_data_neighbor;

	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<NonlinearVariableName> _variables;
	vector<VariableName> _aux_variables;
	int _n_variables;
	int _var_order;
	const Real &_current_elem_volume;
	const Real &_neighbor_elem_volume;
	const Real &_current_side_volume;

	vector<VariableValue*> _uh;
	vector<VariableValue*> _uh_neighbor;
	vector<VariableGradient*> _grad_uh;
	vector<VariableGradient*> _grad_uh_neighbor;

	virtual Real computeQpIntegral();
	void computeIndicator();
	void finalize();
};

template<>
InputParameters validParams<FluxJumpIndicator>();
