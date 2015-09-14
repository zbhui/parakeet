
#pragma once

#include "ElementIntegralPostprocessor.h"

class CFDProblem;

using std::vector;

class ProblemElementalL2Error : public ElementIntegralPostprocessor
{
public:
	ProblemElementalL2Error(const InputParameters &parameters);
	virtual Real getValue();

protected:
	virtual Real computeQpIntegral();

private:
	CFDProblem & _cfd_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	vector<VariableName> _variables;
	int _n_equations;

	std::vector<VariableValue*> _uh;
};

template<>
InputParameters validParams<ProblemElementalL2Error>();
