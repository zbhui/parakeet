#pragma once

#include "ElementIntegralPostprocessor.h"

class CFDProblem;

class ProblemElementalL1Error : public ElementIntegralPostprocessor
{
public:
	ProblemElementalL1Error(const InputParameters & parameters);
	~ProblemElementalL1Error(){}

	virtual Real getValue();

protected:
//	Real _test;
//	bool _is_something;
//	Real _max;

	virtual Real computeQpIntegral();

private:
	CFDProblem & _cfd_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	std::vector<VariableName> _variables;
	int _n_equations;

//	MooseEnum _error_type;

	std::vector<VariableValue*> _uh;


};


template<>
InputParameters validParams<ProblemElementalL1Error>();

