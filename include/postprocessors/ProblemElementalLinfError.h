#pragma once

#include "ElementIntegralPostprocessor.h"

class CFDProblem;

class ProblemElementalLinfError : public ElementIntegralPostprocessor
{
public:
	ProblemElementalLinfError(const InputParameters & parameters);
	~ProblemElementalLinfError(){}

	virtual Real getValue();

protected:

//	bool _is_something;
	Real _max;
	Real _h;

	virtual void initialize();
	virtual Real computeQpIntegral();
	virtual Real computeIntegral();
	virtual void execute();
	virtual void threadJoin(const UserObject & y);

	CFDProblem & _cfd_problem;
	NonlinearSystem &_nl;
	THREAD_ID _tid;
	std::vector<VariableName> _variables;
	int _n_equations;

	std::vector<VariableValue*> _uh;


};


template<>
InputParameters validParams<ProblemElementalLinfError>();

