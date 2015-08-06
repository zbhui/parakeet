
#include "AddMultiKernel.h"
#include "FEProblem.h"

template<>
InputParameters validParams<AddMultiKernel>()
{
	InputParameters params = validParams<MooseObjectAction>();
	return params;
}

AddMultiKernel::AddMultiKernel(InputParameters params) :
	MooseObjectAction(params)
{
}

void AddMultiKernel::act()
{
	std::vector<VariableName> var_name = _problem->getNonlinearSystem().getVariableNames();
	std::vector<NonlinearVariableName> variables;
	for(int i = 0; i< var_name.size(); ++i)
		variables.push_back(var_name[i]);

    _app.parser().extractParams(_name, _moose_object_pars);
	_moose_object_pars.set<NonlinearVariableName>("variable") = variables[0];
	_moose_object_pars.set<std::vector<NonlinearVariableName> >("variables") = variables;

	_problem->addKernel(_type, _name, _moose_object_pars);

	for(int i = 0; i< var_name.size(); ++i)
	{
		std::string time_kernel_name = "TimeDerivative";
		InputParameters params = _factory.getValidParams(time_kernel_name);
		params.set<NonlinearVariableName>("variable") = variables[i];
		_problem->addKernel(time_kernel_name, variables[i] + "_time", params);
	}
}
