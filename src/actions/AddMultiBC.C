
#include "AddMultiBC.h"
#include "FEProblem.h"

template<>
InputParameters validParams<AddMultiBC>()
{
	InputParameters params = validParams<MooseObjectAction>();
	return params;
}

AddMultiBC::AddMultiBC(InputParameters params) :
		MooseObjectAction(params)
{
}

void AddMultiBC::act()
{
	std::vector<VariableName> var_name = _problem->getNonlinearSystem().getVariableNames();
	std::vector<NonlinearVariableName> variables;
	for(int i = 0; i< var_name.size(); ++i)
		variables.push_back(var_name[i]);

	_moose_object_pars.set<NonlinearVariableName>("variable") = variables[0];
	_moose_object_pars.set<std::vector<NonlinearVariableName> >("variables") = variables;

    _problem->addBoundaryCondition(_type, _name, _moose_object_pars);
}
