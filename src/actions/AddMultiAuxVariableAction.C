
#include "AddMultiAuxVariableAction.h"
#include "AddKernelAction.h"
#include "FEProblem.h"
#include "Factory.h"
#include "MooseEnum.h"
#include "EigenSystem.h"

#include "libmesh/libmesh.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/string_to_enum.h"

template<>
InputParameters validParams<AddMultiAuxVariableAction>()
{
  InputParameters params = validParams<AddAuxVariableAction>();
  params += validParams<AddKernelAction>();
  params.addRequiredParam<std::vector<AuxVariableName> >("aux_variables", "非线性变量");
  params.addRequiredParam<std::string>("type", "A string representing the Moose Object that will be built by this Action");

  return params;
}

AddMultiAuxVariableAction::AddMultiAuxVariableAction(const InputParameters &params) :
	AddAuxVariableAction(params),
	_type(getParam<std::string>("type")),
	_aux_variables(getParam<std::vector<AuxVariableName> >("aux_variables"))
{
}

void AddMultiAuxVariableAction::act()
{
	for (int i = 0; i < _aux_variables.size(); ++i)
	{
		std::string var_name = _aux_variables[i];
		if (_current_task == "add_aux_variable")
			addAuxVariable(var_name);
	}
	if(_current_task == "add_aux_kernel")
		addAuxKernel(_type);
}

void AddMultiAuxVariableAction::addAuxKernel(std::string var_name)
{
	InputParameters params = _factory.getValidParams(_type);
	std::vector<VariableName> non_linear_var_name = _problem->getNonlinearSystem().getVariableNames();
	std::vector<NonlinearVariableName> variables;
	for(int i = 0; i< non_linear_var_name.size(); ++i)
		variables.push_back(non_linear_var_name[i]);

	_app.parser().extractParams(_name, params);
	params.set<AuxVariableName>("variable") = _aux_variables[0];
	params.set<std::vector<VariableName> >("variables") = non_linear_var_name;
	_problem->addAuxKernel(_type, var_name, params);
}

void AddMultiAuxVariableAction::addAuxVariable(std::string var_name)
{
	std::set<SubdomainID> blocks = getSubdomainIDs();

	if (_scalar_var)
		_problem->addAuxScalarVariable(var_name, _fe_type.order);

	else
	{
		if (_fe_type.order > 9)
			mooseError("Non-scalar AuxVariables must be CONSTANT, FIRST, SECOND, THIRD, FOURTH, FIFTH, SIXTH, SEVENTH, EIGHTH or NINTH order (" << _fe_type.order << " supplied)");

		if (blocks.empty())
			_problem->addAuxVariable(var_name, _fe_type);
		else
			_problem->addAuxVariable(var_name, _fe_type, &blocks);
	}
}
