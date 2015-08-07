
#include "AddMultiVariable.h"
//#include "FEProblem.h"

template<>
InputParameters validParams<AddMultiVariable>()
{
  InputParameters params = validParams<AddVariableAction>();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "非线性变量");

  return params;
}

AddMultiVariable::AddMultiVariable(InputParameters params) :
	AddVariableAction(params)
{
}

void AddMultiVariable::act()
{
	std::vector<NonlinearVariableName> variables = getParam<std::vector<NonlinearVariableName> >("variables");
	for (int i = 0; i < variables.size(); ++i)
	{
		std::string var_name = variables[i];
		addVariable(var_name);
	}
}
