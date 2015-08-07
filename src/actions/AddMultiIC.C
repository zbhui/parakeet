
#include "AddMultiIC.h"
#include "FEProblem.h"

template<>
InputParameters validParams<AddMultiIC>()
{
  InputParameters params = validParams<MooseObjectAction>();
  return params;
}

AddMultiIC::AddMultiIC(InputParameters params) :
	MooseObjectAction(params)
{
}

void AddMultiIC::act()
{
	std::vector<VariableName> var = _problem->getNonlinearSystem().getVariableNames();
    _app.parser().extractParams(_name, _moose_object_pars);
	for (int i = 0; i < var.size(); ++i)
	{
	    _moose_object_pars.set<VariableName>("variable") = var[i];
	    _moose_object_pars.set<unsigned int>("component") = i;
	    _problem->addInitialCondition(_type, var[i]+"_ic", _moose_object_pars);
	}
}

