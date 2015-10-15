
#include "MultiVariableInterface.h"
#include "Assembly.h"
#include "MooseVariable.h"
#include "FEProblem.h"
#include "SubProblem.h"
#include "SystemBase.h"
#include "MaterialData.h"
#include "ParallelUniqueId.h"

template<>
InputParameters validParams<MultiVariableInterface>()
{
  InputParameters params = emptyInputParameters();
  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");

  return params;
}

MultiVariableInterface::MultiVariableInterface(const InputParameters & parameters) :
    _mv_fe_problem(*parameters.getCheckedPointerParam<FEProblem *>("_fe_problem")),
    _sys(*parameters.get<SystemBase *>("_sys")),
    _tid(parameters.get<THREAD_ID>("_tid")),
	_variables(parameters.get<std::vector<NonlinearVariableName> >("variables")),
	_var_order(_mv_fe_problem.getVariable(_tid, _variables[0]).order())
{
	MooseVariable &val0 = _sys.getVariable(_tid, _variables[0]);
	_n_equation = _variables.size();
	for (size_t i = 0; i < _n_equation; ++i)
	{
		MooseVariable &val = _sys.getVariable(_tid, _variables[i]);
		if(val.feType() != val0.feType()) mooseError("MultiDGKernel中变量的类型不一致");

		_uh.push_back(&val.sln());
		_uh_neighbor.push_back(&val.slnNeighbor());
		_grad_uh.push_back(&val.gradSln());
		_grad_uh_neighbor.push_back(&val.gradSlnNeighbor());
	}
}
