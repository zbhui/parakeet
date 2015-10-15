
#include "CFDInitialCondition.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDInitialCondition>()
{
  InputParameters params = validParams<InitialCondition>();
  params.addRequiredParam<unsigned int>("component", "分量");
  params.addParam<bool>("constant_ic", false, "常数初值");
  return params;
}

CFDInitialCondition::CFDInitialCondition(const InputParameters & parameters) :
    InitialCondition(parameters),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
    _component(getParam<unsigned int>("component")),
	_constant_ic(getParam<bool>("constant_ic"))
{}

void CFDInitialCondition::compute()
{
	InitialCondition::compute();

	if(!_constant_ic)
		return;

	NumericVector<Number> & solution = _var.sys().solution();
	std::vector<dof_id_type> dof_indices;
	dof_indices = _var.dofIndices();

	Number fineval = value(_current_elem->centroid());
	solution.set(dof_indices[0], fineval);
	const unsigned int n_dofs = dof_indices.size();
	for (unsigned int i = 1; i < n_dofs; i++)
	{
		solution.set(dof_indices[i], 0);
		_var.setNodalValue(0, i);
	}
}

Real CFDInitialCondition::value(const Point & p)
{
	return _cfd_problem.initialCondition(p, _component);
//	  return value(_component, p);
}

Real CFDInitialCondition::value(int component, const Point & p)
{
	return _cfd_problem.initialCondition(p, component);
}
