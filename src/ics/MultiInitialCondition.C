
#include "MultiInitialCondition.h"

template<>
InputParameters validParams<MultiInitialCondition>()
{
  InputParameters params = validParams<InitialCondition>();
  return params;
}

MultiInitialCondition::MultiInitialCondition(const InputParameters & parameters) :
    InitialCondition(parameters),
	_constant_ic(getParam<bool>("constant_ic"))
{}

void MultiInitialCondition::compute()
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

Real MultiInitialCondition::value(const Point & p)
{
	return 0;//_cfd_problem.initialCondition(p, _component);
}

Real MultiInitialCondition::value(int component, const Point & p)
{
	return 0;//_cfd_problem.initialCondition(p, component);
}
