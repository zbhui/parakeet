
#include "MultiIntegratedBC.h"

template<>
InputParameters validParams<MultiIntegratedBC>()
{
	InputParameters params = validParams<IntegratedBC>();

	params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");
	return params;
}
MultiIntegratedBC::MultiIntegratedBC(const InputParameters & parameters):
		IntegratedBC(parameters),
		_variables(getParam<std::vector<NonlinearVariableName> >("variables")),
		_var_order(_fe_problem.getVariable(_tid, _variables[0]).order())
{
	MooseVariable &val0 = _sys.getVariable(_tid, _variables[0]);
	_n_equation = _variables.size();
	for (size_t i = 0; i < _n_equation; ++i)
	{
		MooseVariable &val = _sys.getVariable(_tid, _variables[i]);
		if(val.feType() != val0.feType()) mooseError("MultiIntegratedBC中变量的类型不一致");

		_uh.push_back(&val.sln());
		_grad_uh.push_back(&val.gradSln());
	}
}

void MultiIntegratedBC::computeResidual()
{
	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		precalculateResidual();
		for (unsigned int p = 0; p < _n_equation; ++p)
		{
			DenseVector<Number> & re = _assembly.residualBlock(p);
			_local_re.resize(re.size());
			_local_re.zero();

			for (_i = 0; _i < _test.size(); _i++)
				_local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual(p);

			re += _local_re;
		}
	}
}

void MultiIntegratedBC::computeJacobian()
{
	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		precalculateJacobian();

		for (unsigned int p = 0; p < _n_equation; ++p)
		for (unsigned int q = 0; q < _n_equation; ++q)
		{
			DenseMatrix<Number> & ke = _assembly.jacobianBlock(p, q);
			_local_ke.resize(ke.m(), ke.n());
			_local_ke.zero();

			for (_i = 0; _i < _test.size(); _i++)
				for (_j = 0; _j < _phi.size(); _j++)
					_local_ke(_i, _j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian(p, q);

			ke += _local_ke;
		}
	}
}

void MultiIntegratedBC::computeJacobianBlock(unsigned int jvar)
{
	if (jvar == _var.number())
	    computeJacobian();
	else
	{
		return;
		mooseError("MultiIntegratedBC::computeJacobianBlock暂不支持");
	}
}

Real MultiIntegratedBC::computeQpResidual()
{
	mooseError("MultiIntegratedBC::computeQpResidual");
}

Real MultiIntegratedBC::computeQpJacobian()
{
	mooseError("MultiIntegratedBC::computeQpJacobian");
}
