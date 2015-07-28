#include "MultiKernel.h"

template<>
InputParameters validParams<MultiKernel>()
{
	  InputParameters params = validParams<Kernel>();
//	  params += validParams<Kernel>();

	  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");
	  return params;
}

MultiKernel::MultiKernel(const InputParameters & parameters):
		Kernel(parameters),
		_variables(parameters.get<std::vector<NonlinearVariableName> >("variables"))
{
	MooseVariable &val0 = _sys.getVariable(_tid, _variables[0]);
	_n_equation = _variables.size();
	for (size_t i = 0; i < _n_equation; ++i)
	{
		MooseVariable &val = _sys.getVariable(_tid, _variables[i]);
		if(val.feType() != val0.feType()) mooseError("multiKernel中变量的类型不一致");

		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
		_grad_uh.push_back(_is_implicit ? &val.gradSln() : &val.gradSlnOld());
		_uh_dot.push_back(&val.uDot());
		_duh_dot_du.push_back(&val.duDotDu());
	}
}

Real MultiKernel::computeQpResidual()
{
	mooseError("MultiKernel::computeQpResidual不能调用");
	return 0;
}
void MultiKernel::computeResidual()
{
	precalculateResidual();
	for (unsigned int p = 0; p < _n_equation; ++p)
	{
		DenseVector<Number> & re = _assembly.residualBlock(p);
		_local_re.resize(re.size());
		_local_re.zero();
		for (_qp = 0; _qp < _qrule->n_points(); _qp++)
			for (_i = 0; _i < _test.size(); _i++)
				_local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual(p);

		re += _local_re;
	}
}


void MultiKernel::computeJacobian()
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
	      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian(p, q);

	  ke += _local_ke;
	}
}

void MultiKernel::computeOffDiagJacobian(unsigned int jvar)
{
	 if (jvar == _var.number())
	    computeJacobian();
	 else
	 {
		 return;
		 mooseError("MultiKernel::computeOffDiagJacobian");
	 }
}

void MultiKernel::computeOffDiagJacobianScalar(unsigned int jvar)
{
	mooseError("MultiKernel::computeOffDiagJacobianScalar");
}

void MultiKernel::valueAtCellPoint(Real *uh)
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		uh[i] = (*_uh[i])[_qp];
	}
}

void MultiKernel::valueGradAtCellPoint(RealGradient* duh)
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		duh[i] = (*_grad_uh[i])[_qp];
	}
}

void MultiKernel::precalculateJacobian()
{
}
