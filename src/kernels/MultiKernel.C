#include "MultiKernel.h"

template<>
InputParameters validParams<MultiKernel>()
{
	  InputParameters params = validParams<Kernel>();
	  params += validParams<MultiVariableInterface>();
	  return params;
}

MultiKernel::MultiKernel(const InputParameters & parameters):
		Kernel(parameters),
		MultiVariableInterface(parameters)
{
}

void MultiKernel::computeResidual()
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


void MultiKernel::computeJacobian()
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
					_local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian(p, q);

			ke += _local_ke;
		}
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

Real MultiKernel::computeQpResidual()
{
	mooseError("MultiKernel::computeQpResidual不能调用");
	return 0;
}

Real MultiKernel::computeQpJacobian()
{
	mooseError("MultiKernel::computeQpJacobian不能调用");
	return 0;
}
