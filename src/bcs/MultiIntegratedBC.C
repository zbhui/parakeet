
#include "MultiIntegratedBC.h"

template<>
InputParameters validParams<MultiIntegratedBC>()
{
	InputParameters params = validParams<IntegratedBC>();

	params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");
	return params;
}
MultiIntegratedBC::MultiIntegratedBC(const std::string& name, InputParameters parameters):
		IntegratedBC(name, parameters),
		_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
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
	precalculateResidual();
	IntegratedBC::computeResidual();

//	for (_eq = 0; _eq < _n_equation; ++_eq)
//	{
//		DenseVector<Number> & re = _assembly.residualBlock(_sys.getVariable(_tid, _variables[_eq]).number());
//		_local_re.resize(re.size());
//		_local_re.zero();
//
//		for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//			for (_i = 0; _i < _test.size(); _i++)
//				_local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();
//
//		re += _local_re;
//
//		if (_has_save_in)
//		{
//		    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
//		    for(unsigned int i=0; i<_save_in.size(); i++)
//		      _save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
//		}
//	}
}

void MultiIntegratedBC::valueAtLeftFace(Real* ul)
{
	for (size_t eq = 0; eq < _uh.size(); ++eq)
	{
		ul[eq] = (*_uh[eq])[_qp];
	}
}

void MultiIntegratedBC::valueAtRightFace(Real* ur)
{
	valueAtLeftFace(ur);
}

void MultiIntegratedBC::valueGradAtLeftFace(RealGradient* dul)
{
	for (size_t eq = 0; eq < _grad_uh.size(); ++eq)
	{
		dul[eq] = (*_grad_uh[eq])[_qp];
	}
}

void MultiIntegratedBC::valueGradAtRightFace(RealGradient* dur)
{
	valueGradAtLeftFace(dur);
}

void MultiIntegratedBC::computeJacobian()
{
	IntegratedBC::computeJacobian();

//	for (_ep = 0; _ep < _n_equation; ++_ep)
//	for (_eq = 0; _eq < _n_equation; ++_eq)
//	{
//		int var_number_p = _sys.getVariable(_tid, _variables[_ep]).number();
//		int var_number_q = _sys.getVariable(_tid, _variables[_eq]).number();
//		DenseMatrix<Number> & ke = _assembly.jacobianBlock(var_number_p, var_number_q);
//		_local_ke.resize(ke.m(), ke.n());
//		_local_ke.zero();
//
//		for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//			for (_i = 0; _i < _test.size(); _i++)
//				 for (_j = 0; _j < _phi.size(); _j++)
//					 _local_ke(_i, _j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian();
//
//		  ke += _local_ke;
//
//		  if (_has_diag_save_in)
//		  {
//			  unsigned int rows = ke.m();
//			  DenseVector<Number> diag(rows);
//			  for (unsigned int i=0; i<rows; i++)
//				  diag(i) = _local_ke(i,i);
//
//			  Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
//			  for (unsigned int i=0; i<_diag_save_in.size(); i++)
//				  _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
//		  }
//	}
}

void MultiIntegratedBC::computeJacobianBlock(unsigned int jvar)
{
	IntegratedBC::computeJacobianBlock(jvar);

//	mooseError("MultiIntegratedBC::computeJacobianBlock暂不支持");
//
//	for (_eq = 0; _eq < _n_equation; ++_eq)
//	{
//		int var_number = _sys.getVariable(_tid, _variables[_eq]).number();
//		DenseMatrix<Number> & ke = _assembly.jacobianBlock(var_number, jvar);
//
//		for (_qp=0; _qp<_qrule->n_points(); _qp++)
//			for (_i=0; _i<_test.size(); _i++)
//				for (_j=0; _j<_phi.size(); _j++)
//				{
//					if (var_number == jvar)
//						ke(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian();
//					else
//						ke(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpOffDiagJacobian(jvar);
//				}
//		}
}
