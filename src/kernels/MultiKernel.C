#include "MultiKernel.h"

template<>
InputParameters validParams<MultiKernel>()
{
	  InputParameters params = validParams<Kernel>();
//	  params += validParams<Kernel>();

	  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");
	  return params;
}

MultiKernel::MultiKernel(const std::string & name, InputParameters parameters):
		Kernel(name, parameters),
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

void MultiKernel::computeResidual()
{
	precalculateResidual();
	Kernel::computeResidual();
//	for (_eq = 0; _eq < _n_equation; ++_eq) {
//		DenseVector<Number> & re = _assembly.residualBlock(_sys.getVariable(_tid, _variables[_eq]).number());
//		_local_re.resize(re.size());
//		_local_re.zero();
//		for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//			for (_i = 0; _i < _test.size(); _i++)
//				_local_re(_i) += _JxW[_qp] * _coord[_qp] * computeQpResidual();
//
//		re += _local_re;
//
//		if (_has_save_in) {
//			Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
//			for (unsigned int i = 0; i < _save_in.size(); i++)
//				_save_in[i]->sys().solution().add_vector(_local_re,
//						_save_in[i]->dofIndices());
//		}
//	}
}


void MultiKernel::computeJacobian()
{
	precalculateJacobian();
	Kernel::computeJacobian();
//	DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
//	std::cout << ke <<std::endl << std::endl;


//	for (_ep = 0; _ep < _n_equation; ++_ep)
//	for (_eq = 0; _eq < _n_equation; ++_eq)
//	{
//	  int var_number_p = _sys.getVariable(_tid, _variables[_ep]).number();
//	  int var_number_q = _sys.getVariable(_tid, _variables[_eq]).number();
//	  DenseMatrix<Number> & ke = _assembly.jacobianBlock(var_number_p, var_number_q);
//	  _local_ke.resize(ke.m(), ke.n());
//
////	  std::cout << ke.m() << " " << ke.n() <<std::endl;
//	  _local_ke.zero();
//
//	  for (_i = 0; _i < _test.size(); _i++)
//	    for (_j = 0; _j < _phi.size(); _j++)
//	      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//	        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
//
//	  ke += _local_ke;
//
//	  if (_has_diag_save_in)
//	  {
//	    unsigned int rows = ke.m();
//	    DenseVector<Number> diag(rows);
//	    for (unsigned int i=0; i<rows; i++)
//	      diag(i) = _local_ke(i,i);
//
//	    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
//	    for (unsigned int i=0; i<_diag_save_in.size(); i++)
//	      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
//	  }
//	}
}

void MultiKernel::computeOffDiagJacobian(unsigned int jvar)
{
	precalculateJacobian();
	Kernel::computeOffDiagJacobian(jvar);

//	std::cout << _var.number() <<" " << jvar << std::endl;
//	std::cout << ke <<std::endl << std::endl;;
//	mooseError("MultiKernel::computeOffDiagJacobian暂不支持");
//
//	for (_eq = 0; _eq < _n_equation; ++_eq)
//	{
//		unsigned int var_number = _sys.getVariable(_tid, _variables[_eq]).number();
//		DenseMatrix<Number> & ke = _assembly.jacobianBlock(var_number, jvar);
//		if (jvar == var_number)
//		{
//			  _local_ke.resize(ke.m(), ke.n());
//			  _local_ke.zero();
//
//			  for (_i = 0; _i < _test.size(); _i++)
//			    for (_j = 0; _j < _phi.size(); _j++)
//			      for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//			        _local_ke(_i, _j) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
//
//			  ke += _local_ke;
//
//			  if (_has_diag_save_in)
//			  {
//			    unsigned int rows = ke.m();
//			    DenseVector<Number> diag(rows);
//			    for (unsigned int i=0; i<rows; i++)
//			      diag(i) = _local_ke(i,i);
//
//			    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
//			    for (unsigned int i=0; i<_diag_save_in.size(); i++)
//			      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
//			  }
//		}
//		else
//		{
//			for (_i=0; _i<_test.size(); _i++)
//				for (_j=0; _j<_phi.size(); _j++)
//					for (_qp=0; _qp<_qrule->n_points(); _qp++)
//					{
//						ke(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpOffDiagJacobian(jvar);
//					}
//		}
//	}
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
