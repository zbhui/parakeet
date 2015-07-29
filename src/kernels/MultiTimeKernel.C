//
//#include "MultiTimeKernel.h"
//
//template<>
//InputParameters validParams<MultiTimeKernel>()
//{
//  InputParameters params = validParams<MultiKernel>();
//  return params;
//}
//
//MultiTimeKernel::MultiTimeKernel(const InputParameters & parameters) :
//    MultiKernel(parameters)
//{
//}
//
//void MultiTimeKernel::computeResidual()
//{
//	precalculateResidual();
//	for (_eq = 0; _eq < _n_equation; ++_eq)
//	{
//		DenseVector<Number> & re = _assembly.residualBlock(_sys.getVariable(_tid, _variables[_eq]).number(), Moose::KT_TIME);
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
//				_save_in[i]->sys().solution().add_vector(_local_re, _save_in[i]->dofIndices());
//		}
//	}
//}
