//
//#include "MultiTimeDerivative.h"
//
//template<>
//InputParameters validParams<MultiTimeDerivative>()
//{
//  InputParameters params = validParams<MultiTimeKernel>();
//  params.addParam<bool>("lumping", false, "True for mass matrix lumping, false otherwise");
//  return params;
//}
//
//MultiTimeDerivative::MultiTimeDerivative(const InputParameters & parameters) :
//    MultiTimeKernel(parameters),
//    _lumping(getParam<bool>("lumping"))
//{
//}
//
//Real
//MultiTimeDerivative::computeQpResidual()
//{
//  return _test[_i][_qp] * _duh_dt[_qp][_eq];
//}
//
//Real
//MultiTimeDerivative::computeQpJacobian()
//{
//  return _test[_i][_qp]*_phi[_j][_qp]*_du_dot_du[_qp];
//}
//
//void
//MultiTimeDerivative::computeJacobian()
//{
//  if (_lumping)
//  {
//    DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
//
//    for (_i = 0; _i < _test.size(); _i++)
//      for (_j = 0; _j < _phi.size(); _j++)
//        for (_qp = 0; _qp < _qrule->n_points(); _qp++)
//        {
//          ke(_i, _i) += _JxW[_qp] * _coord[_qp] * computeQpJacobian();
//        }
//  }
//  else
//  {
//	  mooseError("erro");
//  }
////    TimeKernel::computeJacobian();
//}
//
//void MultiTimeDerivative::precalculateResidual()
//{
//	for (_qp = 0; _qp < _qrule->n_points(); _qp++) {
//		for (_eq = 0; _eq < _n_equation; ++_eq)
//		{
//			_duh_dt[_qp][_eq] = (*_uh_dot[_eq])[_qp];
//		}
//	}
//}
