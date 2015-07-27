
#include "BurgersCellKernel.h"

template<>
InputParameters validParams<BurgersCellKernel>()
{
  InputParameters params = validParams<FDKernel>();
  params += validParams<BurgersBase>();

  return params;
}

BurgersCellKernel::BurgersCellKernel(const InputParameters & parameters):
		FDKernel(parameters),
		BurgersBase(parameters)
{
}

Real BurgersCellKernel::computeQpResidual()
{
	RealVectorValue inviscous_term;
	inviscousTerm(inviscous_term, _u[_qp]);
	Point xyz = _q_point[_qp];
	Real x = xyz(0);
//	Real source = -2*pi*sin(2*pi*x)*cos(2*pi*x);
	Real source = -2*pi *cos(2*pi*x);
	return -inviscous_term*_grad_test[_i][_qp]+source*_test[_i][_qp];
}

Real BurgersCellKernel::computeQpJacobian()
{
//	RealVectorValue jacobi(_u[_qp], _u[_qp], _u[_qp]);
	Real ds = 1E-08;
	RealVectorValue inviscous_term, inviscous_term_new;
	inviscousTerm(inviscous_term, _u[_qp]);
	inviscousTerm(inviscous_term_new, _u[_qp]+ds*_phi[_j][_qp]);
	RealVectorValue jacobi = (inviscous_term_new - inviscous_term)/ds;
//	std::cout << jacobi <<std::endl;
//	RealVectorValue jacobi(1,1,1);
//	return -jacobi*_grad_test[_i][_qp]*_phi[_j][_qp];

	return -jacobi*_grad_test[_i][_qp];
}
