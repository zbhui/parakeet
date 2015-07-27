
#include "BurgersFaceKernel.h"

template<>
InputParameters validParams<BurgersFaceKernel>()
{
	InputParameters params = validParams<DGKernel>();
	params += validParams<BurgersBase>();

	return params;
}

BurgersFaceKernel::BurgersFaceKernel(const std::string &name, InputParameters parameters):
		DGKernel(name, parameters),
		BurgersBase(name, parameters)
{

}

Real BurgersFaceKernel::computeQpResidual(Moose::DGResidualType type)
{
	Point normal = _normals[_qp];
	Real rflux = flux(_u[_qp], _u_neighbor[_qp], normal);

	switch(type)
	{
	case Moose::Element:
		return  rflux * _test[_i][_qp];
		break;
	case Moose::Neighbor:
		return  -rflux * _test_neighbor[_i][_qp];
		break;
	}
	return 0;
}

Real BurgersFaceKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	  Real r = 0, ds = 1E-08;
	  Real rflux_new;
	  Real jacobi_ee, jacobi_en, jacobi_ne, jacobi_nn;

	  Point normal = _normals[_qp];
	  Real rflux = flux(_u[_qp], _u_neighbor[_qp], normal);

	   rflux_new = flux(_u[_qp]+ds*_phi[_j][_qp], _u_neighbor[_qp], normal);
	   jacobi_ee = (rflux_new-rflux)/ds;
	   jacobi_ne = -(rflux_new-rflux)/ds;

	   rflux_new = flux(_u[_qp], _u_neighbor[_qp]+ds*_phi_neighbor[_j][_qp], normal);
	   jacobi_en = (rflux_new-rflux)/ds;
	   jacobi_nn = -(rflux_new-rflux)/ds;

//	   std::cout << jacobi_ee -(0.5*1*normal(0) +0.1) <<std::endl;
//	   std::cout << jacobi_en - (0.5*_u_neighbor[_qp]*normal(0) - 0.8)<<std::endl;

//	   jacobi_ee = (0.5* 1*normal(0) + 0.1);
//	   jacobi_ne = -0.5* 1*normal(0) - 0.1;
//	   jacobi_en =  (0.5*1*normal(0) - 0.1);
//	   jacobi_nn = - 0.5*1*normal(0) + 0.1;

	   switch (type)
	  {
	  case Moose::ElementElement:
//		  r = jacobi_ee*_phi[_j][_qp]*_test[_i][_qp];
		  r = jacobi_ee*_test[_i][_qp];
	    break;

	  case Moose::ElementNeighbor:
//		  r = jacobi_en*_phi_neighbor[_j][_qp]*_test[_i][_qp];
		  r = jacobi_en*_test[_i][_qp];
	    break;

	  case Moose::NeighborElement:
//		  r = jacobi_ne*_phi[_j][_qp]*_test_neighbor[_i][_qp];
		  r = jacobi_ne*_test_neighbor[_i][_qp];
	    break;

	  case Moose::NeighborNeighbor:
//		  r = jacobi_nn*_phi_neighbor[_j][_qp]*_test_neighbor[_i][_qp];
		  r = jacobi_nn*_test_neighbor[_i][_qp];
	    break;
	  }
	  return r;
}

Real BurgersFaceKernel::flux(Real ul, Real ur, Point normal)
{
	RealVectorValue fl, fr;
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);
	return  0.5*(fl(0)+fr(0))*normal(0)+ 0.1*(ul - ur);
}
