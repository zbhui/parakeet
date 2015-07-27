
#include "ConservationLawBC.h"

template<>
InputParameters validParams<BurgersBC>()
{
	InputParameters params = validParams<IntegratedBC>();
	params += validParams<BurgersBase>();

	return params;
}

BurgersBC::BurgersBC(const std::string &name, InputParameters parameters):
		IntegratedBC(name, parameters),
		BurgersBase(name, parameters)
{
}

Real BurgersBC::computeQpResidual()
{
	Real ul = _u[_qp];
	Real ur = 0;
	Point normal = _normals[_qp];
	Real rflux = flux(ul, ur, normal);

	return  rflux * _test[_i][_qp];
//	return 0.;
}

Real BurgersBC::computeQpJacobian()
{
	 Real r = 0, ds = 1E-010;
	 Real rflux_new;
	 Real jacobi_ee;

	 Real ur = 0.;
	 Point normal = _normals[_qp];
	 Real rflux = flux(_u[_qp], ur, normal);

	 rflux_new = flux(_u[_qp]+ds, ur, normal);
	 jacobi_ee = (rflux_new-rflux)/ds;

//	 jacobi_ee = (0.5* 1*normal(0) + 0.1);

	 return jacobi_ee*_phi[_j][_qp]*_test[_i][_qp];
//	 return 0.;
}

Real BurgersBC::flux(Real ul, Real ur, Point normal)
{
	RealVectorValue fl, fr;
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);
	return  0.5*(fl(0)+fr(0))*normal(0)+ 0.1*(ul - ur);
}
