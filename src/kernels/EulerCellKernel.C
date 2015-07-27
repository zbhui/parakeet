#include "EulerCellKernel.h"

template<>
InputParameters validParams<EulerCellKernel>()
{
  InputParameters params = validParams<MultiKernel>();
  params += validParams<CFDBase>();

  return params;
}

EulerCellKernel::EulerCellKernel(const InputParameters & parameters):
		MultiKernel(parameters),
		CFDBase(parameters)
{
	std::string var_name = _var.name();

	if(var_name == "rho")
		_eq = 0;
	if(var_name == "momentum_x")
		_eq = 1;
	if(var_name == "momentum_y")
		_eq = 2;
	if(var_name == "momentum_z")
		_eq = 3;
	if(var_name == "rhoe")
		_eq = 4;
}

void EulerCellKernel::precalculateResidual()
{
	mooseAssert(_n_equation < 10, "multiKernel方程个数应<10");
	mooseAssert(_qrule->n_points() < 40, "mulitKernel积分点个数应<40");

	Real uh[10];
	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		valueAtCellPoint(uh);
		inviscousTerm(_flux[_qp], uh);
	}
}

void EulerCellKernel::precalculateJacobian()
{
	Real uh[10];
	RealVectorValue flux_vector_new[10], flux_vector[10];

	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		valueAtCellPoint(uh);
		inviscousTerm(flux_vector, uh);
		for (int q = 0; q < _n_equation; ++q)
		{
			uh[q] += _ds;
			inviscousTerm(flux_vector_new, uh);
			for (int p = 0; p < _n_equation; ++p)
				_jacobi_variable[_qp][p][q] = (flux_vector_new[p] - flux_vector[p])/_ds;

			uh[q] -= _ds;
		}
	}
}

Real EulerCellKernel::computeQpJacobian()
{
	Real uh[10];
	Matrix5x5 jacobi[3];

	valueAtCellPoint(uh);
	inviscousJacobian(jacobi, uh);
	for (int p = 0; p < _n_equation; ++p)
	{
		for (int q = 0; q < _n_equation; ++q)
		{
			RealVectorValue tmp(jacobi[0](p,q), jacobi[1](p,q), jacobi[2](p,q));
			_jacobi_variable[_qp][p][q] = tmp;
		}
	}

	return -(_jacobi_variable[_qp][_eq][_eq]*_phi[_j][_qp])*_grad_test[_i][_qp];
//	return -(_jacobi_variable[_qp][_eq][_eq])*_grad_test[_i][_qp];
}

Real EulerCellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	Real uh[10];
	Matrix5x5 jacobi[3];

	valueAtCellPoint(uh);
	inviscousJacobianFD(jacobi, uh);
	for (int p = 0; p < _n_equation; ++p)
	{
		for (int q = 0; q < _n_equation; ++q)
		{
			RealVectorValue tmp(jacobi[0](p,q), jacobi[1](p,q), jacobi[2](p,q));
			_jacobi_variable[_qp][p][q] = tmp;
		}
	}
	return -(_jacobi_variable[_qp][_eq][jvar]*_phi[_j][_qp])*_grad_test[_i][_qp];
}

Real EulerCellKernel::computeQpResidual()
{
	return -_flux[_qp][_eq]*_grad_test[_i][_qp];
}


