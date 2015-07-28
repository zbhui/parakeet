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

}

void EulerCellKernel::precalculateResidual()
{
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

Real EulerCellKernel::computeQpJacobian(unsigned int p, unsigned int q)
{
	return -(_jacobi_variable[_qp][p][q]*_phi[_j][_qp])*_grad_test[_i][_qp];
}

Real EulerCellKernel::computeQpResidual(unsigned int p)
{
	return -_flux[_qp][p]*_grad_test[_i][_qp];
}


