#include "ConservationLawCellKernel.h"

template<>
InputParameters validParams<ConservationLawCellKernel>()
{
	  InputParameters params = validParams<MultiKernel>();

	  return params;
}
ConservationLawCellKernel::ConservationLawCellKernel(const InputParameters &parameters):
		MultiKernel(parameters)
{
}

Real ConservationLawCellKernel::computeQpResidual()
{
	return 0;
//	return -_flux_vector[_qp][_eq]*_grad_test[_i][_qp];
}


