#include "ConservationLawCellKernel.h"

template<>
InputParameters validParams<ConservationLawCellKernel>()
{
	  InputParameters params = validParams<MultiKernel>();

	  return params;
}
ConservationLawCellKernel::ConservationLawCellKernel(const std::string & name, InputParameters parameters):
		MultiKernel(name, parameters)
{
}

Real ConservationLawCellKernel::computeQpResidual()
{
	return -_flux_vector[_qp][_eq]*_grad_test[_i][_qp];
}


