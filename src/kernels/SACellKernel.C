
#include "SACellKernel.h"

template<>
InputParameters validParams<SACellKernel>()
{
	  InputParameters params = validParams<Kernel>();
	  params += validParams<CFDBase>();

	  return params;
}
SACellKernel::SACellKernel(const InputParameters & parameters):
		Kernel(parameters),
		SAModelBase(parameters)
{
}

Real SACellKernel::computeQpResidual()
{
	return _u[_qp];
}

Real SACellKernel::computeQpJacobian()
{
	return 0.;
}

Real SACellKernel::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.;
}

