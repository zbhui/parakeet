
#include "EmptyTimeDerivative.h"

template<>
InputParameters validParams<EmptyTimeDerivative>()
{
  InputParameters params = validParams<TimeDerivative>();
  return params;
}

EmptyTimeDerivative::EmptyTimeDerivative(const InputParameters &parameters) :
	TimeDerivative(parameters)
{
}


void EmptyTimeDerivative::computeOffDiagJacobian(unsigned int jvar)
{
//	Moose::perf_log.push("computeOffDiagJacobian()","TimeDerivative");
//	if (jvar == _var.number())
//		computeJacobian();
	TimeDerivative::computeOffDiagJacobian(jvar);
//	Moose::perf_log.pop("computeOffDiagJacobian()","TimeDerivative");
}
