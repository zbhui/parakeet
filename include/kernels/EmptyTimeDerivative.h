
#pragma once

#include "TimeDerivative.h"

class EmptyTimeDerivative : public TimeDerivative
{
public:
	EmptyTimeDerivative(const InputParameters &parameters);

	virtual void computeOffDiagJacobian(unsigned int jvar);

};

template<>
InputParameters validParams<EmptyTimeDerivative>();
