
#pragma once

#include "CFDProblem.h"

class EulerProblem : public CFDProblem
{
public:
	EulerProblem(const InputParameters &params);
};

template<>
InputParameters validParams<EulerProblem>();
