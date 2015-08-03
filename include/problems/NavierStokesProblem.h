
#pragma once

#include "CFDProblem.h"

class NavierStokesProblem : public CFDProblem
{
public:
	NavierStokesProblem(const InputParameters &params);

};

template<>
InputParameters validParams<NavierStokesProblem>();
