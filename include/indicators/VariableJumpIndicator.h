
#pragma once

#include "FluxJumpIndicator.h"
#include <vector>
using std::vector;

class VariableJumpIndicator : public FluxJumpIndicator
{
public:
	VariableJumpIndicator(const InputParameters &parameters);
	virtual ~VariableJumpIndicator(){};

protected:
	virtual Real computeQpIntegral();
	void computeIndicator();
	void finalize();
};

template<>
InputParameters validParams<VariableJumpIndicator>();
