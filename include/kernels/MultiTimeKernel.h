
#pragma once

#include "MultiKernel.h"

// Forward Declaration
class MultiTimeKernel;

template<>
InputParameters validParams<MultiTimeKernel>();

/**
 * All Multi-time kernels should inherit from this class
 *
 */
class MultiTimeKernel : public MultiKernel
{
public:
	MultiTimeKernel(const std::string & name, InputParameters parameters);
	virtual ~MultiTimeKernel(){};

	virtual void computeResidual();

};

