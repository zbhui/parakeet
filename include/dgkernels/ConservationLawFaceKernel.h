
#pragma once

#include "MultiDGKernel.h"

class ConservationLawFaceKernel;

template<>
InputParameters validParams<ConservationLawFaceKernel>();

class ConservationLawFaceKernel :
public MultiDGKernel
{
public:
	ConservationLawFaceKernel(const std::string &name, InputParameters parameters);
	virtual ~ConservationLawFaceKernel(){}

protected:

	virtual Real computeQpResidual(Moose::DGResidualType type);

	Real _flux[40][10];
};
