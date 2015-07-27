#pragma once

#include "MultiKernel.h"

// 前置声明
class ConservationLawCellKernel;

template<>
InputParameters validParams<ConservationLawCellKernel>();

class ConservationLawCellKernel :
public MultiKernel
{
public:
	ConservationLawCellKernel(const InputParameters & parameters);
	virtual ~ConservationLawCellKernel(){}

protected:
	virtual Real computeQpResidual();

	RealVectorValue _flux_vector[40][10];
protected:

};
