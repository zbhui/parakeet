#pragma once

#include "MultiKernel.h"
#include "CFDBase.h"
#include "CFDDataPack.h"

class EulerCellKernel;

template<>
InputParameters validParams<EulerCellKernel>();

class EulerCellKernel :
public MultiKernel,
public CFDBase
{
public:
	EulerCellKernel(const InputParameters & parameters);
	virtual ~EulerCellKernel(){}

protected:
	CFDDataPack _cfd_data;
	RealVectorValue _flux[10], _flux_old[10];
	RealVectorValue _jacobi_variable[10][10];

	void fluxTerm();
	virtual void precalculateResidual();
	virtual Real computeQpResidual(unsigned int p);

	virtual Real computeQpJacobian(unsigned int p, unsigned int q);
	virtual void precalculateJacobian();

};
