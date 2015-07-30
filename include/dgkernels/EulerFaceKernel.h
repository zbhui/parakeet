
#pragma once

#include "MultiDGKernel.h"
#include "CFDBase.h"
#include "CFDDataPack.h"

class EulerFaceKernel;

template<>
InputParameters validParams<EulerFaceKernel>();

class EulerFaceKernel :
public MultiDGKernel,
public CFDBase
{
public:
	EulerFaceKernel(const InputParameters & parameters);
	virtual ~EulerFaceKernel(){}

protected:
	CFDDataPack _cfd_data, _cfd_data_neighbor;

	Real _flux[10], _flux_old[10];

	Real _jacobi_variable_ee[10][10];
	Real _jacobi_variable_en[10][10];

	virtual Real computeQpResidual(Moose::DGResidualType type, unsigned int p);
	virtual Real computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q);

	void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	void fluxRiemann();
	virtual void precalculateResidual();
	virtual void precalculateJacobian();

};
