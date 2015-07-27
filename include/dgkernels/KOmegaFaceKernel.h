#pragma once


#pragma once

#include "MultiDGKernel.h"
#include "KOmegaModelBase.h"

class KOmegaFaceKernel;

template<>
InputParameters validParams<KOmegaFaceKernel>();

class KOmegaFaceKernel :
public MultiDGKernel,
public KOmegaModelBase
{
public:
	KOmegaFaceKernel(const InputParameters & parameters);
	virtual ~KOmegaFaceKernel(){}

protected:
	Real _flux[40][10];
	RealVectorValue _penalty[40][10];
	RealVectorValue _penalty_neighbor[40][10];

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	virtual void precalculateResidual();

	virtual void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
};





