
#pragma once

#include "KOmegaModelBase.h"
#include "MultiKernel.h"

class KOmegaCellKernel;

template<>
InputParameters validParams<KOmegaCellKernel>();

class KOmegaCellKernel :
public MultiKernel,
public KOmegaModelBase
{
public:
	KOmegaCellKernel(const InputParameters & parameters);
	virtual ~KOmegaCellKernel(){}

protected:
	RealVectorValue _flux_vector[40][10];
	Real _source_term[40][10];

	virtual void precalculateResidual();

	virtual Real computeQpResidual(unsigned int p);
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);
};
