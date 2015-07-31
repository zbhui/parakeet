
#pragma once
#include "CFDBC.h"
#include "CFDDataPack.h"

class EulerBC :
		public CFDBC,
		public CFDBase
{
public:
	EulerBC(const InputParameters & params);
	virtual ~EulerBC(){}


protected:
	CFDDataPack _cfd_data, _cfd_data_neighbor;

	Real _flux[10], _flux_old[10];
	Real _ul[10], _ur[10];
	Real _jacobi_variable[10][10];

	virtual void precalculateResidual();
	virtual Real computeQpResidual(unsigned int p);

	virtual void precalculateJacobian();
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);

	void fluxRiemann();
	virtual void boundaryCondition();
	virtual void wall();
	virtual void farfield();
	virtual void symmetric();
};

template<>
InputParameters validParams<EulerBC>();
