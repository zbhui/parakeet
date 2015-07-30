
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

	virtual void boundaryCondition();

	void fluxRiemann();
	virtual Real computeQpResidual(unsigned int p);
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);

	virtual void precalculateResidual();
	virtual void precalculateJacobian();

	virtual void wallBC(Real *ur);
	virtual void farFieldBC(Real *ur);
	virtual void symmetricBC(Real *ur);

	virtual void wallBC();
	virtual void farFieldBC();
	virtual void symmetricBC();
};

template<>
InputParameters validParams<EulerBC>();
