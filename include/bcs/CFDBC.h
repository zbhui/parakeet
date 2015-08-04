
#pragma once

#include "MooseEnum.h"

#include "MultiIntegratedBC.h"
#include "CFDDataPack.h"
#include "CFDBase.h"

class CFDProblem;

class CFDBC : public MultiIntegratedBC
{
public:
	CFDBC(const InputParameters & params);
	virtual ~CFDBC(){}

protected:
	CFDProblem &_cfd_problem;

	CFDDataPack _cfd_data, _cfd_data_neighbor;

	Real _flux[10], _flux_old[10];
	Real _ul[10], _ur[10];
	Real _jacobi_variable[10][10];
	Real _perturbation;

	virtual void boundaryCondition();

	virtual void precalculateResidual();
	virtual Real computeQpResidual(unsigned int p);

	virtual void precalculateJacobian();
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);
	void fluxRiemann();
};

template<>
InputParameters validParams<CFDBC>();
