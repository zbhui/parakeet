
#pragma once

#include "MooseEnum.h"

#include "MultiIntegratedBC.h"
#include "CFDDataPack.h"

class CFDProblem;
class Attitude;

class CFDBC : public MultiIntegratedBC
{
public:
	CFDBC(const InputParameters & params);
	virtual ~CFDBC(){}

protected:
	CFDProblem &_cfd_problem;
	CFDDataPack _cfd_data, _cfd_data_neighbor, _lift_data;
	MooseEnum _flux_type;

	Real _flux[10], _flux_old[10];
	Real _flux_jacobi_variable[10][10];
	RealVectorValue _flux_jacobi_grad_variable[10][10];
	RealVectorValue _lift[10], _lift_old[10];
	RealVectorValue _lift_jacobi_variable[10][10];
	Real _perturbation;
	Real _penalty;
	Real _gamma;
	Real _mach;
	Attitude &_attitude;

	virtual void boundaryCondition();

	virtual void precalculateResidual();
	virtual Real computeQpResidual(unsigned int p);

	virtual void precalculateJacobian();
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);
	void reinit();
	void fluxRiemann();
	void fluxLaxF();
	void fluxHLLCPV();
	void liftOperator();
};

template<>
InputParameters validParams<CFDBC>();
