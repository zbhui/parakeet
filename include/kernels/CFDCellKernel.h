#pragma once

#include "MultiKernel.h"
#include "CFDDataPack.h"

class CFDProblem;

class CFDCellKernel : public MultiKernel
{
public:
	CFDCellKernel(const InputParameters & parameters);
	virtual ~CFDCellKernel(){}

protected:
	CFDProblem &_cfd_problem;
	CFDDataPack _cfd_data;
	RealVectorValue _flux[10], _flux_old[10];
	RealVectorValue _viscous[10], _viscous_old[10];
	RealVectorValue _artificial_viscous[10], _artificial_viscous_old[10];
	RealVectorValue _flux_jacobi_variable[10][10];
	RealTensorValue _flux_jacobi_grad_variable[10][10];

	Real _source[10], _source_old[10];
	Real _source_jacobi_variable[10][10];
	RealVectorValue _source_jacobi_grad_variable[10][10];

	Real _perturbation;

	virtual void precalculateResidual();
	virtual Real computeQpResidual(unsigned int p);

	virtual void precalculateJacobian();
	virtual Real computeQpJacobian(unsigned int p, unsigned int q);

	void fluxTerm();
	void reinitArtificialViscous();
	void reinit();
	void reinitViscous();

	ErrorVector & getErrorVector(std::string indicator);
	ErrorVector * _error_vector;
	bool _has_artificial_vis;
};


template<>
InputParameters validParams<CFDCellKernel>();
