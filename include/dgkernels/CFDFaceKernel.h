
#pragma once

#include "MultiDGKernel.h"
#include "CFDDataPack.h"

class CFDProblem;

class CFDFaceKernel : public MultiDGKernel
{
public:
	CFDFaceKernel(const InputParameters & parameters);
	virtual ~CFDFaceKernel(){}

protected:
	CFDProblem &_cfd_problem;
	CFDDataPack _cfd_data, _cfd_data_neighbor, _lift_data;
	MooseEnum _flux_type;

	Real _flux[10], _flux_old[10];
	Real _flux_jacobi_variable_ee[10][10];
	Real _flux_jacobi_variable_en[10][10];
	RealVectorValue _flux_jacobi_grad_variable_ee[10][10];
	RealVectorValue _flux_jacobi_grad_variable_en[10][10];
	RealVectorValue _lift[10], _lift_old[10];
	RealVectorValue _lift_jacobi_variable_ee[10][10];
	RealVectorValue _lift_jacobi_variable_en[10][10];
	Real _perturbation;
	Real _penalty;

	virtual void precalculateResidual();
	virtual Real computeQpResidual(Moose::DGResidualType type, unsigned int p);
	virtual void precalculateJacobian();
	virtual Real computeQpJacobian(Moose::DGJacobianType type, unsigned int p, unsigned int q);

	void reinit();
	void fluxRiemann();
	void fluxLaxF();
	void fluxHLLCPV();
	void liftOperator();

};
template<>
InputParameters validParams<CFDFaceKernel>();
