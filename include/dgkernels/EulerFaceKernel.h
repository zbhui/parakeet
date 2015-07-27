
#pragma once

#include "MultiDGKernel.h"
#include "CFDBase.h"

class EulerFaceKernel;

template<>
InputParameters validParams<EulerFaceKernel>();

class EulerFaceKernel :
public MultiDGKernel,
public CFDBase
{
public:
	EulerFaceKernel(const std::string &name, InputParameters parameters);
	virtual ~EulerFaceKernel(){}

protected:
	Real _flux[40][10];

	Real _jacobi_variable_ee[40][10][10];
	Real _jacobi_variable_en[40][10][10];
	Real _jacobi_variable_ne[40][10][10];
	Real _jacobi_variable_nn[40][10][10];

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	virtual void precalculateResidual();
	virtual void precalculateJacobian();

};
