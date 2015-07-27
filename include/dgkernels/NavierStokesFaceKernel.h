#pragma once

#include "EulerFaceKernel.h"

class NavierStokesFaceKernel;

template<>
InputParameters validParams<NavierStokesFaceKernel>();

class NavierStokesFaceKernel :
public EulerFaceKernel
{
public:
	NavierStokesFaceKernel(const InputParameters & parameters);
	virtual ~NavierStokesFaceKernel(){}

protected:
	RealVectorValue _penalty[40][10];
	RealVectorValue _penalty_neighbor[40][10];

	RealVectorValue _jacobi_gradient_ee[40][10][10];
	RealVectorValue _jacobi_gradient_en[40][10][10];
	RealVectorValue _jacobi_gradient_ne[40][10][10];
	RealVectorValue _jacobi_gradient_nn[40][10][10];

	RealVectorValue _jacobi_penalty_variable_ee[40][10][10];
	RealVectorValue _jacobi_penalty_variable_en[40][10][10];
	RealVectorValue _jacobi_penalty_variable_ne[40][10][10];
	RealVectorValue _jacobi_penalty_variable_nn[40][10][10];

	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);
	virtual Real computeQpOffDiagJacobian(Moose::DGJacobianType type, unsigned int jvar);

	virtual void precalculateResidual();

	virtual void fluxTerm(Real *flux, Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	virtual void penaltyTerm(RealVectorValue *penalty, Real *ul, Real *ur);
	virtual void penaltyTermNeighbor(RealVectorValue *penalty_neighbor, Real *ul, Real *ur);

	virtual void computeJacobianVariable(Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);
	virtual void computeJacobianGradient(Real *ul, Real *ur, RealGradient *dul, RealGradient *dur);

};





