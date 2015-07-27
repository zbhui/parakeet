
#pragma once

#include "BurgersBase.h"
#include "DGKernel.h"

class BurgersFaceKernel;

template<>
InputParameters validParams<BurgersFaceKernel>();

class BurgersFaceKernel :
public DGKernel,
public BurgersBase
{
public:
	BurgersFaceKernel(const std::string &name, InputParameters parameters);

protected:
	virtual Real computeQpResidual(Moose::DGResidualType type);
	virtual Real computeQpJacobian(Moose::DGJacobianType type);

	Real flux(Real ul, Real ur, Point normal);
};
