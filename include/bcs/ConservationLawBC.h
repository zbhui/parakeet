
#pragma once

#include "MultiIntegratedBC.h"
#include "BurgersBase.h"

class BurgersBC;

template<>
InputParameters validParams<BurgersBC>();

class BurgersBC :
public IntegratedBC,
public BurgersBase
{
public:
	  BurgersBC(const InputParameters &parameters);
	  virtual ~BurgersBC(){}


protected:
	  virtual Real computeQpResidual();
	  virtual Real computeQpJacobian();

private:
	  Real flux(Real ul, Real ur, Point normal);
//	  Function &_func;
};
