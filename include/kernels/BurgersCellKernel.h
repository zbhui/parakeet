#pragma once

#include "MultiKernel.h"
#include "FDKernel.h"
#include "BurgersBase.h"

// 前置声明
class BurgersCellKernel;

template<>
InputParameters validParams<BurgersCellKernel>();

class BurgersCellKernel :
public FDKernel,
public BurgersBase
{
public:
	BurgersCellKernel(const std::string & name, InputParameters parameters);
	virtual ~BurgersCellKernel(){}

protected:
	virtual Real computeQpResidual();
	virtual Real computeQpJacobian();
};
