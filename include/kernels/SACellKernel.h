
#pragma once

#include "SAModelBase.h"
#include "MultiKernel.h"

class SACellKernel;

template<>
InputParameters validParams<SACellKernel>();

class SACellKernel :
public Kernel,
public SAModelBase
{
public:
	SACellKernel(const std::string & name, InputParameters parameters);
	virtual ~SACellKernel(){}

protected:
	  virtual Real computeQpResidual();
	  virtual Real computeQpJacobian();
	  virtual Real computeQpOffDiagJacobian(unsigned int jvar);
};
