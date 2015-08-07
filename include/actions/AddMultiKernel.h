#pragma once

#include "MooseObjectAction.h"

class AddMultiKernel : public MooseObjectAction
{
public:
	AddMultiKernel(InputParameters params);

	virtual void act();
};

template<>
InputParameters validParams<AddMultiKernel>();
