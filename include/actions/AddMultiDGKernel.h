#pragma once

#include "MooseObjectAction.h"

class AddMultiDGKernel : public MooseObjectAction
{
public:
	AddMultiDGKernel(InputParameters params);

	virtual void act();
};


template<>
InputParameters validParams<AddMultiDGKernel>();
