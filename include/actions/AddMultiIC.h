
#pragma once

#include "MooseObjectAction.h"

class AddMultiIC : public MooseObjectAction
{
public:
	AddMultiIC(InputParameters params);

	virtual void act();
};

template<>
InputParameters validParams<AddMultiIC>();
