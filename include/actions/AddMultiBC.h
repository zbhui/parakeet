#pragma once

#include "MooseObjectAction.h"

class AddMultiBC: public MooseObjectAction
{
public:
	AddMultiBC(InputParameters params);

	virtual void act();
};


template<>
InputParameters validParams<AddMultiBC>();
