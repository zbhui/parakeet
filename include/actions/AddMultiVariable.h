#pragma once

#include "AddVariableAction.h"

class AddMultiVariable : 	public AddVariableAction
{
public:
	AddMultiVariable(InputParameters params);

	virtual void act();
};

template<>
InputParameters validParams<AddMultiVariable>();
