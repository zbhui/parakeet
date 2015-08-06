
#pragma once

#include "Action.h"
#include "MooseObjectAction.h"

class CLawICAction : public MooseObjectAction
{
public:
  CLawICAction(InputParameters params);
  virtual void act();

protected:
};

template<>
InputParameters validParams<CLawICAction>();
