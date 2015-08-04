
#pragma once

#include "Action.h"

class CFDAction;

template<>
InputParameters validParams<CFDAction>();

class CFDAction : public Action
{
public:
  CFDAction(const InputParameters &params);

  virtual void act();
};

