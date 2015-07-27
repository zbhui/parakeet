
#pragma once

#include "MooseObject.h"
#include "MooseVariableBase.h"

class BurgersBase;

template<>
InputParameters validParams<BurgersBase>();


class BurgersBase
{
public:
	BurgersBase(const InputParameters & parameters);
  virtual ~BurgersBase(){};

protected:
  virtual void inviscousTerm(RealVectorValue &flux_term, Real uh);
};

