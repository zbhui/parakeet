#pragma once

#include "Function.h"

class IsoVortexExact;

template<>
InputParameters validParams<IsoVortexExact>();

class IsoVortexExact : public Function
{
public:
  IsoVortexExact(const InputParameters & parameters);

  Real value(Real t, const Point & p);

protected:
};

