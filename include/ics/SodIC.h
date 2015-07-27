
#pragma once

// MOOSE Includes
#include "CFDInitialCondition.h"

// Forward Declarations
class SodIC;

template<>
InputParameters validParams<SodIC>();

/**
 * Sod问题初始条件
 */
class SodIC :
public CFDInitialCondition
{
public:

  SodIC(const InputParameters & parameters);

protected:
  virtual Real density(const Point &p);
  virtual Real x_momentum(const Point &p);
  virtual Real y_momentum(const Point &p);
  virtual Real z_momentum(const Point &p);
  virtual Real total_energy(const Point &p);
};
