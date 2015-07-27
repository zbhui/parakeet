
#pragma once

// MOOSE Includes
#include "CFDInitialCondition.h"

// Forward Declarations
class CFDPassFlowIC;

template<>
InputParameters validParams<CFDPassFlowIC>();

/**
 * 等熵涡初始条件
 */
class CFDPassFlowIC : public CFDInitialCondition
{
public:

  CFDPassFlowIC(const std::string & name, InputParameters parameters);

private:
  Real _velocity;
protected:
  virtual Real density(const Point &p);
  virtual Real x_momentum(const Point &p);
  virtual Real y_momentum(const Point &p);
  virtual Real z_momentum(const Point &p);
  virtual Real total_energy(const Point &p);
};
