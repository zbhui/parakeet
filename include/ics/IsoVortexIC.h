
#pragma once

// MOOSE Includes
#include "CFDInitialCondition.h"

// Forward Declarations
class IsoVortexIC;

template<>
InputParameters validParams<IsoVortexIC>();

/**
 * 等熵涡初始条件
 */
class IsoVortexIC : public CFDInitialCondition
{
public:

  IsoVortexIC(const std::string & name, InputParameters parameters);

private:
//  Point &centre;
protected:
  virtual Real density(const Point &p);
  virtual Real x_momentum(const Point &p);
  virtual Real y_momentum(const Point &p);
  virtual Real z_momentum(const Point &p);
  virtual Real total_energy(const Point &p);
};
