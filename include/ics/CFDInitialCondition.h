
#pragma once

#include "InitialCondition.h"
#include "CFDBase.h"

// Forward Declarations
class CFDInitialCondition;

template<>
InputParameters validParams<CFDInitialCondition>();

/**
 * 为CFD提供初始条件接口
 */
class CFDInitialCondition :
public InitialCondition,
public CFDBase
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
	CFDInitialCondition(const InputParameters & parameters);

  /**
   * The value of the variable at a point.
   *
   * This must be overriden by derived classes.
   */
  virtual Real value(const Point & p);

protected:
  virtual Real density(const Point &p) = 0;
  virtual Real x_momentum(const Point &p) = 0;
  virtual Real y_momentum(const Point &p) = 0;
  virtual Real z_momentum(const Point &p) = 0;
  virtual Real total_energy(const Point &p) = 0;

//  unsigned int _component;
};
