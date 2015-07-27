
#pragma once

#include "InitialCondition.h"
#include "BurgersBase.h"

// Forward Declarations
class SinIC;

template<>
InputParameters validParams<SinIC>();

/**
 *
 */
class SinIC :
public InitialCondition,
public BurgersBase
{
public:

  /**
   * Constructor: Same as the rest of the MOOSE Objects
   */
	SinIC(const std::string & name, InputParameters parameters);

  virtual Real value(const Point & p);

};
