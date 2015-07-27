/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "SodIC.h"

template<>
InputParameters validParams<SodIC>()
{
  InputParameters params = validParams<CFDInitialCondition>();
  return params;
}

SodIC::SodIC(const std::string & name, InputParameters parameters) :
    CFDInitialCondition(name, parameters)
{}

Real SodIC::density(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	if(x < 0.5) return 1;
	else return 0.125;
}

Real SodIC::x_momentum(const Point &p)
{
	return 0.0;
}

Real SodIC::y_momentum(const Point &p)
{
	return 0.0;
}

Real SodIC::z_momentum(const Point &p)
{
	return 0.0;
}

Real SodIC::total_energy(const Point &p)
{
	Real pre = 0.;
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
	if(x < 0.5) pre = 1;
	else pre = 0.1;

	return pre/(_gamma-1);
}
