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

#include "SinIC.h"

template<>
InputParameters validParams<SinIC>()
{
  InputParameters params = validParams<InitialCondition>();
  return params;
}

SinIC::SinIC(const std::string & name, InputParameters parameters) :
    InitialCondition(name, parameters),
    BurgersBase(name, parameters)
{}

Real SinIC::value(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real PI = libMesh::pi;
	return sin(2*PI*(x+y+z));
}

