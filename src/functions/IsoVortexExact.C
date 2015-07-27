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

#include "IsoVortexExact.h"

template<>
InputParameters validParams<IsoVortexExact>()
{
  InputParameters params = validParams<Function>();
//  params += validParams<IsoVortexIC>();
  return params;
}

IsoVortexExact::IsoVortexExact(const InputParameters & parameters) :
    Function(parameters)
//    IsoVortexIC(name, parameters)
{}

Real
IsoVortexExact::value(Real t, const Point & p)
{
	Real x = p(0)-t;
	Real y = p(1)-t;
	Real z = p(2)-t;
//
//	Real pi = 3.1415926535;
//	return std::sin(2*pi*(x+y+z));

	Real gam = 1.4, gamm1 = gam - 1, epi = 5.0;
	Real xb, yb, r2;
	Real rho,T;

	Real PI = libMesh::pi;
	xb = x+5;
	yb = y+5;
	r2 = xb * xb + yb * yb;

	T = 1.0 - gamm1 * epi * epi / ( 8 * gam * PI* PI ) * exp( 1 - r2 );
	rho = pow( T, 1 / gamm1 );
	return rho;
}
