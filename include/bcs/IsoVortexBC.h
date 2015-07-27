
#pragma once

#include "EulerBC.h"

class IsoVortexBC;

template<>
InputParameters validParams<IsoVortexBC>();

class IsoVortexBC :
public EulerBC
{
public:
	IsoVortexBC(const InputParameters & params);
	virtual ~IsoVortexBC(){}


protected:
	virtual void valueAtRightFace(Real *ur);

private:

  Real density(Real t, const Point &p);
  Real x_momentum(Real t, const Point &p);
  Real y_momentum(Real t, const Point &p);
  Real z_momentum(Real t, const Point &p);
  Real total_energy(Real t, const Point &p);
};
