
#pragma once

#include "MooseObject.h"
#include "MooseVariableBase.h"


class ConservationLaw;
/**
 * CFD base class.
 */
class ConservationLaw
{
public:
	ConservationLaw(const InputParameters & parameters){};
  virtual ~ConservationLaw(){};

protected:
 virtual void inviscousTerm(RealVectorValue &inviscous_term, unsigned component, VariableValue &uh) = 0;
 virtual Real numericalFlux(VariableValue &ul, VariableValue &ur, Point &normal, unsigned component) = 0;

};

