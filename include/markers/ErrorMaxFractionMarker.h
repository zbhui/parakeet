
#pragma once

#include "IndicatorMarker.h"

class ErrorMaxFractionMarker : public IndicatorMarker
{
public:
	ErrorMaxFractionMarker(const InputParameters & parameters);
  virtual ~ErrorMaxFractionMarker(){};

  virtual void markerSetup();

protected:
  virtual MarkerValue computeElementMarker();

  Real _coarsen;
  Real _refine;

  Real _max;
  Real _min;
  Real _delta;
  Real _refine_cutoff;
  Real _coarsen_cutoff;
};


template<>
InputParameters validParams<ErrorMaxFractionMarker>();

