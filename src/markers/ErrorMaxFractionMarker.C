
#include "ErrorMaxFractionMarker.h"

template<>
InputParameters validParams<ErrorMaxFractionMarker>()
{
  InputParameters params = validParams<IndicatorMarker>();
  params.addRequiredParam<Real>("coarsen", "Elements with error less than this will be coarsened.");
  params.addRequiredParam<Real>("refine", "Elements with error more than this will be refined.");
  return params;
}


ErrorMaxFractionMarker::ErrorMaxFractionMarker(const InputParameters & parameters) :
    IndicatorMarker(parameters),
    _coarsen(parameters.get<Real>("coarsen")),
    _refine(parameters.get<Real>("refine"))
{
}

void ErrorMaxFractionMarker::markerSetup()
{
  _min = std::numeric_limits<Real>::max();
  _max = 0;

  for (unsigned int i=0; i<_error_vector.size(); i++)
  {
    _min = std::min(_min, static_cast<Real>(_error_vector[i]));
    _max = std::max(_max, static_cast<Real>(_error_vector[i]));
  }

  _delta = _max-_min;
  _refine_cutoff = (_refine)*_max;
  _coarsen_cutoff = _coarsen*_delta + _min;
}

Marker::MarkerValue
ErrorMaxFractionMarker::computeElementMarker()
{
  Real error = _error_vector[_current_elem->id()];

  if (error > _refine_cutoff)
    return REFINE;
  else if (error < _coarsen_cutoff)
    return COARSEN;

  return DO_NOTHING;
}

