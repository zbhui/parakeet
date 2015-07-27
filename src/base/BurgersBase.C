
#include "BurgersBase.h"

template<>
InputParameters validParams<BurgersBase>()
{
  InputParameters params = validParams<MooseObject>();

  return params;
}

BurgersBase::BurgersBase(const InputParameters &parameters)
{
}

void BurgersBase::inviscousTerm(RealVectorValue &flux_term, Real uh)
{
	Real flux = uh * uh/2;
	flux = uh;
	flux_term(0) = flux;
	flux_term(1) = flux;
	flux_term(2) = flux;
}




