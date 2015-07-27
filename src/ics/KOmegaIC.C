
#include "KOmegaIC.h"

template<>
InputParameters validParams<KOmegaIC>()
{
  InputParameters params = validParams<InitialCondition>();
  params += validParams<KOmegaModelBase>();
  params.addRequiredParam<int>("component", "方程分量");
  return params;
}

KOmegaIC::KOmegaIC(const std::string & name, InputParameters parameters) :
		InitialCondition(name, parameters),
		KOmegaModelBase(name, parameters),
	    _component(getParam<int>("component"))
{
}

Real KOmegaIC::value(const Point& p)
{
	Real rho_infty = 1.;
	Real k_infty = 3./2*_tu_infty;
	if(_component == 5)
	{
		return rho_infty*k_infty;
	}
	if(_component == 6)
	{
		return rho_infty*std::log(_reynolds*k_infty/_r_mu);
	}

	return 0.;
}
