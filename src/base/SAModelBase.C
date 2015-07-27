
#include "SAModelBase.h"

template<>
InputParameters validParams<SAModelBase>()
{
	InputParameters params = validParams<CFDBase>();

	return params;
}

SAModelBase::SAModelBase(const std::string& name, InputParameters parameters):
		CFDBase(name, parameters),
		_cb1(0.1355), _cb2(0.622), _sigma(2./3), _kappa(0.41),
		_cw2(0.3), _cw3(2.0), _cv1(7.1),
		_ct1(1.0), _ct2(2.0), _ct3(1.2), _ct4(0.5),
		_prandtl_turb(0.9)

{
	_cw1 = _cb1/_kappa/_kappa + (1+_cb2)/_sigma;
}
