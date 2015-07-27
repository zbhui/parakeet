
#pragma once
#include "CFDBase.h"

class SAModelBase;

template<>
InputParameters validParams<SAModelBase>();

/**
 * SA湍流模型 base class.
 */
class SAModelBase :
public CFDBase
{
public:
	SAModelBase(const InputParameters & parameters);
	virtual ~SAModelBase(){};

protected:
	Real _cb1, _cb2, _sigma, _kappa;
	Real _cw1, _cw2, _cw3, _cv1;
	Real _ct1, _ct2, _ct3, _ct4;
	Real _prandtl_turb;
};


