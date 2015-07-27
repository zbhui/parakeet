
#pragma once
#include "CFDBase.h"

class KOmegaModelBase;

template<>
InputParameters validParams<KOmegaModelBase>();

/**
 * K-Oemga湍流模型 base class.
 */
class KOmegaModelBase :
public CFDBase
{
public:
	KOmegaModelBase(const std::string & name, InputParameters parameters);
	virtual ~KOmegaModelBase(){};

protected:
	Real eddyViscosity(Real *uh);

	virtual void inviscousTerm(RealVectorValue *inviscous_term, Real *uh);
	virtual void viscousTerm(RealVectorValue *viscous_term, Real *uh, RealGradient *duh);
	virtual void sourceTerm(Real *source_term, Real *uh, RealGradient *duh);

protected:
	Real _sigma_k, _sigma_o, _beta_k, _beta_o, _alpha_o;
	Real _prandtl_turb;

	Real _tu_infty;
	Real _r_mu;
};


