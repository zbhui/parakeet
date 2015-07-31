#pragma once

#include "libmesh/tensor_value.h"
#include "MooseObject.h"

class CFDDataPack
{
public:
	CFDDataPack(Real mach, Real reynolds, Real gamma = 1.4, Real prandtl = 0.72);

public:
	Real _mach, _reynolds, _gamma, _prandtl;
	Real r, p, t, h, s, c, m, q, re;
	RealVectorValue vel, mom;
	Real vel_size, vel_div, mom_size;

	RealGradient grad_rho, grad_enthalpy;
	RealTensor grad_mom, grad_vel, tau;

	Real vis;
	RealVectorValue invis_flux[10], vis_flux[10], flux[10];
	Real uh[10];
	RealGradient duh[10];

	void reinit();
//	void reinit(CFDProblem &cfd_problem);

//	virtual void InvisFlux(RealVectorValue* inviscous_term){};
};
