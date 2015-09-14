#pragma once

#include "libmesh/tensor_value.h"
#include "MooseObject.h"

class CFDProblem;

class CFDDataPack
{
public:
	CFDDataPack(CFDProblem &cfd_problem);

public:
	CFDProblem &_cfd_problem;
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
	void reinitInviscous();
	void reinitViscous();
};
