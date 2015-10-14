
#include "CFDDataPack.h"
#include "CFDProblem.h"

CFDDataPack::CFDDataPack(CFDProblem &cfd_problem):
	_cfd_problem(cfd_problem),
	_mach(_cfd_problem._mach),
	_reynolds(_cfd_problem._reynolds),
	_gamma(_cfd_problem._gamma),
	_prandtl(_cfd_problem._prandtl)
{

}

void CFDDataPack::reinit()
{
	reinitInviscous();
	reinitViscous();
}

void CFDDataPack::reinitInviscous()
{
	r = uh[0];
	mom(0) = uh[1]; mom(1) = uh[2]; mom(2) = uh[3]; mom_size = mom.size();
	re = uh[4];
	vel = mom/r; vel_size = vel.size();

	p = (_gamma-1) * (fabs(uh[4]) - 0.5*(uh[1]*uh[1] + uh[2]*uh[2] + uh[3]*uh[3])/uh[0]);
	t = _gamma*_mach*_mach*p/uh[0];
	h = (fabs(uh[4]) + p)/uh[0];
	s = p/pow(r, _gamma);
	c = sqrt(fabs(_gamma*p/r));
	m = vel_size/c;

	invis_flux[0] = r*vel;

	invis_flux[1] = mom(0)*vel;
	invis_flux[1](0) += p;

	invis_flux[2] = mom(1)*vel;
	invis_flux[2](1) += p;

	invis_flux[3] = mom(2)*vel;
	invis_flux[3](2) += p;

	invis_flux[4] = r*h*vel;
}

void CFDDataPack::reinitViscous()
{
	grad_rho = duh[0];
	grad_mom = RealTensor(duh[1], duh[2], duh[3]);

	for (int a = 0; a < 3; ++a)
		for (int b = 0; b < 3; ++b)
			grad_vel(a, b) = (grad_mom(a, b) - vel(a)*grad_rho(b))/r;

	tau = grad_vel + grad_vel.transpose();
	vel_div =  grad_vel.tr();
	tau(0, 0) -= 2./3*vel_div; tau(1, 1) -= 2./3*vel_div; tau(2, 2) -= 2./3*vel_div;

	if(_cfd_problem._vis_type == 0) vis = 0.0;
	else if(_cfd_problem._vis_type == 1) vis = 1.0;
	else if(_cfd_problem._vis_type == 2)
	{
		Real t_ref(288), t_s(110.4);
		vis = pow(t/t_ref, 1.5)*(t_ref+t_s)/(t+t_s);
	}
	else mooseError("不可知的粘性模型");
	tau *= vis/_reynolds;

	grad_enthalpy = (duh[4]-uh[4]/uh[0] * duh[0])/r - grad_vel.transpose() * vel;
	grad_enthalpy *= (vis/_reynolds)*(_gamma/_prandtl);

	invis_flux[0] = r*vel;

	invis_flux[1] = mom(0)*vel;
	invis_flux[1](0) += p;

	invis_flux[2] = mom(1)*vel;
	invis_flux[2](1) += p;

	invis_flux[3] = mom(2)*vel;
	invis_flux[3](2) += p;

	invis_flux[4] = r*h*vel;

	vis_flux[0].zero();
	vis_flux[1] = tau.row(0);
	vis_flux[2] = tau.row(1);
	vis_flux[3] = tau.row(2);
	vis_flux[4] = tau * vel + grad_enthalpy;
}
