
#include "KOmegaBC.h"

template<>
InputParameters validParams<KOmegaBC>()
{
	InputParameters params = validParams<CFDBC>();
	params += validParams<KOmegaModelBase>();
//	params.addParam<Real>("tu_infty", 0.005, "远场湍流度");
//	params.addParam<Real>("r_mu", 1e-05, "mu_infty");
	return params;
}

KOmegaBC::KOmegaBC(const std::string &name, InputParameters parameters):
		CFDBC(name, parameters),
		KOmegaModelBase(name, parameters)
//		_tu_infty(getParam<Real>("tu_infty")),
//		_r_mu(getParam<Real>("r_mu"))
{

}

void KOmegaBC::precalculateResidual()
{
	const unsigned int elem_b_order = static_cast<unsigned int> (_var.getOrder());
	const double h_elem = (_current_elem_volume+_current_elem_volume)/_current_side_volume * 1./std::pow(elem_b_order, 2.)/2.;

	Real ul[10], ur[10], uh_bar[10];
	RealGradient dul[10], dur[10], duh[10];
	RealVectorValue vfl[10], vfr[10];

	mooseAssert(_n_equation < 10, "multiBC方程个数应<10");
	mooseAssert(_qrule->n_points() < 40, "mulitBC积分点个数应<40");

	for (_qp = 0; _qp < _qrule->n_points(); _qp++)
	{
		valueAtLeftFace(ul);
		valueAtRightFace(ur);
		valueGradAtLeftFace(dul);
		valueGradAtRightFace(dur);

		Point normal = _normals[_qp];

		viscousTerm(vfl, ul, dul);
		viscousTerm(vfr, ur, dur);
		fluxRiemann(_flux[_qp], ul, ur, normal);

		for (int eq = 0; eq < _n_equation; ++eq)
		{
			uh_bar[_eq] = (ul[_eq]+ur[_eq])/2.;
			duh[_eq] = (ul[_eq]-ur[_eq])*normal;
		}

		viscousTerm(_penalty[_qp], ul, duh);
		viscousTerm(_penalty_neighbor[_qp], ur, duh);

		for (_eq = 0; _eq < _n_equation; ++_eq)
		{
			_flux[_qp][_eq] -= ((vfl[_eq]+vfr[_eq]) - _sigma/h_elem*(_penalty[_qp][_eq] + _penalty_neighbor[_qp][_eq]))/2.*normal;
		}
	}
}

Real KOmegaBC::computeQpResidual()
{
	return _flux[_qp][_eq] * _test[_i][_qp] + 0.5*_epsilon * _penalty[_qp][_eq]* _grad_test[_i][_qp];
}

void KOmegaBC::wallBC(Real* ur)
{
	Real rho = ur[0];
	Real mu = physicalViscosity(ur);
	Real distance = _current_elem_volume/_current_side_volume;

	Real ui[10];
	valueAtLeftFace(ui);
    Real pre = pressure(ui);
    Real twall = 1.;

    ur[0] = _gamma*_mach*_mach*pre/twall;
    ur[1] = 0.;
    ur[2] = 0.;
    ur[3] = 0.;
    ur[4] = pre/(_gamma-1) + 0.5*( ur[1]*ur[1] + ur[2]*ur[2] + ur[3]*ur[3] )/ur[0];
	ur[5] = rho * 0;
	ur[6] = rho * std::log(60.*mu/(_reynolds*rho*_beta_o*distance));
}

void KOmegaBC::farFieldBC(Real* ur)
{
	Real ui[10];
	valueAtLeftFace(ui);
	Point normal = _normals[_qp];

	Real rhoR, uR, vR, wR, pR;
	Real rhoL, uL, vL, wL, pL;
	Real cR, cL, cb;
	Real vnR, vnL, vnb;
	Real vel, s;
	Real Rp, Rm;

	uR = 1.0 * cos(_attack) * cos(_slide);
	vR = 1.0 * sin(_attack) * cos(_slide);
	wR = 1.0 * sin(_slide);
	rhoR = 1.0;

	pR = 1 / _gamma /_mach / _mach;
	cR = sqrt(fabs(_gamma * pR / rhoR));
	vnR = normal(0) * uR + normal(1) * vR + normal(2) * wR;

	rhoL = ui[0];
	uL = ui[1] / rhoL;
	vL = ui[2] / rhoL;
	wL = ui[3] / rhoL;
	vel = sqrt(uL * uL + vL * vL + wL * wL);
	pL = pressure(ui);
	cL = sqrt(fabs(_gamma * pL / rhoL));
	vnL =  normal(0) * uL + normal(1) * vL + normal(2) * wL;

	if (sqrt(vel) >= cL) {	//超声速
		if (vnL >= 0.0) //exit
		{
			ur[0] = ui[0];
			ur[1] = ui[1];
			ur[2] = ui[2];
			ur[3] = ui[3];
			ur[4] = ui[4];
		}
		else //inlet
		{
			ur[0] = rhoR;
			ur[1] = rhoR * uR;
			ur[2] = rhoR * vR;
			ur[3] = rhoR * wR;
			ur[4] = pR / (_gamma - 1)
					+ 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
		}
	}
	else
	{ 	//  亚声速
		if (vnL >= 0.0)
		{			//exit
			s = pL / pow(rhoL, _gamma);
			Rp = vnL + 2 * cL / (_gamma - 1);
			Rm = vnR - 2 * cR / (_gamma - 1);
			vnb = (Rp + Rm) / 2.0;
			cb = (Rp - Rm) * (_gamma - 1) / 4.0;

			ur[0] = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
			ur[1] = ur[0] * (uL + normal(0) * (vnb - vnL));
			ur[2] = ur[0] * (vL + normal(1) * (vnb - vnL));
			ur[3] = ur[0] * (wL + normal(2) * (vnb - vnL));
			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1)
					+ 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
		}
		else
		{
			s = pR / pow(rhoR, _gamma);
			Rp = -vnR + 2.0 * cR / (_gamma - 1);
			Rm = -vnL - 2.0 * cL / (_gamma - 1);
			vnb = -(Rp + Rm) / 2.0;
			cb = (Rp - Rm) * (_gamma - 1) / 4.0;

			ur[0] = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
			ur[1] = ur[0] * (uR + normal(0) * (vnb - vnR));
			ur[2] = ur[0] * (vR + normal(1) * (vnb - vnR));
			ur[3] = ur[0] * (wR + normal(2) * (vnb - vnR));
			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1)
					+ 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
		}
	}

	Real k_infty = 3./2*_tu_infty;
	ur[5] = rhoR*k_infty;
	ur[6] = rhoR*std::log(_reynolds*k_infty/_r_mu);
}

void KOmegaBC::symmetricBC(Real* ur)
{
	Real ui[10];
	valueAtLeftFace(ui);
    Real pre = pressure(ui);
	Point normal = _normals[_qp];
    Real  vn = ui[1]*normal(0) + ui[2]*normal(1) + ui[3]*normal(2);

    ur[0] = ui[0];
    ur[1] = ui[1] - 2.0 * vn * normal(0);
    ur[2] = ui[2] - 2.0 * vn * normal(1);
    ur[3] = ui[3] - 2.0 * vn * normal(2);
    ur[4] = pre/(_gamma-1) + 0.5*( ur[1]*ur[1] + ur[2]*ur[2] + ur[3]*ur[3] )/ur[0];
	ur[5] = ui[5];
	ur[6] = ui[6];
}

void KOmegaBC::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
{
	RealVectorValue fl[10], fr[10];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real rho, u, v, w, pre;
	rho = (ul[0] + ur[0])/2.;
	u = (ul[1] + ur[1])/rho/2.;
	v = (ul[2] + ur[2])/rho/2.;
	w = (ul[3] + ur[3])/rho/2.;
	pre = (pressure(ul) + pressure(ur))/2.;
	Real lam = fabs(u*normal(0) + v * normal(1) + w * normal(2)) + std::sqrt(_gamma*pre/rho);
	for (int eq = 0; eq < _n_equation; ++eq)
	{
		flux[eq] = 0.5*((fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]));
	}
}

Real KOmegaBC::computeQpJacobian()
{
	return 0.;
}

Real KOmegaBC::computeQpOffDiagJacobian(unsigned int jvar)
{
	return 0.;
}
