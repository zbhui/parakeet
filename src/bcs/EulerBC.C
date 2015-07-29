
#include "EulerBC.h"

// Euler方程边界条件
template<>
InputParameters validParams<EulerBC>()
{
	InputParameters params = validParams<CFDBC>();
	params += validParams<CFDBase>();

	return params;
}

EulerBC::EulerBC(const InputParameters & parameters):
		CFDBC(parameters),
		CFDBase(parameters)
{
}

void EulerBC::precalculateResidual()
{
	valueAtLeftFace(_ul);
	valueAtRightFace(_ur);
	Point normal = _normals[_qp];
	fluxRiemann(_flux, _ul, _ur, normal);
}

Real EulerBC::computeQpResidual(unsigned int p)
{
	return _flux[p] * _test[_i][_qp];
}

Real EulerBC::computeQpJacobian(unsigned int p, unsigned int q)
{
	return _jacobi_variable[p][q]*_phi[_j][_qp]*_test[_i][_qp];
}

void EulerBC::fluxRiemann(Real* flux, Real* ul, Real* ur, Point& normal)
{
	RealVectorValue fl[5], fr[5];
	inviscousTerm(fl, ul);
	inviscousTerm(fr, ur);

	Real rho, u, v, w, pre;
	rho = (ul[0] + ur[0])/2.;
	u = (ul[1] + ur [1])/rho/2;
	v = (ul[2] + ur [2])/rho/2;
	w = (ul[3] + ur [3])/rho/2;
	pre = (pressure(ul) + pressure(ur))/2.;
	Real lam = fabs(u*normal(0) + v * normal(1) + w * normal(2)) + sqrt(_gamma*pre/rho);
	for (int eq = 0; eq < 5; ++eq)
	{
		flux[eq] = 0.5*(fl[eq] + fr[eq])*normal + lam*(ul[eq] - ur[eq]);
	}
}

void EulerBC::precalculateJacobian()
{
	Real flux_new[10], flux[10];
	Real ul[10], ur[10];

	Point normal = _normals[_qp];
	valueAtLeftFace(_ul);
	valueAtRightFace(_ur);
	fluxRiemann(flux, _ul, _ur, normal);
	for (int q = 0; q < _n_equation; ++q)
	{
		_ul[q] += _ds;
		valueAtRightFace(ur);

		fluxRiemann(flux_new, _ul, ur, normal);
		for (int p = 0; p < _n_equation; ++p)
			_jacobi_variable[p][q] = (flux_new[p] - flux[p])/_ds;

		_ul[q] -= _ds;
	}
}

void EulerBC::wallBC(Real* ur)
{
    Real pre = pressure(_ul);
	const Point &normal = _normals[_qp];
    Real  vn = _ul[1]*normal(0) + _ul[2]*normal(1) + _ul[3]*normal(2);

    ur[0] = _ul[0];
    ur[1] = _ul[1] - 2.0 * vn * normal(0);
    ur[2] = _ul[2] - 2.0 * vn * normal(1);
    ur[3] = _ul[3] - 2.0 * vn * normal(2);
    ur[4] = pre/(_gamma-1) + 0.5*( ur[1]*ur[1] + ur[2]*ur[2] + ur[3]*ur[3] )/ur[0];
}

void EulerBC::farFieldBC(Real* ur)
{
	const Point &normal = _normals[_qp];

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

	rhoL = _ul[0];
	uL = _ul[1] / rhoL;
	vL = _ul[2] / rhoL;
	wL = _ul[3] / rhoL;
	vel = sqrt(uL * uL + vL * vL + wL * wL);
	pL = pressure(_ul);
	cL = sqrt(fabs(_gamma * pL / rhoL));
	vnL =  normal(0) * uL + normal(1) * vL + normal(2) * wL;

	if (vel >= cL) {	//超声速
		if (vnL >= 0.0) //exit
		{
			ur[0] = _ul[0];
			ur[1] = _ul[1];
			ur[2] = _ul[2];
			ur[3] = _ul[3];
			ur[4] = _ul[4];
		}
		else //inlet
		{
			ur[0] = rhoR;
			ur[1] = rhoR * uR;
			ur[2] = rhoR * vR;
			ur[3] = rhoR * wR;
			ur[4] = pR / (_gamma - 1) + 0.5 * rhoR * (uR * uR + vR * vR + wR * wR);
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
			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1) + 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
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
			ur[4] = cb * cb * ur[0] / _gamma / (_gamma - 1) + 0.5 * (ur[1] * ur[1] + ur[2] * ur[2] + ur[3] * ur[3]) / ur[0];
		}
	}
}

void EulerBC::symmetricBC(Real* ur)
{
	wallBC(ur);
}
