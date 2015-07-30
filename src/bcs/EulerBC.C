
#include "EulerBC.h"

template<>
InputParameters validParams<EulerBC>()
{
	InputParameters params = validParams<CFDBC>();
	params += validParams<CFDBase>();

	return params;
}

EulerBC::EulerBC(const InputParameters & parameters):
		CFDBC(parameters),
		CFDBase(parameters),
		_cfd_data(_mach, _reynolds, _gamma, _prandtl),
		_cfd_data_neighbor(_mach, _reynolds, _gamma, _prandtl)
{
}

void EulerBC::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
	}
	boundaryCondition();

	_cfd_data.reinit();
	_cfd_data_neighbor.reinit();

	fluxRiemann();
}

Real EulerBC::computeQpResidual(unsigned int p)
{
	return _flux[p] * _test[_i][_qp];
}

Real EulerBC::computeQpJacobian(unsigned int p, unsigned int q)
{
	return _jacobi_variable[p][q]*_phi[_j][_qp]*_test[_i][_qp];
}

void EulerBC::fluxRiemann()
{
	Real lam = fabs(_cfd_data.vel*_normals[_qp]) + _cfd_data.c;
	lam += fabs(_cfd_data_neighbor.vel*_normals[_qp]) + _cfd_data_neighbor.c;
	lam /= 2.;
	for (int p = 0; p < _n_equation; ++p)
	{
		_flux[p] = 0.5*(_cfd_data.invis_flux[p] + _cfd_data_neighbor.invis_flux[p])*_normals[_qp] +
				    lam*(_cfd_data.uh[p] - _cfd_data_neighbor.uh[p]);
	}
}

void EulerBC::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
		_flux_old[q] = _flux[q];

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _ds;
		boundaryCondition();
		_cfd_data.reinit();
		_cfd_data_neighbor.reinit();
		fluxRiemann();
		for (int p = 0; p < _n_equation; ++p)
			_jacobi_variable[p][q] = (_flux[p] - _flux_old[p])/_ds;

		_cfd_data.uh[q] -= _ds;
	}
//	_cfd_data.reinit();
//	_cfd_data_neighbor.reinit();
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

void EulerBC::boundaryCondition()
{
}

void EulerBC::symmetricBC(Real* ur)
{
	wallBC(ur);
}

void EulerBC::wallBC()
{
	const Point &normal = _normals[_qp];
    Real  vn = _cfd_data.mom*normal;

    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
    _cfd_data_neighbor.uh[1] = _cfd_data.uh[1] - 2*vn*normal(0);
    _cfd_data_neighbor.uh[2] = _cfd_data.uh[2] - 2*vn*normal(1);
    _cfd_data_neighbor.uh[3] = _cfd_data.uh[3] - 2*vn*normal(2);
    _cfd_data_neighbor.uh[4] = _cfd_data.uh[4];
}

void EulerBC::farFieldBC()
{
}

void EulerBC::symmetricBC()
{
}
