
#include "EulerBC.h"
#include "Eigen/Eigen"

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

void EulerBC::boundaryCondition()
{
	if(_bc_type == "wall")
	{
		wall();
		return;
	}
	if(_bc_type == "far-field")
	{
		farfield();
		return;
	}
	if(_bc_type == "symmetric")
	{
		symmetric();
		return;
	}
	mooseError("未定义的边界条件类型");
}

void EulerBC::wall()
{
	const Point &normal = _normals[_qp];
    Real  vn = _cfd_data.mom*normal;

    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
    _cfd_data_neighbor.uh[1] = _cfd_data.uh[1] - 2*vn*normal(0);
    _cfd_data_neighbor.uh[2] = _cfd_data.uh[2] - 2*vn*normal(1);
    _cfd_data_neighbor.uh[3] = _cfd_data.uh[3] - 2*vn*normal(2);
    _cfd_data_neighbor.uh[4] = _cfd_data.uh[4];
}

void EulerBC::farfield()
{
	Real rho_inf = 1.0;
	Eigen::Vector3d vel_inf(1, 0, 0);
	Real p_inf = 1/_gamma/_mach/_mach;

	CFDDataPack cfd_data_infity(_mach, _reynolds, _gamma, _prandtl);
	cfd_data_infity.r = rho_inf;
	cfd_data_infity.mom(0) = 1.;
	cfd_data_infity.mom(1) = 0.;
	cfd_data_infity.mom(2) = 0.;
	cfd_data_infity.re = 1/_gamma/_mach/_mach + 0.5*cfd_data_infity.mom.size()/1.;

	Real vnl = _cfd_data.vel*_normals[_qp];
	Real vnr = cfd_data_infity.vel*_normals[_qp];
	if(vnr <= 0)  // 入口
	{
		if (_cfd_data.m > 1) //超音速
		{
		    _cfd_data_neighbor.uh[0] = cfd_data_infity.r;
		    _cfd_data_neighbor.uh[1] = cfd_data_infity.mom(0);
		    _cfd_data_neighbor.uh[2] = cfd_data_infity.mom(1);
		    _cfd_data_neighbor.uh[3] = cfd_data_infity.mom(2);
		    _cfd_data_neighbor.uh[4] = cfd_data_infity.re;
		}
		else	//亚音速
		{
		    _cfd_data_neighbor.uh[0] = cfd_data_infity.r;
		    _cfd_data_neighbor.uh[1] = cfd_data_infity.mom(0);
		    _cfd_data_neighbor.uh[2] = cfd_data_infity.mom(1);
		    _cfd_data_neighbor.uh[3] = cfd_data_infity.mom(2);
		    _cfd_data_neighbor.uh[4] = _cfd_data.p/(_gamma-1) + 0.5 * rho_inf*vel_inf.squaredNorm();
		}
	}
	else  //出口
	{
		if (_cfd_data.m > 1) //超音速
		{
		    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
		    _cfd_data_neighbor.uh[1] = _cfd_data.uh[1];
		    _cfd_data_neighbor.uh[2] = _cfd_data.uh[2];
		    _cfd_data_neighbor.uh[3] = _cfd_data.uh[3];
		    _cfd_data_neighbor.uh[4] = _cfd_data.uh[4];
		}
		else	//亚音速
		{
		    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
		    _cfd_data_neighbor.uh[1] = _cfd_data.uh[1];
		    _cfd_data_neighbor.uh[2] = _cfd_data.uh[2];
		    _cfd_data_neighbor.uh[3] = _cfd_data.uh[3];
		    _cfd_data_neighbor.uh[4] = p_inf/(_gamma - 1) + 0.5*_cfd_data.r*_cfd_data.vel_size;
		}
	}
}

void EulerBC::symmetric()
{
	wall();
}
