#include "FarFieldRiemann.h"
#include "Attitude.h"
using namespace Eigen;

template<>
InputParameters validParams<FarFieldRiemann>()
{
	InputParameters params = validParams<CFDBC>();
	return params;
}

FarFieldRiemann::FarFieldRiemann(const InputParameters & parameters):
		CFDBC(parameters),
		_rho_inf(1.0),
		_vel_inf(1.0),
		_tem_inf(1.0)
{
}

void FarFieldRiemann::boundaryCondition()
{
	Real p_inf = _rho_inf*_tem_inf/_gamma/_mach/_mach;
	Vector3d vel_inf = _vel_inf*(_attitude.earthFromWind()*Vector3d::UnitX());
	if(_mesh.dimension() == 2)
	{
		vel_inf(2) = 0.;
	}
	else if(_mesh.dimension() == 1)
	{
		vel_inf(1) = 0.;
		vel_inf(2) = 0.;
	}


	CFDDataPack infity(_cfd_problem);
	infity.r = _rho_inf;
	infity.mom(0) = _rho_inf*vel_inf(0);
	infity.mom(1) = _rho_inf*vel_inf(1);
	infity.mom(2) = _rho_inf*vel_inf(2);
	infity.re = p_inf/(_gamma-1) + 0.5*_rho_inf*vel_inf.squaredNorm();
	infity.vel = infity.mom/_rho_inf;

	Real vnl = _cfd_data.vel*_normals[_qp];
	Real vnr = infity.vel*_normals[_qp];
	if(vnr <= 0)  // 入口
	{
		if (_cfd_data.m > 1) //超音速
		{
		    _cfd_data_neighbor.uh[0] = infity.r;
		    _cfd_data_neighbor.uh[1] = infity.mom(0);
		    _cfd_data_neighbor.uh[2] = infity.mom(1);
		    _cfd_data_neighbor.uh[3] = infity.mom(2);
		    _cfd_data_neighbor.uh[4] = infity.re;
		}
		else	//亚音速
		{
			Real cl = _cfd_data.c;
			Real cr = sqrt(_gamma*p_inf/_rho_inf);
			Real s = p_inf/pow(_rho_inf, _gamma);
			Real Rp = -vnr+2.0*cr/(_gamma-1);
			Real Rm = -vnl-2.0*cl/(_gamma-1);
			Real vnb = -(Rp + Rm)/2.0;
			Real cb = (Rp - Rm)*(_gamma - 1)/4.0;

			Real rho = pow((cb*cb)/(s*_gamma), 1.0/(_gamma - 1));
			RealVectorValue mom = rho*(Point(vel_inf(0), vel_inf(1), vel_inf(2)) + (vnb-vnr)*_normals[_qp]);
			_cfd_data_neighbor.uh[0] = rho;
			_cfd_data_neighbor.uh[1] = mom(0);
			_cfd_data_neighbor.uh[2] = mom(1);
			_cfd_data_neighbor.uh[3] = mom(2);
			_cfd_data_neighbor.uh[4] = cb*cb*rho/_gamma/(_gamma - 1)+ 0.5*mom.size_sq()/rho;
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
			Real cl = _cfd_data.c;
			Real cr = sqrt(_gamma*p_inf/_rho_inf);
			Real s =_cfd_data.s;
			Real Rp = vnl + 2*cl/(_gamma - 1);
			Real Rm = vnr - 2*cr/(_gamma - 1);
			Real vnb = (Rp + Rm) / 2.0;
			Real cb = (Rp - Rm) * (_gamma - 1) / 4.0;

			vnb = (Rp + Rm) / 2.0;
			cb = (Rp - Rm) * (_gamma - 1) / 4.0;

			Real rho = pow((cb * cb) / (s * _gamma), 1.0 / (_gamma - 1));
			RealVectorValue mom = rho*(_cfd_data.vel+(vnb-vnl)*_normals[_qp]);
			_cfd_data_neighbor.uh[0] = rho;
			_cfd_data_neighbor.uh[1] = mom(0);
			_cfd_data_neighbor.uh[2] = mom(1);
			_cfd_data_neighbor.uh[3] = mom(2);
			_cfd_data_neighbor.uh[4] = cb*cb*rho/_gamma/(_gamma - 1)+ 0.5*mom.size_sq()/rho;
		}
	}
}
