
#include "FarFieldPressure.h"
#include "CFDProblem.h"
using namespace Eigen;

template<>
InputParameters validParams<FarFieldPressure>()
{
	InputParameters params = validParams<CFDBC>();
	return params;
}

FarFieldPressure::FarFieldPressure(const InputParameters & parameters):
	CFDBC(parameters),
	_rho_inf(1.0),
	_vel_inf(1.0),
	_tem_inf(1.0)
{
}

void FarFieldPressure::boundaryCondition()
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
		    _cfd_data_neighbor.uh[0] = infity.r;
		    _cfd_data_neighbor.uh[1] = infity.mom(0);
		    _cfd_data_neighbor.uh[2] = infity.mom(1);
		    _cfd_data_neighbor.uh[3] = infity.mom(2);
		    _cfd_data_neighbor.uh[4] = _cfd_data.p/(_gamma-1) + 0.5 * _rho_inf*vel_inf.squaredNorm();
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
