
#include "FarFieldPressure.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<FarFieldPressure>()
{
	InputParameters params = validParams<CFDBC>();
	return params;
}

FarFieldPressure::FarFieldPressure(const InputParameters & parameters):
	CFDBC(parameters),
	_gamma(_cfd_problem._gamma),
	_mach(_cfd_problem._mach),
	_rho_inf(1.0),
	_vel_inf(1.0)
{
}

void FarFieldPressure::boundaryCondition()
{
	Eigen::Vector3d vel_inf(1, 0, 0);
	Real p_inf = 1/_gamma/_mach/_mach;

	CFDDataPack cfd_data_infity(_cfd_problem);
	cfd_data_infity.r = _rho_inf;
	cfd_data_infity.mom(0) = 1.;
	cfd_data_infity.mom(1) = 0.;
	cfd_data_infity.mom(2) = 0.;
	cfd_data_infity.re = 1/_gamma/_mach/_mach + 0.5*cfd_data_infity.mom.size()/1.;
	cfd_data_infity.vel = cfd_data_infity.mom/_rho_inf;
//	cfd_data_infity.reinit();

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
