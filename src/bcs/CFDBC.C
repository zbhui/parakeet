
#include "CFDBC.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<CFDBC>()
{
	InputParameters params = validParams<MultiIntegratedBC>();
	params.addParam<Real>("perturbation", 1E-08, "有限差分求Jacobian矩阵的变量增量");
	params.addParam<Real>("penalty", 10, "IP 方法的惩罚值");
	return params;
}

CFDBC::CFDBC(const InputParameters & parameters):
		MultiIntegratedBC(parameters),
		_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
		_cfd_data(_cfd_problem),
		_cfd_data_neighbor(_cfd_problem),
		_lift_data(_cfd_problem),
		_flux_type(_cfd_problem._flux_type),
		_perturbation(getParam<Real>("perturbation")),
		_penalty(getParam<Real>("penalty")),
		_gamma(_cfd_problem._gamma),
		_mach(_cfd_problem._mach),
		_attitude(_cfd_problem._attitude)
{
}

void CFDBC::reinit()
{
	_cfd_data.reinit();
	boundaryCondition();
	_cfd_data_neighbor.reinit();
	fluxRiemann();
	liftOperator();

}

void CFDBC::precalculateResidual()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data.duh[i] = (*_grad_uh[i])[_qp];
		_cfd_data_neighbor.duh[i] = (*_grad_uh[i])[_qp];
	}

	reinit();
}

Real CFDBC::computeQpResidual(unsigned int p)
{
	return _flux[p] * _test[_i][_qp] + _lift[p]*_grad_test[_i][_qp];
}

Real CFDBC::computeQpJacobian(unsigned int p, unsigned int q)
{
	Real r(0);
	r = _flux_jacobi_variable[p][q]*_phi[_j][_qp]*_test[_i][_qp];
	r += _flux_jacobi_grad_variable[p][q]*_grad_phi[_j][_qp]*_test[_i][_qp];
	r += _lift_jacobi_variable[p][q]*_grad_test[_i][_qp]*_phi[_j][_qp];

	return r;
}

void CFDBC::fluxRiemann()
{
	if(_flux_type == 0)
		fluxLaxF();
	else if(_flux_type == 2)
		fluxHLLCPV();
	else
	{
		mooseError("Riemann 通量未定义，请选择Lax-F 或者 HLLC-PV");
	}
}

void CFDBC::fluxLaxF()
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

void CFDBC::fluxHLLCPV()
{
	Real gamma = 1.4;

	Real ql = _cfd_data.vel*_normals[_qp];
	Real qr = _cfd_data_neighbor.vel*_normals[_qp];
	Real pl = _cfd_data.p;
	Real pr = _cfd_data_neighbor.p;
	Real cl = _cfd_data.c;
	Real cr = _cfd_data_neighbor.c;
	Real rl = _cfd_data.r;
	Real rr = _cfd_data_neighbor.r;

	Real rho_bar = (rl + rr)/2.;
	Real c_bar =   (cl + cr)/2.;
	Real p_bar =   (pl + pr)/2.;
	Real q_bar =   (ql + qr)/2.;

	Real p_star = p_bar - 0.5*(qr-ql)*rho_bar*c_bar;
	Real gl;
	if(p_star <= pl)
		gl = 1;
	else
		gl = sqrt(1+(gamma+1)/2./gamma*fabs(p_star/pl -1 ));

	Real gr;
	if(p_star <= pr)
		gr = 1;
	else
		gr = sqrt(1+(gamma+1)/2./gamma*fabs(p_star/pr -1 ));

	Real sl = ql-gl*cl;
	Real sr = qr+gr*cr;
	Real sm = rr*qr*(sr-qr)-rl*ql*(sl-ql)+pl-pr; sm /= rr*(sr-qr)-rl*(sl-ql);
	std::vector<Real> ul_star, ur_star;
	ul_star.push_back(rl*(sl-ql)/(sl-sm));
	ur_star.push_back(rr*(sr-qr)/(sr-sm));

	ul_star.push_back((_cfd_data.uh[1]*         (sl-ql)+(p_star-pl)*_normals[_qp](0))/(sl-sm));
	ur_star.push_back((_cfd_data_neighbor.uh[1]*(sr-qr)+(p_star-pr)*_normals[_qp](0))/(sr-sm));

	ul_star.push_back((_cfd_data.uh[2]*         (sl-ql)+(p_star-pl)*_normals[_qp](1))/(sl-sm));
	ur_star.push_back((_cfd_data_neighbor.uh[2]*(sr-qr)+(p_star-pr)*_normals[_qp](1))/(sr-sm));

	ul_star.push_back((_cfd_data.uh[3]*         (sl-ql)+(p_star-pl)*_normals[_qp](2))/(sl-sm));
	ur_star.push_back((_cfd_data_neighbor.uh[3]*(sr-qr)+(p_star-pr)*_normals[_qp](2))/(sr-sm));


	ul_star.push_back((_cfd_data.uh[4]*         (sl-ql)-pl*ql+p_star*sm)/(sl-sm));
	ur_star.push_back((_cfd_data_neighbor.uh[4]*(sr-qr)-pr*qr+p_star*sm)/(sr-sm));

//	std::cout << sl <<std::endl;
//	std::cout << gr <<std::endl;
//	std::cout << sr <<std::endl;

	if(sl > 0)
	{
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux[p] = _cfd_data.invis_flux[p]*_normals[_qp];
		}
	}
	else if(sl<=0 && 0<sm)
	{
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux[p] = _cfd_data.invis_flux[p]*_normals[_qp] + sl*(ul_star[p]-_cfd_data.uh[p]);
		}
	}

	else if(sm<=0 && 0<sr)
	{
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux[p] = _cfd_data_neighbor.invis_flux[p]*_normals[_qp] + sr*(ur_star[p]-_cfd_data_neighbor.uh[p]);
		}
	}

	else if(sr < 0)
	{
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux[p] = _cfd_data_neighbor.invis_flux[p]*_normals[_qp];
		}
	}

	else
	{
		mooseWarning("HLLC flux error. instead of LF flux.");
		fluxRiemann();
	}

}


void CFDBC::liftOperator()
{
	Real penalty = _penalty*(_var_order*_var_order+1)*_current_side_volume / (_current_elem_volume+_current_elem_volume);

	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_lift_data.uh[i] = (_cfd_data.uh[i] + _cfd_data_neighbor.uh[i])/2.;
		_lift_data.duh[i] = penalty*(_cfd_data.uh[i] - _cfd_data_neighbor.uh[i])*_normals[_qp];
	}

	_lift_data.reinit();

	for (int p = 0; p < _n_equation; ++p)
			_flux[p] -= 0.5*((_cfd_data.vis_flux[p] + _cfd_data_neighbor.vis_flux[p]) - _lift_data.vis_flux[p])*_normals[_qp];
}

void CFDBC::precalculateJacobian()
{
	precalculateResidual();
	for (int q = 0; q < _n_equation; ++q)
	{
		_flux_old[q] = _flux[q];
		_lift_old[q] = _lift[q];
	}

	for (int q = 0; q < _n_equation; ++q)
	{
		_cfd_data.uh[q] += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux_jacobi_variable[p][q] = (_flux[p] - _flux_old[p])/_perturbation;
			_lift_jacobi_variable[p][q] = (_lift[p] - _lift_old[p])/_perturbation;
		}
		_cfd_data.uh[q] -= _perturbation;
	}

	for (int q = 0; q < _n_equation; ++q)
	for (int beta = 0; beta < 3; ++beta)
	{
		_cfd_data.duh[q](beta) += _perturbation;
		reinit();
		for (int p = 0; p < _n_equation; ++p)
		{
			_flux_jacobi_grad_variable[p][q](beta) = (_flux[p] - _flux_old[p])/_perturbation;
		}
		_cfd_data.duh[q](beta) -= _perturbation;
	}
}

void CFDBC::boundaryCondition()
{
	for(int p = 0; p < _n_equation; ++p)
		_cfd_data_neighbor.uh[p] = _cfd_problem.boundaryCondition(_t, _q_point[_qp], p);
}
