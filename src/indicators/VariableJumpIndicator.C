
#include "VariableJumpIndicator.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<VariableJumpIndicator>()
{
  InputParameters params = validParams<FluxJumpIndicator>();
  return params;
}

VariableJumpIndicator::VariableJumpIndicator(const InputParameters &parameters) :
	FluxJumpIndicator(parameters)
{
}

void VariableJumpIndicator::computeIndicator()
{
	Real sum = 0;
	for (_qp=0; _qp<_qrule->n_points(); _qp++)
		sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();

	sum /= _current_side_volume;

	Real scale = getParam<Real>("scale");

	{
		Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

		Real hmax, hmax_neighbor, ind, ind_neighbor;
		hmax = _current_elem->hmax();
		hmax_neighbor = _neighbor_elem->hmax();

		ind = scale*fabs(sum)*hmax;
		ind_neighbor = scale*fabs(sum)*hmax_neighbor;

		_solution.add(_field_var.nodalDofIndex(), ind);
		_solution.add(_field_var.nodalDofIndexNeighbor(), ind_neighbor);

	}
}

Real VariableJumpIndicator::computeQpIntegral()
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		_cfd_data.uh[i] = (*_uh[i])[_qp];
		_cfd_data_neighbor.uh[i] = (*_uh_neighbor[i])[_qp];
	}
	_cfd_data.reinitInviscous();
	_cfd_data_neighbor.reinitInviscous();

	Real dp_du[] = {
			(0.5*_cfd_data.vel.size_sq() + 0.5*_cfd_data_neighbor.vel.size_sq()),
			-(_cfd_data.vel(0)+_cfd_data_neighbor.vel(0)),
			-(_cfd_data.vel(1)+_cfd_data_neighbor.vel(1)),
			-(_cfd_data.vel(2)+_cfd_data_neighbor.vel(2)),
			(1+1)
	};

	Real indicator(0);
	for (int eq = 0; eq < 1	; ++eq)
	{
		Real flux_jump= (_cfd_data.uh[eq] - _cfd_data_neighbor.uh[eq])/(_cfd_data.uh[eq] + _cfd_data_neighbor.uh[eq]);
//		flux_jump *= (_cfd_data.vel+_cfd_data_neighbor.vel)*_normals[_qp] + _cfd_data.c+_cfd_data_neighbor.c;
		indicator += (flux_jump);//*(dp_du[eq]);
	}

	Real p = 1;//(_cfd_data.p + _cfd_data_neighbor.p)/_cfd_problem._gamma;

	return (indicator)/p;
}

void VariableJumpIndicator::finalize()
{

}

