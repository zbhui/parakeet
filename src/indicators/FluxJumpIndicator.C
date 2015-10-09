
#include "FluxJumpIndicator.h"
#include "CFDProblem.h"

template<>
InputParameters validParams<FluxJumpIndicator>()
{
  InputParameters params = validParams<InternalSideIndicator>();
  params += validParams<TransientInterface>();
  params.addCoupledVar("variables", "多个求解变量");
  params.addParam<Real>("scale", 1, "人工粘性放大尺度");
  return params;
}


FluxJumpIndicator::FluxJumpIndicator(const InputParameters &parameters) :
	InternalSideIndicator(parameters),
    TransientInterface(parameters, "indicators"),
	_cfd_problem(static_cast<CFDProblem&>(_fe_problem)),
	_cfd_data(_cfd_problem),
	_cfd_data_neighbor(_cfd_problem),
	_nl(_fe_problem.getNonlinearSystem()),
	_tid(parameters.get<THREAD_ID>("_tid")),
//	_variables(getParam<std::vector<NonlinearVariableName> >("variables")),
//	_variables(_nl.getVariableNames()),
//	_var_order(_fe_problem.getVariable(_tid, _variables[0]).order()),
    _current_elem_volume(_assembly.elemVolume()),
    _neighbor_elem_volume(_assembly.neighborVolume()),
	_current_side_volume(_assembly.sideElemVolume())
{
//	_n_variables = _variables.size();
	std::cout << _n_variables;
	for (int ivar = 0; ivar < 5; ++ivar)
	{
//		MooseVariable &val = _fe_problem.getVariable(_tid, _variables[ivar]);
		_uh.push_back(&coupledValue("variables", ivar));
		_uh_neighbor.push_back(&coupledNeighborValue("variables", ivar));

//		_uh.push_back(_is_implicit ? &val.sln() : &val.slnOld());
//		_uh_neighbor.push_back(_is_implicit ? &val.slnNeighbor() : &val.slnOldNeighbor());
//		_grad_uh.push_back(_is_implicit ? &val.gradSln(): &val.gradSlnOld());
//		_grad_uh_neighbor.push_back(_is_implicit ? &val.gradSlnNeighbor(): &val.gradSlnOldNeighbor());
	}
}

void FluxJumpIndicator::computeIndicator()
{
	Real sum = 0;
	for (_qp=0; _qp<_qrule->n_points(); _qp++)
		sum += _JxW[_qp]*_coord[_qp]*computeQpIntegral();

	sum /= _current_side_volume;
	Real scale = getParam<Real>("scale");
	sum *= scale;

	{
		Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
		_solution.add(_field_var.nodalDofIndex(), fabs(sum));
		_solution.add(_field_var.nodalDofIndexNeighbor(), fabs(sum));
	}
}

Real FluxJumpIndicator::computeQpIntegral()
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
	for (int eq = 0; eq < 5; ++eq)
	{
		Real flux_jump= (_cfd_data.invis_flux[eq] - _cfd_data_neighbor.invis_flux[eq]) * _normals[_qp];
		indicator += (flux_jump)*(dp_du[eq]);
	}

	Real p = (_cfd_data.p + _cfd_data_neighbor.p)/_cfd_problem._gamma;

	return (indicator)/p;
}

void FluxJumpIndicator::finalize()
{

}

