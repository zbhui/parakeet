#include "MultiAuxKernel.h"

template<>
InputParameters validParams<MultiAuxKernel>()
{
	InputParameters params = validParams<AuxKernel>();
	params.addRequiredParam<std::vector<AuxVariableName> >("aux_variables", "多个求解变量");
	return params;
}

MultiAuxKernel::MultiAuxKernel(const InputParameters & parameters):
	AuxKernel(parameters),
	_aux_variables(parameters.get<std::vector<AuxVariableName> >("aux_variables")),
	_ivar(0)
{
}

void MultiAuxKernel::compute()
{
	for(_ivar = 0; _ivar < _aux_variables.size(); ++_ivar)
	{
		MooseVariable &var = _aux_sys.getVariable(_tid, _aux_variables[_ivar]);
		if (isNodal())
		{
			if (_var.isNodalDefined())
			{
				_qp = 0;
				Real value = computeValue();
				var.setNodalValue(value);
			}
		}
		else                     /* elemental variables */
		{
			_n_local_dofs = _var.numberOfDofs();

			if (_n_local_dofs==1)  /* p0 */
			{
				Real value = 0;
				for (_qp=0; _qp<_qrule->n_points(); _qp++)
					value += _JxW[_qp]*_coord[_qp]*computeValue();
				value /= (_bnd ? _current_side_volume : _current_elem_volume);
				var.setNodalValue(value);
			}
			else                   /* high-order */
			{
				_local_re.resize(_n_local_dofs);
				_local_re.zero();
				_local_ke.resize(_n_local_dofs, _n_local_dofs);
				_local_ke.zero();

				for (unsigned int i = 0; i < _test.size(); i++)
					for (_qp = 0; _qp < _qrule->n_points(); _qp++)
					{
						Real t = _JxW[_qp] * _coord[_qp] * _test[i][_qp];
						_local_re(i) += t * computeValue();
						for (unsigned int j = 0; j < _test.size(); j++)
							_local_ke(i, j) += t * _test[j][_qp];
					}

				_local_sol.resize(_n_local_dofs);
				_local_ke.cholesky_solve(_local_re, _local_sol);

				var.setNodalValue(_local_sol);
			}
		}

	}
}

