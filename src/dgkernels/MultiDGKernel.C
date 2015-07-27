#include "MultiDGKernel.h"

template<>
InputParameters validParams<MultiDGKernel>()
{
	  InputParameters params = validParams<DGKernel>();

	  params.addRequiredParam<std::vector<NonlinearVariableName> >("variables", "多个求解变量");

//	  params.addParam<bool>("use_displaced_mesh", true, "Whether or not this object should use the displaced mesh for computation.  Note that in the case this is true but no displacements are provided in the Mesh block the undisplaced mesh will still be used.");
//	  params.addParamNamesToGroup("use_displaced_mesh", "Advanced");
	  return params;
}

MultiDGKernel::MultiDGKernel(const InputParameters & parameters):
		DGKernel(parameters),
		_variables(getParam<std::vector<NonlinearVariableName> >("variables"))
{

	MooseVariable &val0 = _sys.getVariable(_tid, _variables[0]);
	_n_equation = _variables.size();
	for (size_t i = 0; i < _n_equation; ++i)
	{
		MooseVariable &val = _sys.getVariable(_tid, _variables[i]);
		if(val.feType() != val0.feType()) mooseError("MultiDGKernel中变量的类型不一致");

		_uh.push_back(&val.sln());
		_uh_neighbor.push_back(&val.slnNeighbor());
		_grad_uh.push_back(&val.gradSln());
		_grad_uh_neighbor.push_back(&val.gradSlnNeighbor());
	}
}

void MultiDGKernel::computeResidual()
{
	precalculateResidual();

	DGKernel::computeElemNeighResidual(Moose::Element);
	DGKernel::computeElemNeighResidual(Moose::Neighbor);
}

void MultiDGKernel::computeElemNeighResidual(Moose::DGResidualType type)
{
	  bool is_elem;
	  if (type == Moose::Element)
	    is_elem = true;
	  else
	    is_elem = false;

	  const VariableTestValue & test_space = is_elem ? _test : _test_neighbor;

	  for (_eq = 0; _eq < _n_equation; ++_eq)
	  {
		  int var_number = _sys.getVariable(_tid, _variables[_eq]).number();
		  DenseVector<Number> & re = is_elem ?
				  	  	  	  	  	  	   _assembly.residualBlock(var_number) :
	                                       _assembly.residualBlockNeighbor(var_number);

		  for (_qp=0; _qp<_qrule->n_points(); _qp++)
			  for (_i=0; _i< test_space.size(); _i++)
				  re(_i) += _JxW[_qp]*_coord[_qp]*computeQpResidual(type);

	  }
}

void MultiDGKernel::computeJacobian()
{
	precalculateJacobian();
	DGKernel::computeElemNeighJacobian(Moose::ElementElement);
	DGKernel::computeElemNeighJacobian(Moose::ElementNeighbor);
	DGKernel::computeElemNeighJacobian(Moose::NeighborElement);
	DGKernel::computeElemNeighJacobian(Moose::NeighborNeighbor);
}

void MultiDGKernel::computeElemNeighJacobian(Moose::DGJacobianType type)
{
  const VariableTestValue & test_space = ( type == Moose::ElementElement || type == Moose::ElementNeighbor ) ?
                                         _test : _test_neighbor;
  const VariablePhiValue & loc_phi = ( type == Moose::ElementElement || type == Moose::NeighborElement ) ?
                                       _phi : _phi_neighbor;

  for (_ep = 0; _ep < _n_equation; ++_ep)
  for (_eq = 0; _eq < _n_equation; ++_eq)
  {
	  int var_number_p = _sys.getVariable(_tid, _variables[_ep]).number();
	  int var_number_q = _sys.getVariable(_tid, _variables[_eq]).number();
	  DenseMatrix<Number> & Kxx = type == Moose::ElementElement ?  _assembly.jacobianBlock(var_number_p, var_number_q) :
                              	  type == Moose::ElementNeighbor ? _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, var_number_p, var_number_q) :
                              	  type == Moose::NeighborElement ? _assembly.jacobianBlockNeighbor(Moose::NeighborElement, var_number_p, var_number_q) :
                              	  	  	  	  	  			   	   _assembly.jacobianBlockNeighbor(Moose::NeighborNeighbor, var_number_p, var_number_q);

	  for (_qp=0; _qp<_qrule->n_points(); _qp++)
		  for (_i=0; _i<test_space.size(); _i++)
			  for (_j=0; _j<loc_phi.size(); _j++)
				  Kxx(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian(type);

  	  }
}

void MultiDGKernel::computeOffDiagJacobian(unsigned int jvar)
{
//	mooseError("MultiDGKernel::computeOffDiagJacobian暂不支持");
	precalculateJacobian();

	DGKernel::computeOffDiagElemNeighJacobian(Moose::ElementElement,jvar);
	DGKernel::computeOffDiagElemNeighJacobian(Moose::ElementNeighbor,jvar);
	DGKernel::computeOffDiagElemNeighJacobian(Moose::NeighborElement,jvar);
	DGKernel::computeOffDiagElemNeighJacobian(Moose::NeighborNeighbor,jvar);
}

void MultiDGKernel::computeOffDiagElemNeighJacobian(Moose::DGJacobianType type,unsigned int jvar)
{
	const VariableTestValue & test_space = ( type == Moose::ElementElement || type == Moose::ElementNeighbor ) ?
                                           _test : _test_neighbor;
	const VariablePhiValue & loc_phi = ( type == Moose::ElementElement || type == Moose::NeighborElement ) ?
                                       	   _phi : _phi_neighbor;


	for (_eq = 0; _eq < _n_equation; ++_eq)
	{
		unsigned int var_number = _sys.getVariable(_tid, _variables[_eq]).number();
		DenseMatrix<Number> & Kxx =
								  type == Moose::ElementElement ?  _assembly.jacobianBlock(var_number, jvar) :
	                              type == Moose::ElementNeighbor ? _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, var_number, jvar) :
	                              type == Moose::NeighborElement ? _assembly.jacobianBlockNeighbor(Moose::NeighborElement, var_number, jvar) :
	                            		  	  	  	  	  	  	   _assembly.jacobianBlockNeighbor(Moose::NeighborNeighbor, var_number, jvar);

		for (_qp=0; _qp<_qrule->n_points(); _qp++)
			for (_i=0; _i<test_space.size(); _i++)
				for (_j=0; _j<loc_phi.size(); _j++)
				{
					if (jvar == var_number)
						Kxx(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian(type);
					else
						Kxx(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpOffDiagJacobian(type, jvar);
				}
	}
}

void MultiDGKernel::valueAtLeftFace(Real* ul)
{
	for (size_t eq = 0; eq < _uh.size(); ++eq)
	{
		ul[eq] = (*_uh[eq])[_qp];
	}
}

void MultiDGKernel::valueAtRightFace(Real* ur)
{
	for (size_t eq = 0; eq < _uh_neighbor.size(); ++eq)
	{
		ur[eq] = (*_uh_neighbor[eq])[_qp];
	}
}

void MultiDGKernel::valueGradAtLeftFace(RealGradient* dul)
{
	for (size_t eq = 0; eq < _grad_uh.size(); ++eq)
	{
		dul[eq] = (*_grad_uh[eq])[_qp];
	}
}

void MultiDGKernel::valueGradAtRightFace(RealGradient* dur)
{
	for (size_t eq = 0; eq < _grad_uh_neighbor.size(); ++eq)
	{
		dur[eq] = (*_grad_uh_neighbor[eq])[_qp];
	}
}

Real MultiDGKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	return 0.;
}
