#include "MultiDGKernel.h"

template<>
InputParameters validParams<MultiDGKernel>()
{
	  InputParameters params = validParams<DGKernel>();
	  params += validParams<MultiVariableInterface>();
	  return params;
}

MultiDGKernel::MultiDGKernel(const InputParameters & parameters):
		DGKernel(parameters),
		MultiVariableInterface(parameters),
	    _fe_problem(*parameters.get<FEProblem *>("_fe_problem"))
{
}

void MultiDGKernel::computeResidual()
{
	for (_qp=0; _qp<_qrule->n_points(); _qp++)
	{
		precalculateResidual();
		computeElemNeighResidual(Moose::Element);
		computeElemNeighResidual(Moose::Neighbor);
	}
}


void MultiDGKernel::computeElemNeighResidual(Moose::DGResidualType type)
{
	bool is_elem;
	if (type == Moose::Element)
		is_elem = true;
	else
		is_elem = false;

	const VariableTestValue & test_space = is_elem ? _test : _test_neighbor;

	for (unsigned int p = 0; p < _n_equation; ++p)
	{
		DenseVector<Number> & re = is_elem ?
				_assembly.residualBlock(p) :
				_assembly.residualBlockNeighbor(p);

		for (_i=0; _i< test_space.size(); _i++)
			re(_i) += _JxW[_qp]*_coord[_qp]*computeQpResidual(type, p);

	}
}

void MultiDGKernel::computeJacobian()
{
	for (_qp=0; _qp<_qrule->n_points(); _qp++)
	{
		precalculateJacobian();
		computeElemNeighJacobian(Moose::ElementElement);
		computeElemNeighJacobian(Moose::ElementNeighbor);
		computeElemNeighJacobian(Moose::NeighborElement);
		computeElemNeighJacobian(Moose::NeighborNeighbor);
	}
}

void MultiDGKernel::computeElemNeighJacobian(Moose::DGJacobianType type)
{
  const VariableTestValue & test_space = ( type == Moose::ElementElement || type == Moose::ElementNeighbor ) ?
                                         _test : _test_neighbor;
  const VariablePhiValue & loc_phi = ( type == Moose::ElementElement || type == Moose::NeighborElement ) ?
                                       _phi : _phi_neighbor;

  for (unsigned int p = 0; p < _n_equation; ++p)
  for (unsigned int q = 0; q < _n_equation; ++q)
  {
	  DenseMatrix<Number> & Kxx = type == Moose::ElementElement ?  _assembly.jacobianBlock(p, q) :
                              	  type == Moose::ElementNeighbor ? _assembly.jacobianBlockNeighbor(Moose::ElementNeighbor, p, q) :
                              	  type == Moose::NeighborElement ? _assembly.jacobianBlockNeighbor(Moose::NeighborElement, p, q) :
                              	  	  	  	  	  			   	   _assembly.jacobianBlockNeighbor(Moose::NeighborNeighbor, p, q);

		  for (_i=0; _i<test_space.size(); _i++)
			  for (_j=0; _j<loc_phi.size(); _j++)
				  Kxx(_i,_j) += _JxW[_qp]*_coord[_qp]*computeQpJacobian(type, p, q);

  	  }
}

void MultiDGKernel::computeOffDiagJacobian(unsigned int jvar)
{
	if (jvar == _var.number())
	    computeJacobian();
	else
	{
		return;
		mooseError("MultiDGKernel::computeOffDiagJacobian");
	}
}

Real MultiDGKernel::computeQpJacobian(Moose::DGJacobianType type)
{
	mooseError("MultiDGKernel::computeQpJacobian");
	return 0.;
}

Real MultiDGKernel::computeQpResidual(Moose::DGResidualType type)
{
	mooseError("MultiDGKernel::computeQpResidual");
	return 0;
}
