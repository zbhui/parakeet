
#include "CFDBC.h"

/// CFDBC abstract class
template<>
InputParameters validParams<CFDBC>()
{
	MooseEnum bc_types("wall, far-field, symmetric, pressure-out, isovortex, exact", "exact");  // 边界条件的类型，可以增加

	InputParameters params = validParams<MultiIntegratedBC>();
	params.addRequiredParam<MooseEnum>("bc_type", bc_types, "边界条件");
	return params;
}

CFDBC::CFDBC(const std::string & name, InputParameters parameters):
		MultiIntegratedBC(name, parameters),
		_bc_type(getParam<MooseEnum>("bc_type"))
{

}

void CFDBC::computeResidual()
{
	precalculateResidual();
	IntegratedBC::computeResidual();
}

void CFDBC::computeJacobian()
{
	precalculateJacobian();
	IntegratedBC::computeJacobian();
}

void CFDBC::computeJacobianBlock(unsigned int jvar)
{
	precalculateJacobian();
	IntegratedBC::computeJacobianBlock(jvar);
}

void CFDBC::valueAtRightFace(Real* ur)
{
	if(_bc_type == "wall")
	{
		wallBC(ur);
		return;
	}
	if(_bc_type == "far-field")
	{
		farFieldBC(ur);
		return;
	}
	if(_bc_type == "symmetric")
	{
		symmetricBC(ur);
		return;
	}
	mooseError("未定义的边界条件类型");
}

