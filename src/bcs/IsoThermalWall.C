
#include "IsoThermalWall.h"

template<>
InputParameters validParams<IsoThermalWall>()
{
	InputParameters params = validParams<CFDBC>();

	return params;
}

IsoThermalWall::IsoThermalWall(const InputParameters & parameters):
		CFDBC(parameters)
{
}

void IsoThermalWall::boundaryCondition()
{

}
