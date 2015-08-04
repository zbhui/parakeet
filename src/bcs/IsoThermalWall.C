
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
    Real twall = 1.;
    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
    _cfd_data_neighbor.uh[1] = 0;
    _cfd_data_neighbor.uh[2] = 0;
    _cfd_data_neighbor.uh[3] = 0;
    _cfd_data_neighbor.uh[4] = _cfd_data.uh[0]*twall/_gamma/(_gamma-1)/_mach/_mach;
}
