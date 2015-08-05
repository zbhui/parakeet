
#include "AdiabaticWall.h"

template<>
InputParameters validParams<AdiabaticWall>()
{
	InputParameters params = validParams<CFDBC>();

	return params;
}

AdiabaticWall::AdiabaticWall(const InputParameters & parameters):
		CFDBC(parameters)
{
}

void AdiabaticWall::boundaryCondition()
{
    Real twall = 1.;
    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
    _cfd_data_neighbor.uh[1] = 0;
    _cfd_data_neighbor.uh[2] = 0;
    _cfd_data_neighbor.uh[3] = 0;
    _cfd_data_neighbor.uh[4] = _cfd_data.uh[4];
    _cfd_data_neighbor._prandtl = 1.0E+10;
}
