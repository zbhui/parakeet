
#include "SlipWall.h"

template<>
InputParameters validParams<SlipWall>()
{
	InputParameters params = validParams<CFDBC>();

	return params;
}

SlipWall::SlipWall(const InputParameters & parameters):
		CFDBC(parameters)
{
}

void SlipWall::boundaryCondition()
{
	const Point &normal = _normals[_qp];
    Real  vn = _cfd_data.mom*normal;

    _cfd_data_neighbor.uh[0] = _cfd_data.uh[0];
    _cfd_data_neighbor.uh[1] = _cfd_data.uh[1] - 2*vn*normal(0);
    _cfd_data_neighbor.uh[2] = _cfd_data.uh[2] - 2*vn*normal(1);
    _cfd_data_neighbor.uh[3] = _cfd_data.uh[3] - 2*vn*normal(2);
    _cfd_data_neighbor.uh[4] = _cfd_data.uh[4];
}
