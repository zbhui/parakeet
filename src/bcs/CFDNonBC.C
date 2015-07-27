#include "CFDNonBC.h"

template<>
InputParameters validParams<CFDNonBC>()
{
	  InputParameters params = validParams<EulerBC>();

	  return params;
}
CFDNonBC::CFDNonBC(const std::string & name, InputParameters parameters):
		EulerBC(name, parameters)
{
}

void CFDNonBC::valueAtRightFace(Real *ur)
{
	for (size_t i = 0; i < _uh.size(); ++i)
	{
		ur[i] = (*_uh[i])[_qp];
	}
}
