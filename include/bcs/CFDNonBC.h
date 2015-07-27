
#pragma once

#include "EulerBC.h"

class CFDNonBC;

template<>
InputParameters validParams<CFDNonBC>();

class CFDNonBC :
public EulerBC
{
public:
	  CFDNonBC(const std::string & name, InputParameters params);
	  virtual ~CFDNonBC(){}

protected:
	virtual void valueAtRightFace(Real *ur);
};
