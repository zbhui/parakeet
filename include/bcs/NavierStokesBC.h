
#pragma once

#include "EulerBC.h"

class NavierStokesBC;

template<>
InputParameters validParams<NavierStokesBC>();

/**
 * NavierStokes BC abstract class
 */
class NavierStokesBC :
public EulerBC
{
public:
	NavierStokesBC(const InputParameters & parameters);
	virtual ~NavierStokesBC(){}


protected:
	  RealVectorValue _penalty[40][10];
	  RealVectorValue _penalty_neighbor[40][10];

	  virtual Real computeQpResidual();
	  virtual Real computeQpJacobian();
	  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

	  virtual void precalculateResidual();
	  virtual void computeJacobianVariable(Real *ul);

	  virtual void wallBC(Real *ur);
	  virtual void farFieldBC(Real *ur);
	  virtual void symmetricBC(Real *ur);
	  virtual void pressureOutBC(Real *ur);
};
