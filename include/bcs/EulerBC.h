
#pragma once
#include "CFDBC.h"

class EulerBC;

template<>
InputParameters validParams<EulerBC>();

/**
 * Euler BC abstract class
 */
class EulerBC :
public CFDBC,
public CFDBase
{
public:
	  EulerBC(const std::string & name, InputParameters params);
	  virtual ~EulerBC(){}


protected:
	  Real _flux[40][10];
	  Real _ul[40][10], _ur[40][10];
	  Real _jacobi_variable[40][10][10];

	  void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	  virtual Real computeQpResidual();
	  virtual Real computeQpJacobian();
	  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

	  virtual void precalculateResidual();
	  virtual void precalculateJacobian();

	  virtual void wallBC(Real *ur);
	  virtual void farFieldBC(Real *ur);
	  virtual void symmetricBC(Real *ur);
};
