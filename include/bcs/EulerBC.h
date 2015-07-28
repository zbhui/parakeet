
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
	  EulerBC(const InputParameters & params);
	  virtual ~EulerBC(){}


protected:
	  Real _flux[40][10];
	  Real _ul[40][10], _ur[40][10];
	  Real _jacobi_variable[40][10][10];

	  void fluxRiemann(Real *flux, Real *ul, Real *ur, Point &normal);
	  virtual Real computeQpResidual(unsigned int p);
	  virtual Real computeQpJacobian(unsigned int p, unsigned int q);

	  virtual void precalculateResidual();
	  virtual void precalculateJacobian();

	  virtual void wallBC(Real *ur);
	  virtual void farFieldBC(Real *ur);
	  virtual void symmetricBC(Real *ur);
};
