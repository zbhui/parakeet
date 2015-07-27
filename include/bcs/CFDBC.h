
#pragma once

#include "MooseEnum.h"

#include "MultiIntegratedBC.h"
#include "CFDBase.h"

class CFDBC;

template<>
InputParameters validParams<CFDBC>();

/**
 *  CFDBC abstract class
 */
class CFDBC :
public MultiIntegratedBC
{
public:
	CFDBC(const std::string & name, InputParameters params);
	virtual ~CFDBC(){}


	virtual void computeResidual();
	virtual void computeJacobian();
	virtual void computeJacobianBlock(unsigned int jvar);
protected:

//	 virtual void precalculateResidual(){};
//	 virtual void precalculateJacobian() {};

	 virtual void valueAtRightFace(Real *ur);

	 virtual void wallBC(Real *ur) = 0;
	 virtual void farFieldBC(Real *ur) = 0;
	 virtual void symmetricBC(Real *ur) = 0;

private:
	MooseEnum _bc_type;

};
