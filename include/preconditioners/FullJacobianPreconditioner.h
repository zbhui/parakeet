
#pragma once

#include "MoosePreconditioner.h"

class FullJacobianPreconditioner;

template<>
InputParameters validParams<FullJacobianPreconditioner>();

/**
 * Single matrix preconditioner.
 */
class FullJacobianPreconditioner : public MoosePreconditioner
{
public:
  FullJacobianPreconditioner(const std::string & name, InputParameters params);
  virtual ~FullJacobianPreconditioner(){};
};

