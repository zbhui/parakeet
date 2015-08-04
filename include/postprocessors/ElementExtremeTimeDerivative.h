
#pragma once

#include "ElementExtremeValue.h"

class ElementExtremeTimeDerivative : public ElementExtremeValue
{
public:
	ElementExtremeTimeDerivative(const InputParameters &parameters);

protected:
  virtual void computeQpValue();
};

template<>
InputParameters validParams<ElementExtremeTimeDerivative>();
