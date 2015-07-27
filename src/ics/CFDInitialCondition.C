
#include "CFDInitialCondition.h"

template<>
InputParameters validParams<CFDInitialCondition>()
{
  InputParameters params = validParams<InitialCondition>();
  params += validParams<CFDBase>();
//  params.addRequiredParam<unsigned int>("component", "方程分量");

  return params;
}

CFDInitialCondition::CFDInitialCondition(const InputParameters & parameters) :
    InitialCondition(parameters),
    CFDBase(parameters)
//    _component(getParam<unsigned int>("component"))
{}

Real CFDInitialCondition::value(const Point & p)
{
	std::string var_name = _var.name();
	if(var_name == "rho")
			return density(p);
	if(var_name == "momentum_x")
			return x_momentum(p);
	if(var_name == "momentum_y")
			return y_momentum(p);
	if(var_name == "momentum_z")
			return z_momentum(p);
	if(var_name == "rhoe")
			return total_energy(p);

	mooseError(var_name << "变量名不存在");
	return 0.;

//switch (_component) {
//	case 0:
//		return density(p);
//		break;
//	case 1:
//		return x_momentum(p);
//		break;
//	case 2:
//		return y_momentum(p);
//		break;
//	case 3:
//		return z_momentum(p);
//	case 4:
//		return total_energy(p);
//		break;
//	default:
//		return 0.0;
//		mooseError("不可用的分量" << _component);
//		break;
//	}
}
