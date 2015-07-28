//
//#include "ConservationLawFaceKernel.h"
//
//template<>
//InputParameters validParams<ConservationLawFaceKernel>()
//{
//	  InputParameters params = validParams<MultiDGKernel>();
//
//	  return params;
//}
//ConservationLawFaceKernel::ConservationLawFaceKernel(const InputParameters & parameters):
//		MultiDGKernel(parameters)
//{
//}
//
//Real ConservationLawFaceKernel::computeQpResidual(Moose::DGResidualType type)
//{
//	if(type == Moose::Element)
//	{
//		return _flux[_qp][_eq] * _test[_i][_qp];
//	}
//	if(type == Moose::Neighbor)
//	{
//		return -_flux[_qp][_eq] * _test_neighbor[_i][_qp];
//	}
//	mooseError("face flux error.");
//	return 0.;
//}
//
