//
//#pragma once
//
//// MOOSE Includes
//#include "CFDInitialCondition.h"
//
//// Forward Declarations
//class CouetteFlowIC;
//
//template<>
//InputParameters validParams<CouetteFlowIC>();
//
///**
// * 等熵涡初始条件
// */
//class CouetteFlowIC : public CFDInitialCondition
//{
//public:
//
//  CouetteFlowIC(const InputParameters & parameters);
//
//private:
////  Point &centre;
//protected:
//  virtual Real density(const Point &p);
//  virtual Real x_momentum(const Point &p);
//  virtual Real y_momentum(const Point &p);
//  virtual Real z_momentum(const Point &p);
//  virtual Real total_energy(const Point &p);
//
//private:
//  Real temperature(const Point &p);
//};
