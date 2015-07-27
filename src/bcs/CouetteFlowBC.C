
#include "CouetteFlowBC.h"

template<>
InputParameters validParams<CouetteFlowBC>()
{
	InputParameters params = validParams<NavierStokesBC>();

	return params;
}

CouetteFlowBC::CouetteFlowBC(const std::string &name, InputParameters parameters):
		NavierStokesBC(name, parameters)
{
}

void CouetteFlowBC::valueAtRightFace(Real *ur)
{
	Real ui[5];
	valueAtLeftFace(ui);
    Point p = _q_point[_qp];
    ur[0] = density(p);
    ur[1] = x_momentum(p);
    ur[2] = y_momentum(p);
    ur[3] = z_momentum(p);
    ur[4] = total_energy(p);
}

Real CouetteFlowBC::density(const Point &p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);

	Real tem = temperature(p);
	Real pre = 1./(_gamma * _mach * _mach);
	return pre * _gamma * _mach * _mach/tem;

}

Real CouetteFlowBC::x_momentum(const Point &p)
{
	Real rho = density(p);
	Real u = p(1)/2.0;
	return rho * u;
}

Real CouetteFlowBC::y_momentum(const Point &p)
{
	Real rho = density(p);
	Real v = 0.;
	return rho * v;
}

Real CouetteFlowBC::z_momentum(const Point &p)
{
	return 0.0;
}

Real CouetteFlowBC::total_energy(const Point &p)
{
	Real rho = density(p);
	Real tem = temperature(p);
	Real pre = 1./(_gamma * _mach * _mach);
	Real u = p(1)/2.0;
	Real v = 0;
	Real w = 0;
	return pre/(_gamma-1)+0.5*rho * (u*u + v*v + w*w);
}

Real CouetteFlowBC::temperature(const Point& p)
{
	Real x = p(0);
	Real y = p(1);
	Real z = p(2);
	Real T2=0.85,T1=0.8;
	return T1 + ( T2 - T1 ) * y / 2 + 0.5 * _prandtl * (_gamma - 1) * _mach * _mach * y / 2 * ( 1 - y / 2 );
}
