#include "geometry.h"
#include <math.h>

/**********************************************************************/

Geometry::Geometry()
{
	metric = M_SCHWARZSCHILD;
	M = 0.1;
	a = 0.0;
}

Geometry::Geometry(int a)
{
	metric = a;
	M = 0.1;
	a = 0.0;
}

Geometry::~Geometry()
{
}

/**********************************************************************/

int Geometry::setMetric(int m)
{
	switch(m)
	{
	case M_SCHWARZSCHILD:
	case M_KERR_EF:
		metric = m;
		return 0;
	default:
		return 1;
	}
}

double Geometry::g_v(vector4 u, vector4 v, vector4 pos)
{
	return toDouble((g(pos)*Vector(u)*Vector(v)).contracted(0,2).contracted(0,1));
}

double Geometry::g_inv_v(vector4 u, vector4 v, vector4 pos)
{
	return toDouble((g_inv(pos)*Covector(u)*Covector(v)).contracted(0,2).contracted(0,1));
}

vector4 Geometry::Gamma_v(vector4 u, vector4 v, vector4 pos)
{
	return toVector4((Gamma(pos)*Vector(u)*Vector(v)).contracted(1,3).contracted(1,2));
}

/**********************************************************************/

Tensor Geometry::g(vector4 x)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return g_schw(x);
	case M_KERR_EF:
		return g_kerr_ef(x);
	default:
		return Tensor();
	}
}

Tensor Geometry::g_inv(vector4 x)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return g_inv_schw(x);
	case M_KERR_EF:
		return g_inv_kerr_ef(x);
	default:
		return Tensor();
	}
}

Tensor Geometry::Gamma(vector4 x)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return Gamma_schw(x);
	case M_KERR_EF:
		return Gamma_kerr_ef(x);
	default:
		return Tensor();
	}
}

vector4 Geometry::convertTo(vector4 x, int c)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return convertTo_schw(x,c);
	case M_KERR_EF:
		return convertTo_kerr_ef(x,c);
	default:
		return vector4();
	}
}

vector4 Geometry::convertFrom(vector4 x, int c)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return convertFrom_schw(x,c);
	case M_KERR_EF:
		return convertFrom_kerr_ef(x,c);
	default:
		return vector4();
	}
}

Tensor Geometry::Jacobian(int c, vector4 x)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return Jacobian_schw(c,x);
	case M_KERR_EF:
		return Jacobian_kerr_ef(c,x);
	default:
		return Tensor();
	}
}

Tensor Geometry::invJacobian(int c, vector4 x)
{
	switch(metric)
	{
	case M_SCHWARZSCHILD:
		return invJacobian_schw(c,x);
	case M_KERR_EF:
		return invJacobian_kerr_ef(c,x);
	default:
		return Tensor();
	}
}

/**********************************************************************/

double Geometry::getM()
{
	return M;
}

void Geometry::setM(double m)
{
	M = m;
}

double Geometry::getA()
{
	return a;
}

void Geometry::setA(double m)
{
	a = m;
}

/**********************************************************************/

Tensor Geometry::g_schw(vector4 pos)
{
	bool a[] = { false, false };
	Tensor res(2, a);
	WRITE_COMPONENT(res, IND(0,0), 1.0-2*M/pos[1]);
	WRITE_COMPONENT(res, IND(1,1), -1.0/(1.0-2*M/pos[1]));
	WRITE_COMPONENT(res, IND(2,2), -pos[1]*pos[1]);
	WRITE_COMPONENT(res, IND(3,3), -pos[1]*pos[1]*sin(pos[2])*sin(pos[2]));
	return res;
}

Tensor Geometry::g_inv_schw(vector4 pos)
{
	bool a[] = {true, true};
	Tensor res(2, a);
	WRITE_COMPONENT(res, IND(0,0), 1.0/(1.0-2*M/pos[1]));
	WRITE_COMPONENT(res, IND(1,1), -1.0+2*M/pos[1]);
	WRITE_COMPONENT(res, IND(2,2), -1.0/(pos[1]*pos[1]));
	WRITE_COMPONENT(res, IND(3,3), -1.0/(pos[1]*pos[1]*sin(pos[2])*sin(pos[2])));
	return res;
}

Tensor Geometry::Gamma_schw(vector4 pos)
{
	bool a[] = {true, false, false};
	Tensor res(3, a);

	WRITE_COMPONENT(res, IND(0,0,1), M/(pos[1]*(pos[1]-2*M)));
	WRITE_COMPONENT(res, IND(0,1,0), M/(pos[1]*(pos[1]-2*M)));

	WRITE_COMPONENT(res, IND(1,0,0), M/(pos[1]*pos[1])*(1-2*M/pos[1]));
	WRITE_COMPONENT(res, IND(1,1,1), -M/(pos[1]*(pos[1]-2*M)));
	WRITE_COMPONENT(res, IND(1,2,2), 2*M-pos[1]);
	WRITE_COMPONENT(res, IND(1,3,3), (2*M-pos[1])*sin(pos[2])*sin(pos[2]));

	WRITE_COMPONENT(res, IND(2,1,2), 1.0/pos[1]);
	WRITE_COMPONENT(res, IND(2,2,1), 1.0/pos[1]);
	WRITE_COMPONENT(res, IND(2,3,3), -sin(pos[2])*cos(pos[2]));

	WRITE_COMPONENT(res, IND(3,1,3), 1.0/pos[1]);
	WRITE_COMPONENT(res, IND(3,3,1), 1.0/pos[1]);
	WRITE_COMPONENT(res, IND(3,2,3), cos(pos[2])/sin(pos[2]));
	WRITE_COMPONENT(res, IND(3,3,2), cos(pos[2])/sin(pos[2]));

	return res;
}

/**********************************************************************/

vector4 Geometry::convertTo_schw(vector4 x, int coords)
{
	vector4 pos;

	switch(coords)
	{
	case COORD_MAIN:
	case COORD_SPHERICAL:
		return x;
	case COORD_CARTESIAN:
		pos[0] = x[0];
		pos[1] = x[1]*sin(x[2])*cos(x[3]);
		pos[2] = x[1]*sin(x[2])*sin(x[3]);
		pos[3] = x[1]*cos(x[2]);
		return pos;
	default:
		return vector4();
	}
}

vector4 Geometry::convertFrom_schw(vector4 x, int coords)
{
	vector4 pos;
	switch(coords)
	{
	case COORD_MAIN:
	case COORD_SPHERICAL:
		return x;
	case COORD_CARTESIAN:
		pos[0] = x[0];
		pos[1] = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
		pos[2] = acos(x[3]/pos[1]);
		pos[3] = atan2(x[2],x[1]);
		return pos;
	default:
		return vector4();
	}
	return pos;
}

Tensor Geometry::Jacobian_schw(int coord, vector4 x)
{
	return Tensor();
}

Tensor Geometry::invJacobian_schw(int coord, vector4 x)
{
	return Tensor();
}

/**********************************************************************/

Tensor Geometry::g_kerr_ef(vector4 pos)
{
	bool a1[] = { false, false };
	double rho2;
	Tensor res(2, a1);

	double r, th, sinth;
	r = pos[1];
	th = pos[2];
	rho2 = r*r + a*a*cos(th)*cos(th);
	sinth = sin(th)*sin(th);

	WRITE_COMPONENT(res, IND(0,0), 1.0-2*M*r/rho2);
	WRITE_COMPONENT(res, IND(0,1), -1.0);
	WRITE_COMPONENT(res, IND(1,0), -1.0);
	WRITE_COMPONENT(res, IND(0,3), 2*M*r*a*sinth/rho2);
	WRITE_COMPONENT(res, IND(3,0), 2*M*r*a*sinth/rho2);
	WRITE_COMPONENT(res, IND(1,3), a*sinth);
	WRITE_COMPONENT(res, IND(3,1), a*sinth);
	WRITE_COMPONENT(res, IND(3,3), -(r*r+a*a+2*M*r*a*a*sinth/rho2)*sinth);
	WRITE_COMPONENT(res, IND(2,2), -rho2);
	return res;
}

Tensor Geometry::g_inv_kerr_ef(vector4 pos)
{
	bool a1[] = {true, true};
	Tensor res(2, a1);

	double r, th, sinth, delta, rho2;
	r = pos[1];
	th = pos[2];
	rho2 = r*r + a*a*cos(th)*cos(th);
	sinth = sin(th)*sin(th);
	delta = r*r + a*a - 2*M*r;

	WRITE_COMPONENT(res, IND(0,0), -a*a/rho2*sinth);
	WRITE_COMPONENT(res, IND(0,1), -(r*r+a*a)/rho2);
	WRITE_COMPONENT(res, IND(1,0), -(r*r+a*a)/rho2);
	WRITE_COMPONENT(res, IND(0,3), -a/rho2);
	WRITE_COMPONENT(res, IND(3,0), -a/rho2);
	WRITE_COMPONENT(res, IND(1,1), -delta/rho2);
	WRITE_COMPONENT(res, IND(1,3), -a/rho2);
	WRITE_COMPONENT(res, IND(3,1), -a/rho2);
	WRITE_COMPONENT(res, IND(3,3), -1.0/rho2/sinth);
	WRITE_COMPONENT(res, IND(2,2), -1.0/rho2);
	return res;
}

Tensor Geometry::Gamma_kerr_ef(vector4 pos)
{
	bool a1[] = {true, false, false};
	Tensor res(3, a1);

	double r = pos[1];
	double th = pos[2];
	double rho2, rho4, rho6, rhobar, sinth, ra, sincos, delta;

	rho2 = r*r + a*a*cos(th)*cos(th);
	rho4 = rho2*rho2;
	rho6 = rho4*rho2;
	ra = r*r + a*a;
	rhobar = r*r - a*a*cos(th)*cos(th);
	sinth = sin(th)*sin(th);
	sincos = sin(th)*cos(th);
	delta = r*r + a*a - 2*M*r;

	WRITE_COMPONENT(res, IND(0,0,0), M*ra*rhobar/rho6);
	WRITE_COMPONENT(res, IND(0,0,2), -2*M*r*a*a/rho4*sincos);
	WRITE_COMPONENT(res, IND(0,2,0), -2*M*r*a*a/rho4*sincos);
	WRITE_COMPONENT(res, IND(0,0,3), -M*a*ra*rhobar*sinth/rho6);
	WRITE_COMPONENT(res, IND(0,3,0), -M*a*ra*rhobar*sinth/rho6);
	WRITE_COMPONENT(res, IND(0,1,2), -a*a/rho2*sincos);
	WRITE_COMPONENT(res, IND(0,2,1), -a*a/rho2*sincos);
	WRITE_COMPONENT(res, IND(0,1,3), a*r/rho2*sinth);
	WRITE_COMPONENT(res, IND(0,3,1), a*r/rho2*sinth);
	WRITE_COMPONENT(res, IND(0,2,2), -ra/r/rho2);
	WRITE_COMPONENT(res, IND(0,2,3), 2*M*r*a*a*a*sinth*sincos/rho4);
	WRITE_COMPONENT(res, IND(0,3,2), 2*M*r*a*a*a*sinth*sincos/rho4);
	WRITE_COMPONENT(res, IND(0,3,3), ra/rho2*(M*rhobar*a*a*sinth*sinth/rho4-r*sinth));

	WRITE_COMPONENT(res, IND(1,0,0), M*rhobar*delta/rho6);
	WRITE_COMPONENT(res, IND(1,0,1), -M*rhobar/rho4);
	WRITE_COMPONENT(res, IND(1,1,0), -M*rhobar/rho4);
	WRITE_COMPONENT(res, IND(1,0,3), -M*rhobar/rho6*delta*a*sinth);
	WRITE_COMPONENT(res, IND(1,3,0), -M*rhobar/rho6*delta*a*sinth);
	WRITE_COMPONENT(res, IND(1,1,2), -a*a/rho2*sincos);
	WRITE_COMPONENT(res, IND(1,2,1), -a*a/rho2*sincos);
	WRITE_COMPONENT(res, IND(1,1,3), (rho2*r+M*rhobar)/rho4*a*sinth);
	WRITE_COMPONENT(res, IND(1,3,1), (rho2*r+M*rhobar)/rho4*a*sinth);
	WRITE_COMPONENT(res, IND(1,2,2), -delta/rho2/r);
	WRITE_COMPONENT(res, IND(1,3,3), delta*sinth/rho6*(M*a*a*rhobar*sinth-r*rho4));

	WRITE_COMPONENT(res, IND(2,0,0), -2*M*r*a*a*sincos/rho6);
	WRITE_COMPONENT(res, IND(2,0,3), 2*M*r*a*ra*sincos/rho6);
	WRITE_COMPONENT(res, IND(2,3,0), 2*M*r*a*ra*sincos/rho6);
	WRITE_COMPONENT(res, IND(2,1,2), r/rho2);
	WRITE_COMPONENT(res, IND(2,2,1), r/rho2);
	WRITE_COMPONENT(res, IND(2,1,3), a*sincos/rho2);
	WRITE_COMPONENT(res, IND(2,3,1), a*sincos/rho2);
	WRITE_COMPONENT(res, IND(2,2,2), -a*a*sincos/rho2);
	WRITE_COMPONENT(res, IND(2,3,3), -sincos/rho6*(rho4*ra+2*M*r*a*a*sinth*(ra+rho2)));

	WRITE_COMPONENT(res, IND(3,0,0), M*a*rhobar/rho6);
	WRITE_COMPONENT(res, IND(3,0,2), -2*M*r*a*cos(th)/sin(th)/rho4);
	WRITE_COMPONENT(res, IND(3,2,0), -2*M*r*a*cos(th)/sin(th)/rho4);
	WRITE_COMPONENT(res, IND(3,0,3), -M*a*a*rhobar*sinth/rho6);
	WRITE_COMPONENT(res, IND(3,3,0), -M*a*a*rhobar*sinth/rho6);
	WRITE_COMPONENT(res, IND(3,1,2), -rho2*a*cos(th)/sin(th));
	WRITE_COMPONENT(res, IND(3,2,1), -rho2*a*cos(th)/sin(th));
	WRITE_COMPONENT(res, IND(3,1,3), r/rho2);
	WRITE_COMPONENT(res, IND(3,3,1), r/rho2);
	WRITE_COMPONENT(res, IND(3,2,2), -a/r/rho2);
	WRITE_COMPONENT(res, IND(3,2,3), cos(th)/sin(th)/rho4*(rho2*r*r+2*M*r*a*a*sinth+rhobar*a*a*cos(th)*cos(th)));
	WRITE_COMPONENT(res, IND(3,3,2), cos(th)/sin(th)/rho4*(rho2*r*r+2*M*r*a*a*sinth+rhobar*a*a*cos(th)*cos(th)));
	WRITE_COMPONENT(res, IND(3,3,3), a*sinth/rho6*(M*a*a*rhobar*sinth-r*rho4));

	return res;
}

/**********************************************************************/

vector4 Geometry::convertTo_kerr_ef(vector4 x, int coord)
{
	vector4 pos,x2;
	double logr, delta;

	switch(coord)
	{
	case COORD_MAIN:
		return x;
	case COORD_SPHERICAL:
		logr = log(fabs((x[1]-M-sqrt(M*M-a*a))/(x[1]-M+sqrt(M*M-a*a))));
		delta = x[1]*x[1] - 2*M*x[1] + a*a;
		pos[0] = x[0] - x[1] - M*log(delta) - M*M/sqrt(M*M-a*a)*logr;
		pos[1] = x[1];
		pos[2] = x[2];
		pos[3] = x[3] - a/2/sqrt(M*M-a*a)*logr;
		return pos;
	case COORD_CARTESIAN:
		x2 = convertTo_kerr_ef(x, COORD_SPHERICAL);
		pos[0] = x2[0];
		pos[1] = x2[1]*sin(x2[2])*cos(x2[3]);
		pos[2] = x2[1]*sin(x2[2])*sin(x2[3]);
		pos[3] = x2[1]*cos(x2[2]);
		return pos;
	default:
		return vector4();
	}
}

vector4 Geometry::convertFrom_kerr_ef(vector4 x, int coord)
{
	vector4 x2, pos;
	double logr, delta;

	switch(coord)
	{
	case COORD_MAIN:
		return x;
	case COORD_SPHERICAL:
		logr = log(fabs((x[1]-M-sqrt(M*M-a*a))/(x[1]-M+sqrt(M*M-a*a))));
		delta = x[1]*x[1] - 2*M*x[1] + a*a;
		pos[0] = x[0] + x[1] + M*log(delta) + M*M/sqrt(M*M-a*a)*logr;
		pos[1] = x[1];
		pos[2] = x[2];
		pos[3] = x[3] + a/2/sqrt(M*M-a*a)*logr;
		return pos;
	case COORD_CARTESIAN:
		x2[0] = x[0];
		x2[1] = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]);
		x2[2] = acos(x[3]/pos[1]);
		x2[3] = atan2(x[2],x[1]);
		return convertFrom_kerr_ef(x2,COORD_SPHERICAL);
	default:
		return vector4();
	}
}

Tensor Geometry::Jacobian_kerr_ef(int coord, vector4 x)
{
	bool a1[] = {true, false};
	Tensor jac(2, a1);
	double logr = log(fabs((x[1]-M-sqrt(M*M-a*a))/(x[1]-M+sqrt(M*M-a*a))));
	double phi2 = x[3] - a/2/sqrt(M*M-a*a)*logr;
	double r = x[1];
	double th = x[2];
	double delta = r*r - 2*M*r + a*a;

	switch(coord)
	{
	case COORD_MAIN:
		WRITE_COMPONENT(jac, IND(0,0), 1.0);
		WRITE_COMPONENT(jac, IND(1,1), 1.0);
		WRITE_COMPONENT(jac, IND(2,2), 1.0);
		WRITE_COMPONENT(jac, IND(3,3), 1.0);
		break;
	case COORD_SPHERICAL:
		WRITE_COMPONENT(jac, IND(0,0), 1.0);
		WRITE_COMPONENT(jac, IND(0,1), (r*r+a*a)/delta);
		WRITE_COMPONENT(jac, IND(1,1), 1.0);
		WRITE_COMPONENT(jac, IND(2,2), 1.0);
		WRITE_COMPONENT(jac, IND(3,1), a/delta);
		WRITE_COMPONENT(jac, IND(3,3), 1.0);
		break;
	case COORD_CARTESIAN:
		WRITE_COMPONENT(jac, IND(0,0), 1.0);
		WRITE_COMPONENT(jac, IND(0,1), (r*r+a*a)/delta*sin(th)*cos(phi2));
		WRITE_COMPONENT(jac, IND(0,2), (r*r+a*a)/delta*sin(th)*sin(phi2));
		WRITE_COMPONENT(jac, IND(0,3), (r*r+a*a)/delta*cos(th));
		WRITE_COMPONENT(jac, IND(1,1), sin(th)*cos(phi2));
		WRITE_COMPONENT(jac, IND(1,2), sin(th)*sin(phi2));
		WRITE_COMPONENT(jac, IND(1,3), cos(th));
		WRITE_COMPONENT(jac, IND(2,1), cos(th)*cos(phi2)/r);
		WRITE_COMPONENT(jac, IND(2,2), cos(th)*sin(phi2)/r);
		WRITE_COMPONENT(jac, IND(2,3), -sin(th)/r);
		WRITE_COMPONENT(jac, IND(3,1), a*cos(phi2)*sin(th)/delta-sin(phi2)/r/sin(th));
		WRITE_COMPONENT(jac, IND(3,2), cos(phi2)/r/sin(th)+a*sin(phi2)*sin(th)/delta);
		WRITE_COMPONENT(jac, IND(3,3), a*cos(th)/delta);
		break;
	}

	return jac;
}

Tensor Geometry::invJacobian_kerr_ef(int coord, vector4 x)
{
	bool a1[] = {true, false};
	Tensor jac(2, a1);
	double logr = log(fabs((x[1]-M-sqrt(M*M-a*a))/(x[1]-M+sqrt(M*M-a*a))));
	double phi2 = x[3] - a/2/sqrt(M*M-a*a)*logr;
	double r = x[1];
	double th = x[2];
	double delta = r*r - 2*M*r + a*a;

	switch(coord)
	{
	case COORD_MAIN:
		WRITE_COMPONENT(jac, IND(0,0), 1.0);
		WRITE_COMPONENT(jac, IND(1,1), 1.0);
		WRITE_COMPONENT(jac, IND(2,2), 1.0);
		WRITE_COMPONENT(jac, IND(3,3), 1.0);
		break;
	case COORD_SPHERICAL:
		WRITE_COMPONENT(jac, IND(0,0), 1.0);
		WRITE_COMPONENT(jac, IND(0,1), -(r*r+a*a)/(r*r-2*M*r+a*a));
		WRITE_COMPONENT(jac, IND(1,1), 1.0);
		WRITE_COMPONENT(jac, IND(2,2), 1.0);
		WRITE_COMPONENT(jac, IND(3,1), -a/(r*r-2*M*r+a*a));
		WRITE_COMPONENT(jac, IND(3,3), 1.0);
		break;
	case COORD_CARTESIAN:
		WRITE_COMPONENT(jac, IND(0,0), 1.0);
		WRITE_COMPONENT(jac, IND(0,1), -(r*r+a*a)/(r*r-2*M*r+a*a));
		WRITE_COMPONENT(jac, IND(1,1), sin(th)*cos(phi2)+a/delta*r*sin(th)*sin(phi2));
		WRITE_COMPONENT(jac, IND(1,2), r*cos(th)*cos(phi2));
		WRITE_COMPONENT(jac, IND(1,3), -r*sin(th)*sin(phi2));
		WRITE_COMPONENT(jac, IND(2,1), sin(th)*sin(phi2)-a/delta*r*sin(th)*cos(phi2));
		WRITE_COMPONENT(jac, IND(2,2), r*cos(th)*sin(phi2));
		WRITE_COMPONENT(jac, IND(2,3), r*sin(th)*cos(phi2));
		WRITE_COMPONENT(jac, IND(3,1), cos(th));
		WRITE_COMPONENT(jac, IND(3,2), -r*sin(th));
		break;
	}

	return jac;
}
