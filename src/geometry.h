#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include "tensor.h"
#include "vector3.h"
#include <math.h>

//Metrics

#define M_SCHWARZSCHILD 0
#define M_KERR_EF 1

#define COORD_MAIN 0
#define COORD_SPHERICAL 1
#define COORD_CARTESIAN 2

class Geometry
{
	double M,a;
	int metric;

	//"Schwarzschild" packet
	Tensor g_schw(vector4);  //component of a metric
	Tensor g_inv_schw(vector4);
	Tensor Gamma_schw(vector4);
	vector4 convertFrom_schw(vector4, int);
	vector4 convertTo_schw(vector4, int);
	Tensor Jacobian_schw(int, vector4); //coords, position
	Tensor invJacobian_schw(int, vector4);

	//"Kerr" (Eddington-Finkelstein) packet
	Tensor g_kerr_ef(vector4);  //component of a metric
	Tensor g_inv_kerr_ef(vector4);
	Tensor Gamma_kerr_ef(vector4);
	vector4 convertFrom_kerr_ef(vector4, int);
	vector4 convertTo_kerr_ef(vector4, int);
	Tensor Jacobian_kerr_ef(int, vector4); //coords, position
	Tensor invJacobian_kerr_ef(int, vector4);

public:
	Geometry();
	Geometry(int);
	~Geometry();

	int setMetric(int);

	double getM();
	void setM(double);
	double getA();
	void setA(double);

	Tensor g(vector4);
	Tensor g_inv(vector4);
	Tensor Gamma(vector4);
	vector4 convertFrom(vector4, int);
	vector4 convertTo(vector4, int);
	Tensor Jacobian(int, vector4);		//d(main)/d(other)
	Tensor invJacobian(int, vector4);	//d(other)/d(main)

	double g_v(vector4, vector4, vector4);  //dot product
	double g_inv_v(vector4, vector4, vector4);  //dot product
	vector4 Gamma_v(vector4, vector4, vector4);
};

#endif
