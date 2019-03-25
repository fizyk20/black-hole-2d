#ifndef __BODY_H__
#define __BODY_H__

#include "geometry.h"

#define VEC_POS 0
#define VEC_4VEL 1
#define VEC_RIGHT 2
#define VEC_FRONT 3
#define VEC_UP 4

class CBody
{
	Geometry* g;
	double tau;
	
	vector4 u[5]; //position, 4-velocity, local x (right), local y (front), local z (up)
	
	vector3 F;  //local force
	vector3 angVel; //local angular velocity;
	
	Tensor omega;  //describes the acceleration

	void constructOmega(vector4[]);
	Tensor conversionMatrix(vector4[]);
	void orthonormalize();
public:
	CBody(vector4, vector4, Geometry*);
	~CBody();
	void propagate(double, bool board = false);

	vector4 getVector(int);
	void setVector(int,vector4);
	void generateBasis();
	
	vector3 getCartesianVel();
	vector3 getCartesianDir(vector4);
	
	double getTau();
	void setTau(double);
	
	vector3 getForce();

	void applyForce(vector3); //prograde, right
	void setForce(vector3); //prograde, right
	void setAngVel(vector3);
};

#endif