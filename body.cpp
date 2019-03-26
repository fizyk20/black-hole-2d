#include "body.h"
#include <math.h>
#include <QMessageBox>

CBody::CBody(vector4 init_pos, vector4 init_vel4, Geometry* g1)
{
	g = g1;

	u[0] = init_pos;
	u[1] = init_vel4;
	generateBasis();

	F = vector3();
	angVel = vector3();
	bool a1[] = {true, false};
	omega = Tensor(2, a1);
	tau=0;
}

CBody::~CBody()
{
}

void CBody::propagate(double dt, bool board) //RK4
{
	double h;
	if(board) h = dt;
	else
	{
		Tensor invjac = g->invJacobian(COORD_SPHERICAL,u[0]);
		vector4 u2 = toVector4((invjac*Vector(u[1])).contracted(1,2));
		h = dt/u2[0];
	}
	tau += h;

	int i,j;

	vector4 k[4][5];
	for(j=0; j<4; j++)
	{
		//generate k[j]
		vector4 new_y[5];
		for(i=0; i<5; i++)
		{
			new_y[i] = u[i];
			if(j>0 && j<3) new_y[i] += k[j-1][i]*h*0.5;
			if(j==3) new_y[i] += k[j-1][i]*h;
		}

		k[j][0] = new_y[1];
		constructOmega(new_y);
		for(i=1; i<5; i++)
		{
			vector4 a;
			a[i-1] = 1.0;
			vector4 f = toVector4((omega*Vector(a)).contracted(1,2));
			k[j][i] = f - g->Gamma_v(new_y[1],new_y[i],new_y[0]); //vi[mu]' + Gamma[mu,nu,sigma]u[nu]vi[sigma] = omega[mu,i];
		}
	}

	double coeff[4] = {1.0, 2.0, 2.0, 1.0};
	for(j=0; j<4; j++)
		for(i=0; i<5; i++)
			u[i] += coeff[j]*k[j][i]*h/6.0;

	//normalize phi=x[2]
	while(u[0][3]>M_PI) u[0][3]-= 2*M_PI;
	while(u[0][3]<-M_PI) u[0][3]+= 2*M_PI;
	if(u[0][2]<0) u[0][2] = -u[0][2];
	if(u[0][2]>M_PI) u[0][2] = 2*M_PI - u[0][2];

	orthonormalize();
}

vector4 CBody::getVector(int i)
{
	if(i<0 || i>4) return vector4();
	return u[i];
}

void CBody::setVector(int i, vector4 v)
{
	if(i<0 || i>4) return;
	u[i]=v;
	orthonormalize();
}

vector3 CBody::getCartesianVel()
{
	vector4 v_cart = toVector4((g->invJacobian(COORD_CARTESIAN,u[0])*Vector(u[1])).contracted(1,2));
	vector3 res;
	double gxx, gyy, gzz;
	Tensor jac = g->Jacobian(COORD_CARTESIAN, u[0]);
	Tensor g_cart = (g->g(u[0])*jac*jac).contracted(0,2).contracted(0,2);
	READ_COMPONENT(gxx, g_cart, IND(1,1));
	READ_COMPONENT(gyy, g_cart, IND(2,2));
	READ_COMPONENT(gzz, g_cart, IND(3,3));
	res[0] = v_cart[1]/v_cart[0]*sqrt(gxx);
	res[1] = v_cart[2]/v_cart[0]*sqrt(gyy);
	res[2] = v_cart[3]/v_cart[0]*sqrt(gzz);
	return res;
}

vector3 CBody::getCartesianDir(vector4 v)
{
	vector3 res;
	Tensor invjac = g->invJacobian(COORD_CARTESIAN, u[0]);
	Tensor jac = g->Jacobian(COORD_CARTESIAN, u[0]);
	Tensor gc = (g->g(u[0])*jac*jac).contracted(0,2).contracted(0,2);
	vector4 dt = vector4(1,0,0,0);
	vector4 v2 = toVector4((invjac*Vector(v)).contracted(1,2));

	v2 -= dt*toDouble((gc*Vector(v2)*Vector(dt)).contracted(0,2).contracted(0,1));
	double nv2 = sqrt(v2[1]*v2[1]+v2[2]*v2[2]+v2[3]+v2[3]);
	res[0] = v2[1]/nv2;
	res[1] = v2[2]/nv2;
	res[2] = v2[3]/nv2;
	return res;
}

double CBody::getTau()
{
	return tau;
}

void CBody::setTau(double t)
{
	tau = t;
}

vector3 CBody::getForce()
{
	return F;
}

void CBody::applyForce(vector3 f)
{
	F += f;
}

void CBody::setForce(vector3 f)
{
	F = f;
}

void CBody::setAngVel(vector3 w)
{
	angVel = w;
}

Tensor CBody::conversionMatrix(vector4 new_y[])
{
	bool a1[] = {true, false};
	Tensor A(2,a1);

	int i,j;
	for(i=0; i<4; i++)
		for(j=0; j<4; j++)
			WRITE_COMPONENT(A, IND(i,j), new_y[j+1][i]);
	return A;
}

void CBody::constructOmega(vector4 new_y[])
{
	Tensor A = conversionMatrix(new_y);
	bool a1[] = {true, false};
	Tensor omega1(2, a1);

	WRITE_COMPONENT(omega1, IND(0,1), F[0]);
	WRITE_COMPONENT(omega1, IND(1,0), F[0]);
	WRITE_COMPONENT(omega1, IND(0,2), F[1]);
	WRITE_COMPONENT(omega1, IND(2,0), F[1]);
	WRITE_COMPONENT(omega1, IND(0,3), F[2]);
	WRITE_COMPONENT(omega1, IND(3,0), F[2]);

	WRITE_COMPONENT(omega1, IND(1,2), -angVel[2]);
	WRITE_COMPONENT(omega1, IND(2,1), angVel[2]);
	WRITE_COMPONENT(omega1, IND(1,3), angVel[1]);
	WRITE_COMPONENT(omega1, IND(3,1), -angVel[1]);
	WRITE_COMPONENT(omega1, IND(2,3), -angVel[0]);
	WRITE_COMPONENT(omega1, IND(3,2), angVel[0]);

	omega = (A*omega1).contracted(1,2);		//omega(mu,sigma) = A(mu,nu)*omega1(nu, sigma);
}

void CBody::orthonormalize()
{
	int i,j;
	for(i=1; i<5; i++)
	{
		for(j=1; j<i; j++)
			u[i] -= g->g_v(u[i],u[j],u[0])*u[j]/(g->g_v(u[j],u[j],u[0]));
		u[i] /= sqrt(fabs(g->g_v(u[i],u[i],u[0])));
	}
}

void CBody::generateBasis()
{
	//we assume that position and four-velocity are already set
	u[2] = toVector4((g->Jacobian(COORD_SPHERICAL,u[0])*Vector(vector4(0,1,0,0))).contracted(1,2));
	u[3] = toVector4((g->Jacobian(COORD_SPHERICAL,u[0])*Vector(vector4(0,0,0,1))).contracted(1,2));
	u[4] = toVector4((g->Jacobian(COORD_SPHERICAL,u[0])*Vector(vector4(0,0,1,0))).contracted(1,2));
	orthonormalize();
}
