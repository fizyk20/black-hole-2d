#include "vector3.h"

vector3::vector3()
{
	v[0] = v[1] = v[2] = 0.0;
}

vector3::vector3(double a, double b, double c)
{
	v[0] = a;
	v[1] = b;
	v[2] = c;
}

vector3::~vector3()
{
}

vector3& vector3::operator=(const vector3& arg)
{
	v[0] = arg.v[0];
	v[1] = arg.v[1];
	v[2] = arg.v[2];

	return (*this);
}

vector3 vector3::operator+(const vector3 arg)
{
	vector3 result(v[0]+arg.v[0], v[1]+arg.v[1], v[2]+arg.v[2]);
	return result;
}

vector3 vector3::operator+=(const vector3 arg)
{
	v[0] += arg.v[0];
	v[1] += arg.v[1];
	v[2] += arg.v[2];
	return (*this);
}

vector3 vector3::operator-(const vector3 arg)
{
	vector3 result(v[0]-arg.v[0], v[1]-arg.v[1], v[2]-arg.v[2]);
	return result;
}

vector3 vector3::operator-=(const vector3 arg)
{
	v[0] -= arg.v[0];
	v[1] -= arg.v[1];
	v[2] -= arg.v[2];

	return (*this);
}

vector3 vector3::operator*(const double arg)
{
	vector3 result(v[0]*arg, v[1]*arg, v[2]*arg);
	return result;
}

vector3 vector3::operator*=(const double arg)
{
	v[0] *= arg;
	v[1] *= arg;
	v[2] *= arg;

	return (*this);
}

vector3 vector3::operator/(const double arg)
{
	vector3 result(v[0]/arg, v[1]/arg, v[2]/arg);
	return result;
}

vector3 vector3::operator/=(const double arg)
{
	v[0] /= arg;
	v[1] /= arg;
	v[2] /= arg;

	return (*this);
}

vector3 operator*(const double a, const vector3 v)
{
	vector3 res(v.v[0]*a, v.v[1]*a, v.v[2]*a);
	return res;
}

double& vector3::operator[](int i)
{
	if(i<0) return v[0];
	if(i>2) return v[2];
	return v[i];
}
