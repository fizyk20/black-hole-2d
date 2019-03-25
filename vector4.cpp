#include "vector4.h"

vector4::vector4()
{
	v[0] = v[1] = v[2] = v[3] = 0.0;
}

vector4::vector4(double a, double b, double c, double d)
{
	v[0] = a;
	v[1] = b;
	v[2] = c;
	v[3] = d;
}

vector4::~vector4()
{
}

vector4& vector4::operator=(const vector4& arg)
{
	v[0] = arg.v[0];
	v[1] = arg.v[1];
	v[2] = arg.v[2];
	v[3] = arg.v[3];
	
	return (*this);
}

vector4 vector4::operator+(const vector4 arg)
{
	vector4 result(v[0]+arg.v[0], v[1]+arg.v[1], v[2]+arg.v[2], v[3]+arg.v[3]);
	return result;
}

vector4 vector4::operator+=(const vector4 arg)
{
	v[0] += arg.v[0];
	v[1] += arg.v[1];
	v[2] += arg.v[2];
	v[3] += arg.v[3];
	
	return (*this);
}

vector4 vector4::operator-(const vector4 arg)
{
	vector4 result(v[0]-arg.v[0], v[1]-arg.v[1], v[2]-arg.v[2], v[3]-arg.v[3]);
	return result;
}

vector4 vector4::operator-=(const vector4 arg)
{
	v[0] -= arg.v[0];
	v[1] -= arg.v[1];
	v[2] -= arg.v[2];
	v[3] -= arg.v[3];
	
	return (*this);
}

vector4 vector4::operator*(const double arg)
{
	vector4 result(v[0]*arg, v[1]*arg, v[2]*arg, v[3]*arg);
	return result;
}

vector4 vector4::operator*=(const double arg)
{
	v[0] *= arg;
	v[1] *= arg;
	v[2] *= arg;
	v[3] *= arg;
	
	return (*this);
}

vector4 vector4::operator/(const double arg)
{
	vector4 result(v[0]/arg, v[1]/arg, v[2]/arg, v[3]/arg);
	return result;
}

vector4 vector4::operator/=(const double arg)
{
	v[0] /= arg;
	v[1] /= arg;
	v[2] /= arg;
	v[3] /= arg;
	
	return (*this);
}

vector4 operator*(const double a, const vector4 v)
{
	vector4 res(v.v[0]*a, v.v[1]*a, v.v[2]*a, v.v[3]*a);
	return res;
}

double& vector4::operator[](int i)
{
	if(i<0) return v[0];
	if(i>3) return v[3];
	return v[i];
}
