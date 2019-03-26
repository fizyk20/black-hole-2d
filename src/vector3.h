#ifndef __VECTOR3__
#define __VECTOR3__

class vector3
{
	double v[3];
public:
	vector3();
	vector3(double, double, double);
	~vector3();

	vector3& operator=(const vector3&);
	vector3 operator+(const vector3);
	vector3 operator+=(const vector3);
	vector3 operator-(const vector3);
	vector3 operator-=(const vector3);
	vector3 operator*(const double);
	vector3 operator*=(const double);
	vector3 operator/(const double);
	vector3 operator/=(const double);
	friend vector3 operator*(const double, const vector3);

	double& operator[](int);
};

#endif
