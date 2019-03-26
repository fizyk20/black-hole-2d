#ifndef __VECTOR4__
#define __VECTOR4__

class vector4
{
	double v[4];
public:
	vector4();
	vector4(double, double, double, double);
	~vector4();

	vector4& operator=(const vector4&);
	vector4 operator+(const vector4);
	vector4 operator+=(const vector4);
	vector4 operator-(const vector4);
	vector4 operator-=(const vector4);
	vector4 operator*(const double);
	vector4 operator*=(const double);
	vector4 operator/(const double);
	vector4 operator/=(const double);
	friend vector4 operator*(const double, const vector4);

	double& operator[](int);
};

#endif
