#ifndef __TENSOR_H__
#define __TENSOR_H__

#include <vector>
#include "vector4.h"
using namespace std;

#define DIM 4   //dimension
#define IND(args...) { args }
#define WRITE_COMPONENT(t,x,y) { int __comp[] = x; \
	(t).component(__comp) = (y); }
#define READ_COMPONENT(y,t,x) { int __comp[] = x; \
	y = (t).component(__comp); }

class Tensor
{
	int nInd;
	int n;
	bool* bInd;
	double* components;

	int toComponent(int[]) const;
	vector<int> toIndex(int) const;
public:
	Tensor();
	Tensor(int n_Ind, bool b_Ind[]);
	Tensor(const Tensor&);
	~Tensor();

	bool isValid();

	Tensor& operator=(const Tensor&);

	Tensor operator+(const Tensor);
	Tensor operator+=(const Tensor);
	Tensor operator-(const Tensor);
	Tensor operator-=(const Tensor);
	Tensor operator*(const Tensor);
	Tensor operator*=(const Tensor);
	Tensor operator*(const double);
	Tensor operator*=(const double);
	friend Tensor operator*(const double, const Tensor);
	Tensor contracted(int, int);
	void contract(int, int);

	Tensor operator[](const int);
	//Tensor at(const int, int index=0);
	double& component(int[]);

	friend double toDouble(Tensor);
	friend vector4 toVector4(Tensor);
	friend Tensor Vector(vector4);
	friend Tensor Covector(vector4);
};

double toDouble(Tensor);
vector4 toVector4(Tensor);
Tensor Vector(vector4);
Tensor Covector(vector4);

#endif
