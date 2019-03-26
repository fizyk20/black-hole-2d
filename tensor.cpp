#include "tensor.h"

Tensor::Tensor()
{
	nInd = -1;
	bInd = NULL;
	components = NULL;
}

Tensor::Tensor(int n_Ind, bool b_Ind[])
{
	int i;
	nInd = n_Ind;
	bInd = new bool[nInd];
	for(i=0; i<nInd; i++) bInd[i] = b_Ind[i];
	n = 1;
	for(i=0; i<nInd; i++) n *= DIM;
	components = new double[n];
	for(i=0; i<n; i++) components[i] = 0.0;
}

Tensor::Tensor(const Tensor& t)
{
	int i;
	nInd = t.nInd;
	bInd = new bool[nInd];
	n = t.n;
	components = new double[n];

	for(i=0; i<nInd; i++)
		bInd[i] = t.bInd[i];
	for(i=0; i<n; i++)
		components[i] = t.components[i];
}

Tensor::~Tensor()
{
	if(components != NULL) delete[] components;
	if(bInd != NULL) delete[] bInd;
}

bool Tensor::isValid()
{
	return (components == NULL);
}

/***************************************************************************/

int Tensor::toComponent(int ind[]) const
{
	int a = 1;
	int res = 0;
	int i;
	for(i=0; i<nInd; i++)
	{
		if(ind[i] < 0 || ind[i] > DIM) return -1;
		res += a*ind[i];
		a *= DIM;
	}
	return res;
}

vector<int> Tensor::toIndex(int x) const
{
	vector<int> res;
	if(x < 0 || x >= n) return res;

	int i;
	for(i=0; i<nInd; i++)
	{
		res.push_back(x % DIM);
		x /= DIM;
	}
	return res;
}

/***************************************************************************/

Tensor& Tensor::operator=(const Tensor& t)
{
	if(components != NULL) delete[] components;
	if(bInd != NULL) delete[] bInd;

	int i;
	nInd = t.nInd;
	bInd = new bool[nInd];
	n = t.n;
	components = new double[n];

	for(i=0; i<nInd; i++)
		bInd[i] = t.bInd[i];
	for(i=0; i<n; i++)
		components[i] = t.components[i];

	return (*this);
}

Tensor Tensor::operator+(const Tensor t)
{
	int i;
	if(nInd != t.nInd) return Tensor();
	for(i=0; i<nInd; i++)
		if(bInd[i] != t.bInd[i]) return Tensor();

	Tensor res(nInd, bInd);
	for(i=0; i<n; i++)
		res.components[i] = components[i] + t.components[i];

	return res;
}

Tensor Tensor::operator+=(const Tensor t)
{
	int i;
	if(nInd != t.nInd) return Tensor();
	for(i=0; i<nInd; i++)
		if(bInd[i] != t.bInd[i]) return Tensor();

	for(i=0; i<n; i++)
		components[i] += t.components[i];

	return (*this);
}

Tensor Tensor::operator-(const Tensor t)
{
	int i;
	if(nInd != t.nInd) return Tensor();
	for(i=0; i<nInd; i++)
		if(bInd[i] != t.bInd[i]) return Tensor();

	Tensor res(nInd, bInd);
	for(i=0; i<n; i++)
		res.components[i] = components[i] - t.components[i];

	return res;
}

Tensor Tensor::operator-=(const Tensor t)
{
	int i;
	if(nInd != t.nInd) return Tensor();
	for(i=0; i<nInd; i++)
		if(bInd[i] != t.bInd[i]) return Tensor();

	for(i=0; i<n; i++)
		components[i] -= t.components[i];

	return (*this);
}

Tensor Tensor::operator*(const Tensor t)
{
	int i,j;

	int n_Ind = nInd + t.nInd;
	bool* b_Ind = new bool[n_Ind];
	for(i=0; i<n_Ind; i++)
		if(i<nInd) b_Ind[i] = bInd[i];
		else b_Ind[i] = t.bInd[i-nInd];

	Tensor res(n_Ind, b_Ind);

	delete[] b_Ind;

	int *ind1, *ind2;
	ind1 = new int[nInd];
	ind2 = new int[t.nInd];

	for(i=0; i<res.n; i++)
	{
		vector<int> ind = res.toIndex(i);
		for(j=0; j<res.nInd; j++)
			if(j<nInd) ind1[j] = ind[j];
			else ind2[j-nInd] = ind[j];
		res.components[i] = components[toComponent(ind1)]*t.components[t.toComponent(ind2)];
	}

	delete[] ind2;
	delete[] ind1;

	return res;
}

Tensor Tensor::operator*=(const Tensor t)
{
	(*this) = (*this)*t;
	return (*this);
}

Tensor Tensor::operator*(const double a)
{
	Tensor res = *this;
	for(int i=0; i<res.n; i++)
		res.components[i] *= a;
	return res;
}

Tensor Tensor::operator*=(const double a)
{
	for(int i=0; i<n; i++)
		components[i] *= a;
	return (*this);
}

Tensor operator*(const double a, const Tensor t)
{
	Tensor res = t;
	for(int i=0; i<res.n; i++)
		res.components[i] *= a;
	return res;
}

Tensor Tensor::contracted(int a, int b)
{
	if(a >= nInd || b >= nInd || a<0 || b<0) return Tensor();
	if(!(bInd[a] || bInd[b]) || (bInd[a] && bInd[b])) return Tensor();  //if not one upper and one lower index

	bool *b_Ind = new bool[nInd - 2];
	int i, j, k, l, sum;
	j = 0;
	for(i=0; i<nInd; i++)
	{
		if(i == a || i == b) continue;
		b_Ind[j] = bInd[i];
		j++;
	}

	Tensor res(nInd-2, b_Ind);
	delete[] b_Ind;

	int* ind1 = new int[nInd];

	for(i=0; i<res.n; i++)
	{
		vector<int> ind = res.toIndex(i);
		l = 0;
		for(k = 0; k < nInd; k++)
		{
			if(k == a || k == b) l++;
			else ind1[k] = ind[k-l];
		}

		res.components[i] = 0;
		for(sum = 0; sum < DIM; sum++)
		{
			ind1[a] = ind1[b] = sum;
			res.components[i] += components[toComponent(ind1)];
		}
	}

	delete[] ind1;

	return res;
}

void Tensor::contract(int a, int b)
{
	(*this) = this->contracted(a,b);
}

Tensor Tensor::operator[](const int a)
{
	if(a<0 || a>DIM) return Tensor();

	Tensor res(nInd-1, &(bInd[1]));
	int* ind = new int[nInd];
	int i,j;
	ind[0] = a;

	for(i=0; i<res.n; i++)
	{
		vector<int> ind1 = res.toIndex(i);
		for(j=1; j<nInd; j++) ind[j] = ind1[j-1];
		res.components[i] = components[toComponent(ind)];
	}

	delete[] ind;

	return res;
}

double& Tensor::component(int ind[])
{
	return components[toComponent(ind)];
}

/***************************************************************************/

double toDouble(Tensor t)
{
	if(t.nInd != 0) return 0.0;
	return t.components[0];
}

vector4 toVector4(Tensor t)
{
	if(t.nInd != 1) return vector4();
	return vector4(t.components[0], t.components[1], t.components[2], t.components[3]);
}

Tensor Vector(vector4 v)
{
	bool b[] = { true };
	Tensor res(1, b);
	for(int i=0; i<4; i++)
		res.components[i] = v[i];
	return res;
}

Tensor Covector(vector4 v)
{
	bool b[] = { false };
	Tensor res(1, b);
	for(int i=0; i<4; i++)
		res.components[i] = v[i];
	return res;
}
