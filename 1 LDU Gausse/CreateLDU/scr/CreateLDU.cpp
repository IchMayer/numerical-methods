#include "Includes.h"
#include "matrix.cpp"

#define type double

void Ori(vector<double> &D)
{
	vector<double> B;
	B.assign(D.begin(), D.end());
	B[2] = 3;
}

int main()
{
	vector<double> D(3, 2);
	Ori(D);
	if (0)//LDU
	{
		matrixA<type> A;
		A.OpenFile(TypeMatrix::Prof);
		matrixLDU<type> LDU;
		if (!LDU.A2LDU(A, TypeMatrix::Prof))
			return 0;
		LDU.SaveLDU(TypeMatrix::Prof);
		LDU.LFX();
		return 0;
	}
	else//GAUS
	{
		matrixA<type> A;
		A.OpenFile(TypeMatrix::Prof);
		A.Prof2Full();
		A.HighElemGausStraight();
		A.HighElemGausBack();
		return 0;
	}
}
