#pragma once

#include "Includes.h"

enum class TypeMatrix
{
	Full,
	Prof,
	Tape
};

template<class T>
class matrixA
{
public:

	vector<vector<T>> A;
	vector<T> al, au, di;
	vector<size_t> ia;

	vector<T> F;
	vector<T> X;
	matrixA();
	matrixA(size_t n);
	~matrixA();

	bool SaveFile();
	bool HighElemGausStraight();
	bool HighElemGausBack();
	bool Prof2Full();

	void PushBack(vector<T> Str);

	bool OpenFile(TypeMatrix typefile);
private:
};

template<class T>
class matrixLDU
{
public:
	vector<T> L;
	vector<T> D;
	vector<T> U;

	vector<size_t> ia;
	vector<vector<T>> A;

	vector<T> X;
	vector<T> Y;
	vector<T> Z;

	vector<T> F;

	matrixLDU();
	~matrixLDU();

	bool A2LDU(matrixA<T> A, TypeMatrix Type);
	bool LFX();
	bool SaveFile();
	bool SaveLDU(TypeMatrix Type);
private:
};