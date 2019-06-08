#pragma once

#include "stdafx.h"

enum class TypeMatrix
{
	Full,
	Prof,
	Tape
};

enum {
	NoErrors = 0,
	ErrorInputFile = 1,
	ErrorOutputFile = 2,
	MatrixNotCount = 3
};

template<class T>
class matrixA
{
public:

	vector<vector<T>> A;	//������� � ������ �������
	vector<T> al, au, di;	//������� � ���������� �������
	vector<size_t> ia;		//������� � ���������� �������

	vector<T> F;			//������ ������ ������� ���� Ax = F
	vector<T> X;			//��������� ���� 
	size_t Error;			//���������� �� �������

	matrixA();
	~matrixA();				

	void SaveFile();				//���������� ������� ���� � ���������� ��������
	void CountUpTrianMatrix();		//��������� ������� ����������� �������
	void GausBack();				//�������� ����� ������
	void Prof2Full();			    //������� �� ����������� ������� � ������

	void OpenFiles(TypeMatrix, char *);    //�������� ������
	void OpenFiles(TypeMatrix, char *, int k); 
private:
	double CompareNumber;
	size_t n; //����������� �������

	template <class U>
	bool OpenFile(const char name[], char *path, vector<U> &a, size_t len);
};

template<class T>
class matrixLDU
{
public:
	matrixLDU();
	~matrixLDU();

	void CountLDU(vector<T> &L, vector<T> &D, vector<T> &U, vector<size_t> &ia, TypeMatrix type, size_t &Error);//������� �� ����������� ������� matrixA � ���������� ������ LDU
	void CountX(vector<T> &L, vector<T> &D, vector<T> &U, vector<T> &F, vector<T> &X, vector<size_t> & ia);		//����� ������� ��������� LDUx = F
	void SaveFile(vector<T> &);																	//���������� ���������� ���� LDU x = F � ����
	void SaveLDU(vector<T> &L, vector<T> &D, vector<T> &U, TypeMatrix Type);						//���������� LDU � ����
private:
	double CompareNumber;
	size_t n; //����������� �������
	size_t m; //���������� �������������� ��������� � �������
	size_t *Error; //���������� �� �������
};