#pragma once
#include "includs.h"
#include "matrix (2).cpp"

template<class T>
class Matrix
{
public:
	Matrix() { maxiter = 10;  d = k = n = m = Error = 0; e = 0.001; w = 1; }
	~Matrix(){}

	//������� A, ���������� ��� ���������� ������
	T matrix(size_t i,size_t j);
	//� ��������� Ax = F
	vector<T> result;
	//����� ����� ��������� �� = F, (F - f) - ������
	vector<T> F, f;
	//������� ���������� (����� ������)
	vector<int> Index;
	size_t k;
	//������ � �����
	T e, w;
	//������������ ���������� ��������
	size_t maxiter;

	//��������� ���������� �� ������
	size_t GetLastError();

	void GetResult();
	void SetMatrix();

	//�������� ������
	void OpenMatrix();
	//���������� ����������
	void SaveResult();

	//����� �����
	void JacobiMethod();
	//������ ������
	void GaussSeidelMethod();
	//����� ���������� � ���������� ������ �������� NBlock
	void BlockRelax(int Nblock);
private:
	//������������
	vector<T> Fractorization(int Nblock);
	//����� ������� (L2)
	inline T NormVector(vector<T>);
	//��������� ������� matrix �� ������ X
	vector<T> MultMatrixAandVector(vector<T> X);
	//�������� �������� a � b
	vector<T> ResidualVectors(vector<T> a, vector<T> b);
	//���������� ���������� ��� ������ ������
	void Gaussnewvector(vector<T> &f, T x, int i);
	//���������� ���������� ��� ������ ����������
	void Blocknewvector(vector<T> &f, T x, int i, int Max);
	//�������� �� �����
	inline bool CheckEnd(); 
	//������� ��������� �������������
	T RelaxParam(T x2, T x1);

#ifdef _DEBUG
	T a[100000]; //������������ ������ 1000 �� 100
	const int MAX = 100000;
#else
	T a[100000]; //������������ ������ 1000 �� 100
	const int MAX = 100000;
#endif // DEBUG
	vector<T> buf;
	size_t n, m;
	size_t Error;
	size_t d;
};
/*Errors:
  0 - not error
  1 - matrix overflow
  2 - end method
*/
