#pragma once
#include "includes.h"

template<class T>

///<summary> Matrix - class for lab 3 4m</summary>
///<remarks>conjugate gradient metod and local - optimal scheme 4 asymmetric matrix</remarks>
class Matrix
{
public:
	Matrix() {}
	~Matrix() {}

	void OpenFile()
	{
		ifstream in("kuslau.txt");
		in >> N >> Maxiter >> e;
		in.close();

		ig.resize(N + 1);
		f.resize(N);
		di.resize(N);
		result.resize(N);

		in.open("ig.txt");
		for (size_t i = 0; i <= N; i++)
			in >> ig[i];
		in.close();

		M = ig[N];
		jg.resize(M);
		ggl.resize(M);
		ggu.resize(M);

		in.open("jg.txt");
		for (size_t i = 0; i < M; i++)
			in >> jg[i];
		in.close();

		in.open("ggl.txt");
		for (size_t i = 0; i < M; i++)
			in >> ggl[i];
		in.close();

		in.open("ggu.txt");
		for (size_t i = 0; i < M; i++)
			in >> ggu[i];
		in.close();

		in.open("di.txt");
		for (size_t i = 0; i < N; i++)
			in >> di[i];
		in.close();

		in.open("pr.txt");
		for (size_t i = 0; i < N; i++)
			in >> f[i];
		in.close();
	}

	//Set matrix dimension
	///<param name = "n"> New matrix dimension</param>
	void SetN(size_t n) { N = n; }
	//Get matrix dimention
	size_t n() { return N; }

	//Возращает резултат операций
	vector<T> Result() { return result; }
	//Возращает количество операций
	size_t K() { return k; }

	//Set max iteration
	///<param name = "Max"> New max ieration </param>
	void Setmaxiter(size_t Max) { Maxiter = Max; }
	//Get max iteration
	size_t maxiter() { return Maxiter; }

	///<param name = "i"> i = 0 без предобуславливония
	///i = 1 диагональное
	///i = 2 Холесского</param>
	void preconditioning(int i)
	{
		switch (i)
		{
			//Без предобуславливония
		case 0:

			break;
			//Диагональное
		case 1:
			//			factorizationLU(false);
			factorizationLLT();
			break;
			//Холесского
	//	case 2:
		//	factorizationLLT(true);
			//break;
			//LU
		case 3:
			factorizationLU(true);
			break;
		default:
			break;
		}
	}

	//i = 0 без предобуславливония
	//i = 1 диагональное
	//i = 2 Холесского</param>
	void preconditioning()
	{
		preconditioning(3);
	}

	void ConjugateGradient(int metod)
	{
		NormF = Norm(f);
		switch (metod)
		{
			//Без
		case 0:
		{
			vector<T> Best = result;
			f = MultTMatrixVector(f);
			r = Residual(f, MultTMatrixVector(MultMatrixVector(Best)));
			z = r;
			T rLastScolar = Scolar(r);
			T CoudA;
			T min = 0.;
			while (checkEnd(rLastScolar, CoudA))
			{
				vector<T> d = MultTMatrixVector(MultMatrixVector(z));
				a = rLastScolar / Scolar(d, z);
				Best = Sum(Best, MultVectorOnT(z, a));
				rLast = r;
				r = Residual(r, MultVectorOnT(d, a));
				T lol = Scolar(r);
				b = lol / rLastScolar;
				rLastScolar = lol;
				z = Sum(r, MultVectorOnT(z, b));
				if (min > CoudA || !min)
				{
					min = CoudA;
					result = Best;
				}
				k++;
			}
			break;
		}
		// Диаганальное
		case 1:
			//LLT
		case 2:
		{
			vector<T> d = Residual(f, MultMatrixVector(result)); //r = U^-t * A ^t * L^ -T * L-1(f - Ax)
			d = Down(L, D, d);// L d1 = d
			d = Up(L, D, d);	//L^t d1 = d
			d = MultTMatrixVector(d); // A^t * d = d1
			r = Down(L, D, d);	// U^t r = d
			z = r;

			x.resize(N);
#pragma region Умножение_x-_=_Ux
			for (size_t i = 0; i < N; i++)
			{
				for (size_t j = ig[i]; j < ig[i + 1]; j++)
					x[i] += L[j] * result[jg[j]];
				x[i] += D[i] * result[i];
			}
#pragma endregion

			k = 0;
			//xbest = x;

			vector<T> Best = x;
			T rLastScolar = Scolar(r);
			T CoudA;
			T min = 0.;

			while (checkEnd(rLastScolar, CoudA))
			{

				vector<T> d = Up(L, D, z);
				d = MultMatrixVector(d);
				d = Down(L, D, d);
				d = Up(L, D, d);
				d = MultTMatrixVector(d);
				d = Down(L, D, d);
				T s = Scolar(d, z);
				if (!s)
					break;
				a = rLastScolar / s;
				x = Sum(x, MultVectorOnT(z, a));
				rLast = r;
				r = Residual(r, MultVectorOnT(d, a));
				T rScolar = Scolar(r);
				b = rScolar / rLastScolar;
				z = Sum(r, MultVectorOnT(z, b));
				rLastScolar = rScolar;
				k++;
				if (min > CoudA || !min)
				{
					min = CoudA;
					Best = x;
				}
			}
			result = Up(L, D, Best);
			break;
		}
		//LU
		case 3:
		{
			vector<T> d = Residual(f, MultMatrixVector(result)); //r = U^-t * A ^t * L^ -T * L-1(f - Ax)
			d = Down(L, d);// L d1 = d
			d = Up(L, d);	//L^t d1 = d
			d = MultTMatrixVector(d); // A^t * d = d1
			r = Down(U, D, d);	// U^t r = d
			z = r;

			x.resize(N);
#pragma region Умножение_x-_=_Ux
			for (size_t i = 0; i < N; i++)
			{
				for (size_t j = ig[i]; j < ig[i + 1]; j++)
					x[i] += U[j] * result[jg[j]];
				x[i] += D[i] * result[i];
			}
#pragma endregion

			k = 0;
			xbest = x;

			vector<T> Best = x;
			T rLastScolar = Scolar(r);
			T Coud;
			T min = 0.;
			while (checkEnd(rLastScolar, Coud))
			{

				vector<T> d = Up(U, D, z);
				d = MultMatrixVector(d);
				d = Down(L, d);
				d = Up(L, d);
				d = MultTMatrixVector(d);
				d = Down(U, D, d);
				T s = Scolar(d, z);
				if (!s)
					break;
				a = rLastScolar / s;
				x = Sum(x, MultVectorOnT(z, a));
				r = Residual(r, MultVectorOnT(d, a));
				T rScolar = Scolar(r);
				b = rScolar / rLastScolar;
				rLastScolar = rScolar;
				z = Sum(r, MultVectorOnT(z, b));

				if (min > Coud || !min)
				{
					min = Coud;
					Best = x;
				}

				k++;
			}
			result = Up(U, D, Best);
			break;
		}
		default:
			break;
		}
	}
	void LocalOptimalScheme(int metod)
	{
		NormF = Norm(f);
		switch (metod)
		{
		case 0:
		{
			r = Residual(f, MultMatrixVector(result));
			z = r;
			p = MultMatrixVector(z);
			k = 0;
			T ScolarP = Scolar(p);
			T ScolarR = Scolar(r);
			//while (checkEnd(ScolarR, ScolarP, ScolarR))
			while (checkEnd(r))
			{
				if (!ScolarP)
					break;
				a = Scolar(p, r) / ScolarP;
				result = Sum(result, MultVectorOnT(z, a));
				r = Residual(r, MultVectorOnT(p, a));
				T gg = Scolar(r);
				vector<T> d = MultMatrixVector(r);
				b = -1.0 * Scolar(p, d) / ScolarP;
				z = Sum(r, MultVectorOnT(z, b));
				p = Sum(d, MultVectorOnT(p, b));
				ScolarP = Scolar(p);
				k++;
			}
			break;
		}
		//DLLT
		case 1:
			//LLT
		case 2:
		{
			r = Down(L, D, Residual(f, MultMatrixVector(result)));
			z = Up(L, D, r);
			p = Down(L, D, MultMatrixVector(z));

			k = 0;
			while (checkEnd(r))
			{
				T s = Scolar(p);
				if (!s)
					break;
				a = Scolar(p, r) / s;
				result = Sum(result, MultVectorOnT(z, a));
				r = Residual(r, MultVectorOnT(p, a));
				vector<T> d = Up(L, D, r);
				d = MultMatrixVector(d);
				d = Down(L, D, d);
				b = -1.0 * Scolar(p, d) / s;
				z = Sum(Up(L, D, r), MultVectorOnT(z, b));
				p = Sum(d, MultVectorOnT(p, b));
				k++;
			}
			break;
		}
		//LU
		case 3:
		{
			r = Down(L, Residual(f, MultMatrixVector(result)));
			z = Up(U, D, r);
			p = Down(L, MultMatrixVector(z));

			k = 0;
			while (checkEnd(r))
			{
				T s = Scolar(p);
				if (!s)
					break;
				a = Scolar(p, r) / s;
				result = Sum(result, MultVectorOnT(z, a));
				r = Residual(r, MultVectorOnT(p, a));
				vector<T> d = Up(U, D, r);
				d = MultMatrixVector(d);
				d = Down(L, d);
				b = -1.0 * Scolar(p, d) / s;
				z = Sum(Up(U, D, r), MultVectorOnT(z, b));
				p = Sum(d, MultVectorOnT(p, b));
				k++;
			}
			break;
		}
		default:
			break;
		}
	}

private:
	vector<T> f; // Массив правой части
	vector<T> di;	//Диагональные элементы матрицы А
	vector<T> ggu;	//Верхний треугольник матрицы А в разреженном формате
	vector<T> ggl;	//Нижний треугольний матрицы А в разреженном формате
	vector<T> L; // Матрица L
	vector<T> U; // НЕ диагональыне элементы  матрицы U
	vector<T> D; // Диагональные элементы матрицы U
	vector<size_t> ig;// Массив индексов
	vector<size_t> jg; // Другой массив индексов

	vector<T> result;	//Результат x
	vector<T> r, rLast, z, p, x, xbest;

	size_t N, M, Maxiter, k;
	T e, a, b, NormF;

	//Проверка окончания по невязки и по maxiter

	///<param name = "ScolarR"> Входной параметр. Сколяр вектора R</param>
	///<param name = "Return"> Выходной параметр. Невязка</param>
	//Проверка конца у МСГ
	bool checkEnd(T ScolarR, T &Return)
	{
		if (k == Maxiter)
			return false;
		Return = sqrt(ScolarR) / NormF;
		if (e > Return)
			return false;
		return true;
	}


	///<param name = "ScolarR"> Входной параметр. Сколяр вектора R</param>
	///<param name = "ScolarP"> Входной параметр. Сколяр вектора P</param>
	///<param name = "ReturnR"> Выходной параметр. Сколяр вектора R</param>
	//Проверка конца у ЛОС
	bool checkEnd(T ScolarR, T ScolarP, T &ReturnR)
	{
		if (k == Maxiter)
			return false;

		T R = sqrt(ScolarR) / NormF;
		T g = Scolar(r);
		if (e > R)
			return false;
		a = Scolar(p, r) / Scolar(p);
		T h = pow(a, 2) * ScolarP;
		ReturnR = ScolarR - h;
		return true;
	}

	bool checkEnd(vector<T> r)
	{
		if (k == Maxiter)
			return false;

		if (e > Norm(r) / NormF)
			return false;
		return true;
	}

	//Прямой ход
	vector<T> Down(vector<T> Matrix, vector<T> Diagonal, vector<T> R)
	{
		vector<T> x(N);
		for (size_t i = 0; i < N; i++)
		{
			T sum = 0;
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
				sum += Matrix[j] * x[jg[j]];
			x[i] = (R[i] - sum) / Diagonal[i];
		}
		return x;
	}
	vector<T> Down(vector<T> Matrix, vector<T> R)
	{
		vector<T> x(N);
		for (size_t i = 0; i < N; i++)
		{
			T sum = 0;
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
				sum += Matrix[j] * x[jg[j]];
			x[i] = R[i] - sum;
		}
		return x;
	}

	//Обратный ход
	vector<T> Up(vector<T> Matrix, vector<T> Diagonal, vector<T> R)
	{
		for (int i = N - 1; i >= 0; i--)
		{
			R[i] /= Diagonal[i];
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				size_t p = jg[j];
				R[p] -= Matrix[j] * R[i];
			}
		}
		return R;
	}
	vector<T> Up(vector<T> Matrix, vector<T> R)
	{
		for (int i = N - 1; i >= 0; i--)
		{
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				size_t p = jg[j];
				R[p] -= Matrix[j] * R[i];
			}
		}
		return R;
	}

	//Умноженик числа на вектор
	vector<T> MultVectorOnT(vector<T> a, T b)
	{
		for (size_t i = 0; i < a.size(); i++)
			a[i] *= b;
		return a;
	}

	//Mult matrix and vector f
	///<param name = "f"> Matrix be multtiplication on thix vector.</param>
	vector<T> MultMatrixVector(vector<T> f)
	{
		vector<T> v(N);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				v[i] += ggl[j] * f[jg[j]];
				v[jg[j]] += ggu[j] * f[i];
			}
			v[i] += di[i] * f[i];
		}
		return v;
	}
	vector<T> MultTMatrixVector(vector<T> f)
	{
		vector<T> v(N);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				v[i] += ggu[j] * f[jg[j]];
				v[jg[j]] += ggl[j] * f[i];
			}
			v[i] += di[i] * f[i];
		}
		return v;
	}

	//Residual a - b
	vector<T> Residual(vector<T> a, vector<T> b)
	{
		vector<T> v;
		v.resize(a.size());
		for (size_t i = 0; i < a.size(); i++)
			v[i] = a[i] - b[i];
		return v;
	}
	//Sum a + b
	vector<T> Sum(vector<T> a, vector<T> b)
	{
		vector<T> v;
		v.resize(a.size());
		for (size_t i = 0; i < a.size(); i++)
		{
			v[i] = a[i] + b[i];
		}
		return v;
	}

	//Норма вектора а
	inline T Norm(vector<T> a)
	{
		return sqrt(Scolar(a));
	}
	//Скалярное произведение 1 элемента
	inline T Scolar(vector<T> a)
	{
		T sum = 0;
		for (size_t i = 0; i < N; i++)
			sum += pow(a[i], 2);
		return sum;
	}
	//Скалярное произведение 2 элемента
	inline T Scolar(vector<T> a, vector<T> b)
	{
		T sum = 0;
		for (size_t i = 0; i < N; i++)
			sum += a[i] * b[i];
		return sum;
	}

	//factorizationLU is a method in the Matrix class. МБ даже работает
	///<param name = "h"> if this param true, then LU, else diagonal LU </param>
	void factorizationLU(bool h)
	{
		if (h)
		{
			L.resize(ggl.size());
			D.resize(di.size());
			U.resize(ggu.size());
			for (size_t i = 0; i < N; i++)
			{
				T sum = 0;
				for (size_t j = ig[i]; j < ig[i + 1]; j++)
				{
					T sum1 = 0, sum2 = 0;
					int jj = jg[j];
					for (size_t k = ig[i], k2 = ig[jj]; k < j && k2 < ig[jj + 1];)
					{
						int p1 = jg[k];
						int p2 = jg[k2];
						if (p1 == p2)
						{
							sum1 += L[k] * U[k2];
							sum2 += L[k2] * U[k];
							k++; k2++;
						}
						else
						{
							if (p1 < p2)
								k++;
							else
								k2++;
						}
					}
					L[j] = (ggl[j] - sum1) / D[jg[j]];
					U[j] = ggu[j] - sum2;
					sum += L[j] * U[j];
				}
				D[i] = di[i] - sum;
			}
		}
		else
		{
			L.resize(ggl.size());
			U.resize(ggu.size());
			D = di;
		}
	}
	//factorizationLU is a method in the Matrix class. МБ даже работает
	///<param name = "h"> if this param true, then Holisskigo, else diagonal Holisskigo </param>
	void factorizationLLT()
	{
		L.resize(ggl.size());
		D.resize(N);
		for (size_t i = 0; i < N; i++)
			D[i] = sqrt(di[i]);
	}
};

