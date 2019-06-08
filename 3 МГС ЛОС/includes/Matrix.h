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
			vector<T> Buf(N), Buf1(N), d(N), Best = result;
			r.resize(N);
			MultTMatrixVector(f, Buf);
			f = Buf;
			MultMatrixVector(Best, Buf);
			MultTMatrixVector(Buf, Buf1);
			Subtraction(f, Buf1, r);
			z = r;
			T rLastScolar = Scolar(r);
			T CoudA;
			T min = 0.;
			while (checkEnd(rLastScolar, CoudA))
			{
				MultMatrixVector(z, Buf);
				MultTMatrixVector(Buf, d);
				a = rLastScolar / Scolar(d, z);
				MultVectorOnT(z, a, Buf);
				Sum(Best, Buf, Buf1);
				Best = Buf1;
				rLast = r;
				MultVectorOnT(d, a, Buf);
				Subtraction(r, Buf, Buf1);
				r = Buf1;
				T lol = Scolar(r);
				b = lol / rLastScolar;
				rLastScolar = lol;
				MultVectorOnT(z, b, Buf);
				Sum(r, Buf, z);
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
			vector<T> Buf(N), Buf1(N), d(N);
			r.resize(N);
			MultMatrixVector(result, Buf);
			Subtraction(f, Buf, d); //r = U^-t * A ^t * L^ -T * L-1(f - Ax)
			Diagonal(D, d, Buf);// L d1 = d
			Diagonal(D, Buf, Buf1);	//L^t d1 = d
			MultTMatrixVector(Buf1, d); // A^t * d = d1
			Diagonal(D, d, r);	// U^t r = d
			z = r;
			
			x.resize(N);
#pragma region Умножение_x-_=_Ux
			for (size_t i = 0; i < N; i++)
			{
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

				Diagonal(D, z, Buf);
				MultMatrixVector(Buf, Buf1);
				Diagonal(D, Buf1, Buf);
				Diagonal(D, Buf, Buf1);
				MultTMatrixVector(Buf1, Buf);
				Diagonal(D, Buf, d);
				T s = Scolar(d, z);
				if (!s)
					break;
				a = rLastScolar / s;
				MultVectorOnT(z, a, Buf);
				Sum(x, Buf, Buf1);
				x = Buf1;
				rLast = r;
				MultVectorOnT(d, a, Buf);
				Subtraction(r, Buf, Buf1);
				r = Buf1;
				T rScolar = Scolar(r);
				b = rScolar / rLastScolar;
				MultVectorOnT(z, b, Buf);
				Sum(r, Buf, z);
				rLastScolar = rScolar;
				k++;
				if (min > CoudA || !min)
				{
					min = CoudA;
					Best = x;
				}
			}
			Diagonal(D, Best, result);
			break;
		}
		//LU
		case 3:
		{
			vector<T> Buf(N), Buf1(N), d(N);
			r.resize(N);
			MultMatrixVector(result, Buf);
			Subtraction(f, Buf, Buf1); //r = U^-t * A ^t * L^ -T * L-1(f - Ax)
			Down(L, Buf1, Buf);// L d1 = d
			Up(L, Buf, Buf1);	//L^t d1 = d
			MultTMatrixVector(Buf1, d); // A^t * d = d1
			Down(U, D, d, r);	// U^t r = d
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

				Up(U, D, z, Buf);
				MultMatrixVector(Buf, Buf1);
				Down(L, Buf1, Buf);
				Up(L, Buf, Buf1);
				MultTMatrixVector(Buf1, Buf);
				Down(U, D, Buf, d);
				T s = Scolar(d, z);
				if (!s)
					break;
				a = rLastScolar / s;
				MultVectorOnT(z, a, Buf);
				Sum(x, Buf, x);
				MultVectorOnT(d, a, Buf);
				Subtraction(r, Buf, Buf1);
				r = Buf1;
				T rScolar = Scolar(r);
				b = rScolar / rLastScolar;
				rLastScolar = rScolar;
				MultVectorOnT(z, b, Buf);
				Sum(r, Buf, z);

				if (min > Coud || !min)
				{
					min = Coud;
					Best = x;
				}

				k++;
			}
			Up(U, D, Best, result);
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
			vector<T> d(N), Buf(N), Buf1(N);
			r.resize(N);
			p.resize(N);
			MultMatrixVector(result, Buf);
			Subtraction(f, Buf, r);
			z = r;
			MultMatrixVector(z, p);
			k = 0;
			T ScolarP = Scolar(p);
			T ScolarR = Scolar(r);
			//while (checkEnd(ScolarR, ScolarP, ScolarR))
			while (checkEnd(r))
			{
				if (!ScolarP)
					break;
				a = Scolar(p, r) / ScolarP;
				MultVectorOnT(z, a, Buf);
				Sum(result, Buf, Buf1);
				result = Buf1;
				MultVectorOnT(p, a, Buf);
				Subtraction(r, Buf, Buf1);
				r = Buf1;
				T gg = Scolar(r);
				MultMatrixVector(r, d);
				b = -1.0 * Scolar(p, d) / ScolarP;
				MultVectorOnT(z, b, Buf);
				Sum(r, Buf, z);
				MultVectorOnT(p, b, Buf);
				Sum(d, Buf, p);
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
			vector<T> d(N), Buf(N), Buf1(N);
			r.resize(N);
			z.resize(N);
			p.resize(N);
			MultMatrixVector(result, Buf);
			Subtraction(f, Buf, Buf1);
			Diagonal(D, Buf1, r);
			Diagonal(D, r, z);
			MultMatrixVector(z, Buf);
			Diagonal(D, Buf, p);
			k = 0;
			while (checkEnd(r))
			{
				T s = Scolar(p);
				if (!s)
					break;
				a = Scolar(p, r) / s;
				MultVectorOnT(z, a, Buf);
				Sum(result, Buf, Buf1);
				result = Buf1;
				MultVectorOnT(p, a, Buf);
				Subtraction(r, Buf, Buf1);
				r = Buf1;
				Diagonal(D, r, Buf);
				MultMatrixVector(Buf, Buf1);
				Diagonal(D, Buf1, d);
				b = -1.0 * Scolar(p, d) / s;
				MultVectorOnT(z, b, Buf);
				Diagonal(D, r, Buf1);
				Sum(Buf1, Buf, z);
				MultVectorOnT(p, b, Buf);
				Sum(d, Buf, p);
				k++;
			}
			break;
		}
		//LU
		case 3:
		{
			vector<T> d(N), Buf(N), Buf1(N);
			r.resize(N);
			p.resize(N);
			MultMatrixVector(result, Buf);
			Subtraction(f, Buf, Buf1);
			Down(L, Buf1, r);
			Up(U, D, r, z);
			MultMatrixVector(z, Buf);
			Down(L, Buf, p);
			k = 0;
			while (checkEnd(r))
			{
				T s = Scolar(p);
				if (!s)
					break;
				a = Scolar(p, r) / s;
				MultVectorOnT(z, a, Buf);
				Sum(result, Buf, Buf1);
				result = Buf1;
				MultVectorOnT(p, a, Buf);
				Subtraction(r, Buf, Buf1);
				r = Buf1;
				Up(U, D, r, Buf);
				MultMatrixVector(Buf, Buf1);
				Down(L, Buf1, d);
				b = -1.0 * Scolar(p, d) / s;
				Up(U, D, r, Buf);
				MultVectorOnT(z, b, Buf1);
				Sum(Buf, Buf1, z);
				MultVectorOnT(p, b, Buf);
				Sum(d, Buf, p);
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
	void Down(vector<T> Matrix, vector<T> Diagonal, vector<T> R, vector<T> &x)
	{
		for (size_t i = 0; i < N; i++)
		{
			T sum = 0;
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
				sum += Matrix[j] * x[jg[j]];
			x[i] = (R[i] - sum) / Diagonal[i];
		}
	}
	void Down(vector<T> Matrix, vector<T> R, vector<T> &x)
	{
		for (size_t i = 0; i < N; i++)
		{
			T sum = 0;
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
				sum += Matrix[j] * x[jg[j]];
			x[i] = R[i] - sum;
		}
	}

	//Обратный ход
	void Up(vector<T> Matrix, vector<T> Diagonal, vector<T> R, vector<T> &v)
	{
		v = R;
		for (int i = N - 1; i >= 0; i--)
		{
			v[i] /= Diagonal[i];
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				size_t p = jg[j];
				v[p] -= Matrix[j] * v[i];
			}
		}
	}
	void Up(vector<T> Matrix, vector<T> R, vector<T> &v)
	{
		v = R;
		for (int i = N - 1; i >= 0; i--)
		{
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				size_t p = jg[j];
				v[p] -= Matrix[j] * v[i];
			}
		}
	}

	void Diagonal(vector<T> Diagonalm, vector<T> R, vector<T> &v)
	{
		for (size_t i = 0; i < R.size(); i++)
			v[i] = R[i] / Diagonalm[i];
	}

	//Умноженик числа на вектор
	inline void MultVectorOnT(vector<T> a, T b, vector<T> &c)
	{
		for (size_t i = 0; i < a.size(); i++)
			c[i] = a[i] * b;
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

	void MultMatrixVector(vector<T> f, vector<T> &v)
	{
		v.clear();
		v.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				v[i] += ggl[j] * f[jg[j]];
				v[jg[j]] += ggu[j] * f[i];
			}
			v[i] += di[i] * f[i];
		}
	}
	void MultTMatrixVector(vector<T> f, vector<T> &v)
	{
		v.clear();
		v.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			for (size_t j = ig[i]; j < ig[i + 1]; j++)
			{
				v[i] += ggu[j] * f[jg[j]];
				v[jg[j]] += ggl[j] * f[i];
			}
			v[i] += di[i] * f[i];
		}
	}

	//Residual a - b
	inline void Subtraction(vector<T> a, vector<T> b, vector<T> &v)
	{
		for (size_t i = 0; i < a.size(); i++)
			v[i] = a[i] - b[i];
	}
	//Sum a + b
	inline void Sum(vector<T> a, vector<T> b, vector<T> &v)
	{
		for (size_t i = 0; i < a.size(); i++)
			v[i] = a[i] + b[i];
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
		D.resize(N);
		for (size_t i = 0; i < N; i++)
			D[i] = sqrt(di[i]);
	}
};