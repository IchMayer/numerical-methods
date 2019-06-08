#pragma once
#include "stdafx.h"

//Метод гаусса для полной матрицы
typedef void(*Gause)(std::vector<std::vector<double>> &A, std::vector<double> &F, std::vector<double> &X);

using namespace std;

//#define Anali
#define Test5

#pragma region Tests

#ifdef Test1
//n = 2 m = 3;

inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0], 2) + pow(x[1], 2) - 4; };
	//Func[1] = [](const vector<double> &x) {return pow(x[0] - 2, 2) / -2 + x[1]; };
	//Func[2] = [](const vector<double> &x) {return x[1] - 2; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0]; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1]; };
	//AnalJ[1][0] = [](const vector<double> &x) {return 2 - x[0]; };
	//AnalJ[1][1] = [](const vector<double> &x) {return 1; };
	//AnalJ[2][0] = [](const vector<double> &x) {return 0; };
	//AnalJ[2][1] = [](const vector<double> &x) {return 1; };
}
#endif // Test1

#ifdef Test2
//n = 3 m = 2
inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0], 2) + pow(x[1], 2) + pow(x[2],2) - 4; };
	Func[1] = [](const vector<double> &x) {return x[0]-x[1] + x[2]; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0]; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1]; };
	AnalJ[0][2] = [](const vector<double> &x) {return 2 * x[2];	};
	AnalJ[1][0] = [](const vector<double> &x) {return 1; };
	AnalJ[1][1] = [](const vector<double> &x) {return -1; };
	AnalJ[1][2] = [](const vector<double> &x) {return 1; };
}

#endif 

#ifdef Test3
//n = 2 m = 2
inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0] - 2, 2) + pow(x[1] - 1, 2) - 0.5; };
	Func[1] = [](const vector<double> &x) {return pow(x[0] - 4, 2) + pow(x[1] - 2, 2) - 2; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0]; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1]; };
	AnalJ[1][0] = [](const vector<double> &x) {return 2 * x[0] + 6; };
	AnalJ[1][1] = [](const vector<double> &x) {return 2 * x[1]; };
}

#endif 

#ifdef Test4
//n = 2 m = 2
inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0] - 2, 2) + pow(x[1] - 1, 2) - 0.5; };
	Func[1] = [](const vector<double> &x) {return pow(x[0] - 4, 2) + pow(x[1] - 2, 2) - 2; };
	Func[2] = [](const vector<double> &x) {return 0.5*(pow(x[0] - 3, 2) + pow(x[1] + 1, 2) - 1); };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0] - 4; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1] - 2; };
	AnalJ[1][0] = [](const vector<double> &x) {return 2 * x[0] - 8; };
	AnalJ[1][1] = [](const vector<double> &x) {return 2 * x[1] - 4; };
	AnalJ[2][0] = [](const vector<double> &x) {return x[0] - 3; };
	AnalJ[2][1] = [](const vector<double> &x) {return x[1] + 1; };
}

#endif 

#ifdef Test5
//n = 2 m = 2
inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0], 2) + pow(x[1], 2) - 4; };
	Func[1] = [](const vector<double> &x) {return pow(x[0] + 1, 2) + pow(x[1], 2) - 1; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0]; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1]; };
	AnalJ[1][0] = [](const vector<double> &x) {return 2 * x[0] + 2; };
	AnalJ[1][1] = [](const vector<double> &x) {return 2 * x[1]; };
}

#endif 

#ifdef Test6
//n = 2 m = 2
inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0], 2) + pow(x[1], 2) - 4; };
	Func[1] = [](const vector<double> &x) {return pow(x[0], 2) + pow(x[1], 2) - 1; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0]; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1]; };
	AnalJ[1][0] = [](const vector<double> &x) {return 2 * x[0]; };
	AnalJ[1][1] = [](const vector<double> &x) {return 2 * x[1]; };
}

#endif 

#ifdef Test7
//n = 2 m = 3
inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return (2* x[0] + x[1] + 2)/100; };
	Func[1] = [](const vector<double> &x) {return x[0]; };
	Func[2] = [](const vector<double> &x) {return x[1]; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 1; };
	AnalJ[0][1] = [](const vector<double> &x) {return 1; };
	AnalJ[1][0] = [](const vector<double> &x) {return 0; };
	AnalJ[1][1] = [](const vector<double> &x) {return 1; };
	AnalJ[2][0] = [](const vector<double> &x) {return 1; };
	AnalJ[2][1] = [](const vector<double> &x) {return 0; };
}

#endif 

#ifdef Test8

inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0] - 2, 2) + pow(x[1] - 2, 2) - 4; };
	Func[1] = [](const vector<double> &x) {return pow(x[0] - 5, 2) + pow(x[1] - 2, 2) - 1; };
	Func[2] = [](const vector<double> &x) {return x[1] - 2; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0] - 4; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1] - 4; };
	AnalJ[1][0] = [](const vector<double> &x) {return 2 * x[0] - 10; };
	AnalJ[1][1] = [](const vector<double> &x) {return 2 * x[1] - 4; };
	AnalJ[2][0] = [](const vector<double> &x) {return 1; };
	AnalJ[2][1] = [](const vector<double> &x) {return 0; };
}
#endif // Test8

#ifdef Test9

inline void InItFunction(function<double(vector<double>)> Func[])
{
	Func[0] = [](const vector<double> &x) {return pow(x[0] - 2, 2) + pow(x[1] - 1, 2) - 4; };
	Func[1] = [](const vector<double> &x) {return pow(x[0] - 7, 2) + pow(x[1] - 1, 2) - 4; };
	Func[2] = [](const vector<double> &x) {return x[1] - 2; };
}

inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
{
	AnalJ[0][0] = [](const vector<double> &x) {return 2 * x[0] - 4; };
	AnalJ[0][1] = [](const vector<double> &x) {return 2 * x[1] - 2; };
	AnalJ[1][0] = [](const vector<double> &x) {return 2 * x[0] - 14; };
	AnalJ[1][1] = [](const vector<double> &x) {return 2 * x[1] - 2; };
	AnalJ[2][0] = [](const vector<double> &x) {return 1; };
	AnalJ[2][1] = [](const vector<double> &x) {return 0; };
}
#endif // Test8

//inline void InItFunction(function<double(vector<double>)> Func[])
//{
//	Func[0] = [](const vector<double> &x) {return x[0] + x[1]; };
//	Func[1] = [](const vector<double> &x) {return x[0] - 2; };
//}
//
//inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
//{
//	AnalJ[0][0] = [](const vector<double> &x) {return 1; };
//	AnalJ[0][1] = [](const vector<double> &x) {return 1; };
//	AnalJ[1][0] = [](const vector<double> &x) {return 1; };
//	AnalJ[1][1] = [](const vector<double> &x) {return 0; };
//}

//inline void InItFunction(function<double(vector<double>)> Func[])
//{
//	Func[0] = [](const vector<double> &x) {return x[0]; };
//	Func[1] = [](const vector<double> &x) {return x[0]- 2; };
//	Func[2] = [](const vector<double> &x) {return x[1] - 2; };
//}

//inline void InItJacobi(function<double(vector<double>)> AnalJ[3][2])
//{
//	AnalJ[0][0] = [](const vector<double> &x) {return 1; };
//	AnalJ[0][1] = [](const vector<double> &x) {return 0; };
//	AnalJ[1][0] = [](const vector<double> &x) {return 2; };
//	AnalJ[1][1] = [](const vector<double> &x) {return 0; };
//	AnalJ[2][0] = [](const vector<double> &x) {return 0; };
//	AnalJ[2][1] = [](const vector<double> &x) {return 1; };
//}

#pragma endregion

class chm4
{
public:
	chm4()
	{
		InItFunction(Func);
		InItJacobi(AnalJ); 
		h = k = 0;
		Move.resize(3);
		Move[0] = D2D1::Point2F(500, 375);
		Move[1] = D2D1::Point2F(-10, 10);
		Move[2] = D2D1::Point2F(7.5, -7.5);
	}
	~chm4() {}
	int maxiter, m, n;
	vector<double> x, dx;
	double e1, e2, h;
	double NormF0, NormF;
	double FMax;
	vector<vector<double>> Jacobi;
	vector<double> F;
	vector<D2D1_POINT_2F> Target;
	vector<D2D1_POINT_2F> Move;
	double b;
	int k;
	void TwoCicle()
	{
		auto hinstLib = LoadLibrary(TEXT("LDU.dll"));
		Gause gause = (Gause)GetProcAddress(hinstLib, "Gause");
		k = 0;
		int p = min(m, n);
		F.resize(p);
		dx.resize(n);
		Jacobi.resize(p);
		for (size_t i = 0; i < p; i++)
		{
			Jacobi[i].resize(p);
			F[i] = -Func[i](x);
		}
		NormF = NormF0 = Norm(F);
		b = 1;
		D2D1_POINT_2F o;
		if (n == 2)
		{
			o = D2D1::Point2F(x[0], x[1]);
			Target.push_back(o);
		}
		log.open("log.txt", ios::app);
		log << "------------------------------------------" << endl;
		PrintLog();
		while (CheckEnd())
		{
			Fbrake();
			NormF = Norm(F);
			gause(Jacobi, F, dx);
			if (Norm(dx) == 0)
			{
				k++; break;
			}
			FindB();
			Sum(x, dx, b);
			if (n == 2)
			{
				o = D2D1::Point2F(x[0], x[1]);
				Target.push_back(o);
			}
			k++;
			PrintLog();
		}
		log.close();
		Shadow();
	}
	void Init()
	{
		ifstream in("info.txt");
		in >> e1 >> e2 >> maxiter;
		in >> n >> m;
		in.close();
		in.open("v0.txt");
		x.resize(n);
		for (size_t i = 0; i < n; i++)
			in >> x[i];
		FMax = GetNorm(-13.0, -13.0);
		if (!FMax || FMax == -1)
			FMax = GetNorm(-13.0, 13.0);
	}
	D2D1_POINT_2F Result()
	{
		if (n != 2)
			return D2D1::Point2F(0.0f, 0.0f);
		else
			return D2D1::Point2F(x[0], x[1]);
	}
	double GetNorm(double X, double Y)
	{
		double sum = 0;
		vector<double> r(2);
		r[0] = X;
		r[1] = Y;
		for (size_t i = 0; i < m; i++)
		{
			double l = pow(Func[i](r), 2);
			if (l < 1e-3)
				return -1;
			else
				sum += l;
		}
		return  sqrt(sum);
	}
private:
	//
	function<double(vector<double>)> Func[3];
	function<double(vector<double>)> AnalJ[3][2];
	vector<int> index;
	ofstream log;
	struct br
	{
		double f;
		int i;
	};
	double Norm(vector<double> x)
	{
		double sum = 0;
		for (size_t i = 0; i < x.size(); i++)
			sum += pow(x[i], 2);
		return  sqrt(sum);
	}
	bool CheckEnd()
	{
		return k < maxiter && NormF / NormF0 >= e1 && b >= e2;
	}
	void FindB()
	{
		b = 2;
		double norm = NormF;
		vector<double> F(m);
		vector<double> a;
		do
		{
			b /= 2;
			for (size_t i = 0; i < m; i++)
			{
				a = x;
				Sum(a, dx, b);
				F[i] = Func[i](a);
			}
			norm = Norm(F);
		} while (norm > NormF && b > e2);
	}
	void Sum(vector<double> &a, vector<double> &b, double c)
	{
		if (b.size() < a.size())
		{
			for (size_t i = 0; i < index.size(); i++)
				a[index[i]] += c * b[i];
			return;
		}
		for (size_t i = 0; i < a.size(); i++)
			a[i] += c * b[i];
	}
	double ChisJ(int i, int j, vector<double> x)
	{
		double oo = 1e-2;
		x[j] -= oo/2;
		double X1 = Func[i](x);
		x[j] += oo;
		double X2 = Func[i](x);
		return (X2 - X1) / oo;
	}
	void Fbrake()
	{
		index.resize(0);
		if (n == m)
		{
			for (size_t i = 0; i < m; i++)
			{
				for (size_t j = 0; j < n; j++)
#ifdef  Anali
					Jacobi[i][j] = AnalJ[i][j](x);
#else
					Jacobi[i][j] = ChisJ(i, j, x);
#endif //  Anali
				F[i] = -Func[i](x);
			}
		}
		else
		{
			if (m > n)
			{
				vector<br> brf(m - n);
				for (size_t i = 0; i < m - n; i++)
					brf[i].f = Func[i](x);
				sort(brf.begin(), brf.end(), [](br i, br j) -> bool { return (i.f < j.f); });
				for (size_t i = m - n, i2 = i - 1, i3 = 0; i < m; i++, i3++)
				{
					if (abs(Func[i](x)) < abs(brf[i2].f))
					{
						F[i3] = -brf[i2].f;
						for (size_t j = 0; j < n; j++)
#ifdef Anali
							Jacobi[i3][j] = AnalJ[brf[i2].i][j](x);
#else
							Jacobi[i3][j] = ChisJ(brf[i2].i, j, x);
#endif
						brf[i2].f = Func[i](x);
						brf[i2].i = i;
						sort(brf.begin(), brf.end(), [](br i, br j) -> bool { return (i.f < j.f); });
					}
					else
					{
						F[i3] = -Func[i](x);
						for (size_t j = 0; j < n; j++)
#ifdef  Anali
							Jacobi[i3][j] = AnalJ[i][j](x);
#else
							Jacobi[i3][j] = ChisJ(i, j, x);
#endif 
					}
				}
			}
			else
			{
				vector<br> brf(n - m);
				for (size_t i = 0; i < n - m; i++)
					brf[i].f = MaxJacobi(i);
				sort(brf.begin(), brf.end(), [](br i, br j) -> bool { return (i.f < j.f); });
				for (size_t i = n - m, i2 = i - 1, i3 = 0; i < n; i++, i3++)
				{
					double p = MaxJacobi(i);
					double ggr = brf[i2].f;
					if (p <= ggr)
					{
						F[i3] = -Func[i3](x);
						index.push_back(brf[i2].i);
						for (size_t j = 0; j < m; j++)
#ifdef Anali
							Jacobi[j][i3] = AnalJ[j][brf[i2].i](x);
#else
							Jacobi[j][i3] = ChisJ(j,brf[i2].i, x);
#endif
						brf[i2].f = p;
						brf[i2].i = i;
						sort(brf.begin(), brf.end(), [](br i, br j) -> bool { return (i.f < j.f); });
					}
					else
					{
						F[i3] = -Func[i3](x);
						index.push_back(i);
						for (size_t j = 0; j < m; j++)
#ifdef Anali
							Jacobi[j][i3] = AnalJ[j][i](x);
#else
							Jacobi[j][i3] = ChisJ(j, i, x);
#endif // Anali
					}
				}
			}
		}
	}
	double MaxJacobi(int j)
	{
		if (!n)
			return 0;
#ifdef Anali
		double max = abs(AnalJ[0][j](x));
#else
		double max = abs(ChisJ(0, j, x));
#endif // Anali
		for (size_t i = 1; i < m; i++)
		{
#ifdef Anali
			double p = abs(AnalJ[i][j](x));
#else
			double p = abs(ChisJ(i, j, x));
#endif // Anali
			if (max < p)
				max = p;
		}
		return max;
	}
	void Shadow()
	{
		if (n != 2)
			return;
		else
			Sum(dx, x, 1);
	}
	void PrintLog()
	{

		log << "k: ";
	//	auto qq = log.precision(15);
		log << k << " " << fixed;
		for (size_t i = 0; i < x.size(); i++)
		{
			log << x[i]<<" ";
		}
		log<< scientific << "b: " << b << " F: "<<NormF << endl;
	//	log.precision(qq);
	}
};