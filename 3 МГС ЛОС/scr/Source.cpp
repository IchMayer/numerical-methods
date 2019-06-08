#include "includes.h"
#include "Matrix.h"
#include "Windows.h"

enum {
	Cc = 0,
	DLUSQ = 1,
	//	LLT = 2,
	LU = 3
};



#define TIME
const int metod = DLUSQ;


#ifdef TIME
SYSTEMTIME operator-(const SYSTEMTIME& pSr, const SYSTEMTIME& pSl)
{
	SYSTEMTIME t_res;
	FILETIME v_ftime;
	ULARGE_INTEGER v_ui;
	__int64 v_right, v_left, v_res;
	SystemTimeToFileTime(&pSr, &v_ftime);
	v_ui.LowPart = v_ftime.dwLowDateTime;
	v_ui.HighPart = v_ftime.dwHighDateTime;
	v_right = v_ui.QuadPart;

	SystemTimeToFileTime(&pSl, &v_ftime);
	v_ui.LowPart = v_ftime.dwLowDateTime;
	v_ui.HighPart = v_ftime.dwHighDateTime;
	v_left = v_ui.QuadPart;

	v_res = v_right - v_left;

	v_ui.QuadPart = v_res;
	v_ftime.dwLowDateTime = v_ui.LowPart;
	v_ftime.dwHighDateTime = v_ui.HighPart;
	FileTimeToSystemTime(&v_ftime, &t_res);
	return t_res;
}
#endif // TIME


int main()
{
#ifdef TIME
	SYSTEMTIME start, end;
	GetSystemTime(&start);
	for (size_t i = 0; i < 1; i++)
	{
		Matrix<double> g;
		g.OpenFile();
		g.preconditioning(metod);
		g.ConjugateGradient(metod);
		//g.LocalOptimalScheme(metod);
	}
	GetSystemTime(&end);

	cout << "\t T =  " << (end - start).wSecond << "." << (end - start).wMilliseconds;
	int d;
	cin >> d;
#else
	Matrix<double> g;
	g.OpenFile();
	g.preconditioning(metod);
	//g.ConjugateGradient(metod);
	g.LocalOptimalScheme(metod);
	vector<double> out = g.Result();
	int d;
	cout << g.K() << endl;
	cout.setf(ios::scientific);
	cout.precision(15);
	for (size_t i = 0; i < out.size(); i++)
		cout << out[i] << "|" << i - out[i] + 1 << "\n";
	cout.unsetf(ios::scientific);
	cin >> d;
#endif // TIME
	return 0;
}