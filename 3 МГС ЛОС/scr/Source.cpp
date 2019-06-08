#include "includes.h"
#include "Matrix.h"
#include "chrono"
#include "ctime"

enum {
	Cc = 0,
	DLUSQ = 1,
//	LLT = 2,
	LU = 3
};



#define TIME
#define NOper 100
const int metod = DLUSQ;


int main()
{
#ifdef TIME
	std::chrono::time_point<std::chrono::system_clock> start, end;
	Matrix<double> p;
	p.OpenFile();
	start = chrono::system_clock::now();
	for (size_t i = 0; i < NOper; i++)
	{
		Matrix<double> g = p;
		g.preconditioning(metod);
		//g.ConjugateGradient(metod);
		g.LocalOptimalScheme(metod);
	}
	end = chrono::system_clock::now();
	auto d = end - start;
	auto delt = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "\t Time =  " << delt/NOper << " microseconds\n";
	system("Pause");
#else
	Matrix<double> g;
	g.OpenFile();
	g.preconditioning(metod);
	g.ConjugateGradient(metod);
	//g.LocalOptimalScheme(metod);
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