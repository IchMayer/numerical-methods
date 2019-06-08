#include "Includes.h"
#include "matrix.h"

template<class T>
matrixLDU<T>::matrixLDU()
{
}

template<class T>
matrixLDU<T>::~matrixLDU()
{
}

template<class T>
bool matrixLDU<T>::A2LDU(matrixA<T> MatrixA, TypeMatrix Type)
{
	switch (Type)
	{
	case TypeMatrix::Full:
	{
#pragma region Не_смотреть
		/*		ЭТОГО НЕТ !!!!!!!!!!!!!!!!!!и не будет!!!!!!!!!!!!!!!!!!!!!!!!!!!
		for (size_t i = 0; i < MatrixA.A.size(); i++)
		{
			for (size_t j = 0; j < i; j++)
			{
				/*T Увеличивает точность->*//* double sum = 0;
				for (size_t k = 0; k < j; k++)
					sum += A[i][k] * A[k][k] * A[k][j];
				A[i][j] = (MatrixA.A[i][j] - sum) / A[j][j];
			}

			/*T Увеличивает точность->*//* double sum = 0;
			for (size_t k = 0; k < i; k++)
				sum += A[i][k] * A[k][k] * A[k][i];
			A[i][i] = MatrixA.A[i][i] - sum;

			for (size_t j = i + 1; j < MatrixA.A.size(); j++)
			{
				/*T Увеличивает точность->*//*/ double sum = 0;
				for (size_t k = 0; k < i; k++)
					sum += A[i][k] * A[k][k] * A[k][j];
				A[i][j] = (MatrixA.A[i][j] - sum) / A[i][i];
			}
		}*/
#pragma endregion
		break;
	}
	case TypeMatrix::Prof:
	{
		ia = MatrixA.ia;
		L.resize(MatrixA.al.size());
		U.resize(MatrixA.au.size());
		D.resize(MatrixA.di.size());
		for (size_t i = 0; i < MatrixA.di.size(); i++)
		{
			for (size_t j = ia[i]; j < ia[i + 1]; j++)
			{
				/*T Увеличивает точность->*/ double sum1 = 0, sum2 = 0;
				for (size_t k = min(j - ia[i], ia[i - ia[i + 1] + j + 1] - ia[i - ia[i + 1] + j]); k > 0; k--)
				{
					sum1 += L[j - k] * D[ia[j + 1] - k] * U[ia[j + 1] - k];
					sum2 += L[ia[j + 1] - k] * D[ia[j + 1] - k] * U[j - k];
				}
				L[j] = (MatrixA.al[j] - sum1) / D[i - ia[i + 1] + j];
				U[j] = (MatrixA.au[j] - sum2) / D[i - ia[i + 1] + j];
			}

			/*T Увеличивает точность->*/ double sum = 0;
			for (size_t k = ia[i]; k < ia[i + 1]; k++)
				sum += L[k] * D[i - ia[i + 1] + k] * U[k];
			D[i] = MatrixA.di[i] - sum;
			if(!D[i])
			{
				int a;
				cout << "Error: matrix don't have single solution" << endl;
				cin >> a;
				return false;
			}
		}
		break;
	}
	default:
		break;
	}
	return true;
}

template<class T>
bool matrixLDU<T>::LFX()
{
	ifstream iF("F.txt");
	if (!iF)
	{
		cout << "Error File F" << endl;
		return false;
	}
	F.resize(D.size());
	Y.resize(F.size());
	for (size_t i = 0; i < F.size(); i++)
		iF >> F[i];
	iF.close();

	for (size_t i = 0; i < Y.size(); i++)
	{
		for (size_t m = ia[i]; m < ia[i + 1]; m++)
		{
			size_t j = i - ia[i + 1] + m;
			Y[i] -= Y[j] * L[m];
		}
		Y[i] += F[i];
	}

	Z.resize(Y.size());
	for (size_t i = 0; i < Z.size(); i++)
		Z[i] = Y[i] / D[i];

	X.resize(Z.size());
	for (int i = X.size() - 1; i >= 0; i--)
	{
		for (size_t j = i + 1; j < X.size(); j++)
		{
			if (ia[j + 1] - ia[j] >= j - i)
			{
				double l1 = X[j];
				double l2 = U[ia[j + 1] - j + i];
				X[i] -= X[j] * U[ia[j + 1] - j + i];
			}
		}
		X[i] += Z[i];
	}

	return true;
}

template<class T>
bool matrixLDU<T>::SaveFile()
{
	ofstream oX("X.txt");
	if (!oX)
	{
		cout << "Error: File X";
		return false;
	}
	for (size_t i = 0; i < X.size(); i++)
		oX << X[i] << endl;
	return true;
}

template<class T>
bool matrixLDU<T>::SaveLDU(TypeMatrix Type)
{
	switch (Type)
	{
	case TypeMatrix::Full:
	{
	ifstream ifs("input.txt");
		if (!ifs) {
			cout << "File error." << endl;
			return 0;
		}
		size_t n;
		ifs >> n;
		A.reserve(n);
		for (size_t i = 0; i < n; i++)
		{
			A[i].reserve(n);
		for (size_t j = 0; j < n; j++)
				ifs >> A[i][j];
		}
		ifs.close();
		break;
	}
	case TypeMatrix::Prof:
	{

		ofstream oL("L.txt");
		for (size_t i = 0; i < L.size(); i++)
			oL << L[i] << endl;
		oL.close();

		ofstream oU("U.txt");
		for (size_t i = 0; i < U.size(); i++)
			oU << U[i] << endl;
		oU.close();

		ofstream oD("D.txt");
		for (size_t i = 0; i < D.size(); i++)
			oD << D[i] << endl;
		oD.close();

		break;
	}
	default:
		return 0;
	}
	return 1;
}

template<class T>
matrixA<T>::matrixA()
{
	A.reserve(0);
	al.reserve(0);
	au.reserve(0);
	ia.reserve(0);
	di.reserve(0);
}

template<class T>
inline matrixA<T>::matrixA(size_t n)
{
	A.reserve(n);
	for (size_t i = 0; i < n; A[i++].reserve(n));
}

template<class T>
matrixA<T>::~matrixA()
{
}


/*
	HighElemGausStraight - Функция поиска решения СЛАУ методом Гаусса с выбором главного коэфициента.
	Матрица будет приводиться к верхнему треугольному виду
*/
template<class T>
bool matrixA<T>::HighElemGausStraight()
{
	ifstream iF("F.txt");
	if (!iF)
	{
		cout << "Error File F" << endl;
		return false;
	}
	F.resize(A.size());
	X.resize(F.size());
	for (size_t i = 0; i < F.size(); i++)
		iF >> F[i];
	iF.close();

	for (size_t i = 0; i < A.size(); i++)
	{
		//Опредиление главного Элемента в столбце
		T Max = A[i][i];
		size_t index = i;
		for (size_t j = i+1; j < A.size(); j++)
		{
			if (abs(A[j][i]) > abs(Max))
			{
				Max = A[j][i];
				index = j;
			}
		}

		//Проверка единственность
		if (!Max)
		{
			int a;
			cout << "Error: matrix don't have single solution" << endl;
			cin >> a;
			return false;
		}

		//Свап строки
		if (i != index)
		{
			A[index].swap(A[i]);
			T buf = F[i];
			F[i] = F[index];
			F[index] = buf;
		}

		//Постройка верхней треугольной матрицы
		for (size_t j = i + 1; j < A.size(); j++)
		{
			if (A[j][i] == 0)
				continue;
			double c = A[j][i] / Max;
//			A[j][i] = 0;
			for (size_t l = i; l < A.size(); l++)
				A[j][l] -= A[i][l] * c;
			double l1 = F[i] * c;
			F[j] -= F[i] * c;
		}
	}
	return true;
}

template<class T>
bool matrixA<T>::HighElemGausBack()
{
	for (int i = A.size() - 1; i >= 0; i--)
	{
		double sum = 0;
		for (size_t j = i + 1; j < A.size(); j++)
			sum += A[i][j] * X[j];
		X[i] = T(double((F[i] - sum)/A[i][i]));
	}
	return true;
}

template<class T>
bool matrixA<T>::Prof2Full()
{
	if (!A.size())
		if (!di.size())
			return false;
		else
		{
			A.resize(di.size());
			for (size_t i = 0; i < di.size(); A[i++].resize(di.size()));
		}
		
	for (size_t i = 0; i < di.size(); i++)
		A[i][i] = di[i];
		
	for (size_t i = 0; i < di.size(); i++)
	{
		for (size_t m = ia[i]; m < ia[i + 1]; m++)
		{
			size_t j = i - ia[i + 1] + m;
			A[i][j] = al[m];
			A[j][i] = au[m];
		}
	}
	return true;
}

template<class T>
void matrixA<T>::PushBack(vector<T> Str)
{
	A.push_back(Str);
}

template<class T>
bool matrixA<T>::OpenFile(TypeMatrix typefile)
{
	switch (typefile)
	{
	case TypeMatrix::Full:
	{
		ifstream ifs("input.txt");
		if (!ifs) {
			cout << "File error." << endl;
			return 0;
		}
		size_t n;
		ifs >> n;

		A.reserve(n);

		for (size_t i = 0; i < n; i++)
		{
			A[i].reserve(n);
			for (size_t j = 0; j < n; j++)
				ifs >> A[i][j];
		}

		ifs.close();
		break;
	}
	case TypeMatrix::Prof:
	{
		size_t n, m;

		ifstream info("info.txt");
		info >> n >> m;
		info.close();

		ifstream ial("al.txt");
		al.resize(m);
		for (size_t i = 0; i < m; i++)
			ial >> al[i];
		ial.close();

		ifstream iau("au.txt");
		au.resize(m);
		for (size_t i = 0; i < m; i++)
			iau >> au[i];
		iau.close();

		ifstream idi("di.txt");
		di.resize(n);
		for (size_t i = 0; i < n; i++)
			idi >> di[i];
		idi.close();

		ifstream iia("ia.txt");
		ia.resize(n + 1);
		for (size_t i = 0; i < n; i++)
			iia >> ia[i];
		iia >> ia[n];
		iia.close();
		break;
	}
	default:
		return 0;
	}
	return 1;
}

template<class T>
bool matrixA<T>::SaveFile()
{
	ofstream oX("X.txt");
	if (!oX)
	{
		cout << "Error: File X";
		return false;
	}
	for (size_t i = 0; i < X.size(); i++)
		oX << X[i] << endl;
	return true;
}