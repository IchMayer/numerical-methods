// LDU.cpp : Определяет экспортированные функции для приложения DLL.
//

#include "stdafx.h"


namespace Matrix
{
	enum class Errors
	{
		noError,
		ErrorOpenFile,
		MatrixNotCount
	};

	enum class Format
	{
		LDU,
		Profile
	};

	template<class T>
	class Matrix
	{
	public:
		Matrix() {}
		~Matrix() {}

		std::vector<T> al, au, di, f, x;
		std::vector<size_t> ia;

		size_t n, m;
		Format format;
		Errors error;
	private:

	};

	namespace Gause
	{
		template<class T>
		class Gause : Matrix<T>
		{
		public: 
			std::vector<std::vector<T>> A;
			std::vector<double> f, x;

			Gause()
			{
				this->format = Format::Profile;
				error = Errors::noError;
			}
			~Gause() {}
			Errors error;

		private:
			size_t n; //Размерность матрицы
		};

		template<class T>
		void GausBack(Gause<T> &matrix)
		{
			std::vector<T> &X = matrix.x;
			std::vector<T> &F = matrix.f;
			std::vector<std::vector<T>> &A = matrix.A;
			auto n = F.size();
			X.resize(n);
			if (matrix.error != Errors::noError)
				return;
			for (int i = n - 1; i >= 0; i--)
			{
				double sum = 0;
				for (size_t j = i + 1; j < n; j++)
					sum += A[i][j] * X[j];
				X[i] = (F[i] - sum) / A[i][i];
			}
		}

		template<class T>
		void CountUpTrianMatrix(Gause<T> &matrix)
		{
			if (matrix.error != Errors::noError)
				return;
			std::vector<std::vector<T>> &A = matrix.A;
			std::vector<T> &F = matrix.f;
			auto n = F.size();
			for (size_t i = 0; i < n; i++)
			{
				//Опредиление главного Элемента в столбце
				T Max = A[i][i];
				size_t index = i;
				for (size_t j = i + 1; j < n; j++)
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
					matrix.error = Errors::MatrixNotCount;
					return;
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
				for (size_t j = i + 1; j < n; j++)
				{
					if (!A[j][i])
						continue;
					T c = A[j][i] / Max;
					A[j][i] = 0.;
					for (size_t l = i + 1; l < n; l++)
						A[j][l] -= A[i][l] * c;
					F[j] -= F[i] * c;
				}
			}
		}
		
		template<class T>
		Gause<T> Create(std::vector<double> &al, std::vector<double> &au, std::vector<double> &di, std::vector<size_t> &ia)
		{
			Gause<T> m;
			std::vector<std::vector<T>> &A = m.A;
			if (!di.size())
			{
			 	m.Error = Errors::MatrixNotCount;
				return;
			}
			else
			{
				A.resize(m.n);
				for (size_t i = 0; i < m.n; A[i++].resize(m.n));
			}
			
			m.n = di.size();
			auto n = m.n;
			for (size_t i = 0; i < n; i++)
				A[i][i] = di[i];

			for (size_t i = 0; i < n; i++)
			{
				for (size_t m = ia[i], j = i - ia[i + 1] + m; m < ia[i + 1]; m++, j++)
				{
					A[i][j] = al[m];
					A[j][i] = au[m];
				}
			}
		}

	}

	namespace LDU
	{
		template<class T>
		class LDU: Matrix<T>
		{
		public:
			LDU<T>()
			{
				this->format = Format::LDU;
			}
			~LDU<T>() {}
			std::vector<T> al, au, di, f, x;
			std::vector<size_t> ia;

			size_t n, m;
			Format format;
			Errors error;
		};

		template<class T>
		void SaveFile(LDU<T> &matrix)
		{

		}

		template<class T>
		void CountX(LDU<T> &matrix, std::vector<T> F)
		{
			size_t n = matrix.n;
			std::vector<T> &D = matrix.di;
			std::vector<T> &U = matrix.au;
			std::vector<T> &L = matrix.al;
			std::vector<T> &X = matrix.x;
			Errors Error = matrix.error;
			std::vector<size_t> &ia = matrix.ia;

			if (Error != Errors::noError)
				return;
			//Поиск Y
			std::vector<T> &Y = X;
			for (size_t i = 0; i < n; i++)
			{
				double sum = 0;
				for (size_t m = ia[i],
					j = i - ia[i + 1] + m;
					m < ia[i + 1]; m++)
					sum += Y[j++] * L[m];
				Y[i] += F[i] - sum;
			}

			//Поиск Z
			std::vector<T> &Z = X;
			for (size_t i = 0; i < n; i++)
				Z[i] /= D[i];

			//Поиск X
			for (int i = n - 1; i >= 0; i--)
				for (int j = ia[i + 1] - ia[i] - 1,
					j1 = i - j - 1,
					j2 = ia[i + 1] - j - 1;
					j >= 0; j--)
					X[j1++] -= X[i] * U[j2++];
		}

		template<class T>
		void CountLDU(LDU<T> &matrix)
		{
			size_t n = matrix.n;
			size_t m = matrix.m;
			auto &D = matrix.di;
			auto &U = matrix.au;
			auto &L = matrix.al;
			auto &Error = matrix.error;
			auto &ia = matrix.ia;

			double CompareNumber;
			if (std::is_same<T, double>::value)
				CompareNumber = pow(10, 15);
			else
				CompareNumber = pow(10, 7);

			for (size_t i = 0; i < n; i++)
			{
				double sum = 0;
				for (size_t j = ia[i], j1 = i - ia[i + 1] + j; j < ia[i + 1]; j++)
				{
					double sum1 = 0, sum2 = 0;
					for (size_t k = min(j - ia[i], ia[i - ia[i + 1] + j + 1] - ia[i - ia[i + 1] + j]),
						i1 = j - k,
						i2 = i - ia[i + 1] + j - k,
						i3 = ia[i - ia[i + 1] + j + 1] - k;
						k > 0; k--)
					{
						sum1 += L[i1] * D[i2] * U[i3];
						sum2 += L[i3++] * D[i2++] * U[i1++];
					}
					L[j] = (L[j] - sum1) / D[j1];
					U[j] = (U[j] - sum2) / D[j1];
					sum += L[j] * D[j1++] * U[j];
				}
				//Если сделать проверка на inf после деления на 0, то можем получить не то что нам нужно
				//Т.к. при разность (D[i] - sum) может получиться не 0, а число близкое к нему(разность нецелого типа)
				//и в таком случае при делении числа на эту разность может получиться не inf.
				if (abs((D[i] - sum)) < abs(D[i] / CompareNumber))
				{
					Error = Errors::MatrixNotCount;
					return;
				}
				D[i] -= sum;
			}
		}

		template <class T>
		void CountLDU(LDU<T> &matrix, std::vector<T> &al, std::vector<T> &au, std::vector<T> &di, std::vector<size_t> &ia)
		{
			size_t n = matrix.n;
			size_t m = matrix.m;
			auto &D = matrix.di;
			auto &U = matrix.au;
			auto &L = matrix.al;
			auto &Error = matrix.error;
			ia._Copy_alloc(matrix.ia.get_allocator());

			double CompareNumber;
			if (std::is_same<T, double>::value)
				CompareNumber = pow(10, 15);
			else
				CompareNumber = pow(10, 7);

			for (size_t i = 0; i < n; i++)
			{
				double sum = 0;
				for (size_t j = ia[i], j1 = i - ia[i + 1] + j; j < ia[i + 1]; j++)
				{
					double sum1 = 0, sum2 = 0;
					for (size_t k = min(j - ia[i], ia[i - ia[i + 1] + j + 1] - ia[i - ia[i + 1] + j]),
						i1 = j - k,
						i2 = i - ia[i + 1] + j - k,
						i3 = ia[i - ia[i + 1] + j + 1] - k;
						k > 0; k--)
					{
						sum1 += L[i1] * D[i2] * U[i3];
						sum2 += L[i3++] * D[i2++] * U[i1++];
					}
					L[j] = (al[j] - sum1) / D[j1];
					U[j] = (au[j] - sum2) / D[j1];
					sum += L[j] * D[j1++] * U[j];
				}
				//Если сделать проверка на inf после деления на 0, то можем получить не то что нам нужно
				//Т.к. при разность (D[i] - sum) может получиться не 0, а число близкое к нему(разность нецелого типа)
				//и в таком случае при делении числа на эту разность может получиться не inf.
				if (abs((D[i] - sum)) < abs(D[i] / CompareNumber))
				{
					Error = Errors::MatrixNotCount;
					return;
				}
				D[i] = di[i] - sum;
			}
		}

		template<class U>
		bool Open(std::string path, std::vector<U> &a, size_t len)
		{
			std::ifstream ia;
			ia.exceptions(std::ifstream::failbit | std::ifstream::badbit);
			try
			{
				ia.open(path);
				a.resize(len);
				for (size_t i = 0; i < len; i++)
					ia >> a[i];
				ia.close();
			}
			catch (std::ifstream::failure e)
			{
				std::cout << "Error Open File: " << path << std::endl;
				return false;
			}
			return true;
		}

		template<class T>
		void OpenFile(LDU<T> &matrix, std::string path)
		{
			int error = 0;
			size_t n, m;
			if (path.size())
				path += "/";
			std::vector<size_t> info(2);
			error += 
				!Open(path + "info.txt", info, 2);
			n = info[0], m = info[1];
			matrix.x.resize(n);
			error +=
				!Open(path + "al.txt", matrix.al, m) +
				!Open(path + "au.txt", matrix.au, m) +
				!Open(path + "di.txt", matrix.di, n) +
				!Open(path + "ia.txt", matrix.ia, n + 1) +
				!Open(path + "F.txt", matrix.f, n);
			if(error)
				matrix.error = Errors::ErrorOpenFile;
		}

		template<class T>
		LDU<T> Create(std::string path)
		{
			LDU<T> m;
			OpenFile(m, path);
			CountLDU(m);
			return m;
		}

		template<class T>
		LDU<T> Create(std::vector<double> &al, std::vector<double> &au, std::vector<double> &di, std::vector<size_t> &ia)
		{
			LDU<T> m;
			CountLDU(m, al, au, di, ia);
			return m;
		}
	}
}


