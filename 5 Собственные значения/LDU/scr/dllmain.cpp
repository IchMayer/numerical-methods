// dllmain.cpp : Определяет точку входа для приложения DLL.
#include "stdafx.h"
#include "LDU.cpp"

Matrix::Gause::Gause<double> gause;
Matrix::LDU::LDU<double> ldu;
Matrix::LU::LU<double> lu;
Matrix::LUFull::LUFull<double> luf;

extern "C"
{
	_declspec(dllexport) void Gause(std::vector<std::vector<double>> &A, std::vector<double> &F, std::vector<double> &X)
	{
		gause.A = A;
		gause.f = F;
		Matrix::Gause::CountUpTrianMatrix(gause);
		Matrix::Gause::GausBack(gause);
		X = gause.x;
	}

	_declspec(dllexport) void LDUInit(std::vector<std::vector<double>> &A)
	{
		Matrix::LDU::descript(ldu, A);
		Matrix::LDU::CountLDU(ldu);
	}

	_declspec(dllexport) void LDU(std::vector<double> &F, std::vector<double> &X)
	{
		Matrix::LDU::CountX(ldu, F);
		X = ldu.x;
	}

	_declspec(dllexport) void LUInit(std::vector<std::vector<double>> &A)
	{
		Matrix::LU::descript(lu, A);
		Matrix::LU::CountLU(lu);

		}

	_declspec(dllexport) void LU(std::vector<double> &F, std::vector<double> &X)
	{
		Matrix::LU::CountX(lu, F);
		X = lu.x;
	}

	_declspec(dllexport) void LUFInit(std::vector<std::vector<double>> &A)
	{
		Matrix::LUFull::descriptLU(luf, A);
	}

	_declspec(dllexport) void LUF(std::vector<double> &F, std::vector<double> &X)
	{
		Matrix::LUFull::CountX(luf, F);
		X = luf.x;
	}
}

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
	setlocale(LC_ALL, "Russian");
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

