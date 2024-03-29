// dllmain.cpp : Определяет точку входа для приложения DLL.
#include "stdafx.h"
#include "LDU.cpp"

Matrix::Gause::Gause<double> gause;

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

