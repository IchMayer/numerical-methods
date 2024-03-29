// SLAU.cpp: Определяет точку входа для приложения.
//

#include "stdafx.h"
#include "SLAU.h"
#include "matrix.cpp"

#pragma region Gaos

template <class T>
void Gaus(HWND hDlg)
{
	matrixA<T> A;
	wchar_t path[255];
	GetDlgItemText(hDlg, IDC_EDIT2, path, sizeof(path));
	_bstr_t p(path);
	A.OpenFiles(TypeMatrix::Prof, p);
	A.Prof2Full();
	A.CountUpTrianMatrix();
	A.GausBack();
	if (A.Error)
		Error2Edit(A.Error, hDlg, IDC_EDIT3);
	else
	{
		vector2edit(A.X, hDlg, IDC_EDIT3);
		X<T>.X = A.X;
	}
}

template <class T>
void Gaus(HWND hDlg, int k)
{
	matrixA<T> A;
	wchar_t path[255];
	GetDlgItemText(hDlg, IDC_EDIT2, path, sizeof(path));
	_bstr_t p(path);
	if (k == -1) {
		A.OpenFiles(TypeMatrix::Prof, p);
	}else
		A.OpenFiles(TypeMatrix::Prof, p, k);
	A.Prof2Full();
	A.CountUpTrianMatrix();
	A.GausBack();
	if(A.Error)
		Error2Edit(A.Error, hDlg, IDC_EDIT3);
	else
	{
		vector2edit(A.X, hDlg, IDC_EDIT3, k, p);
		X<T>.X = A.X;
	}
}

#pragma endregion
#pragma region LDU

template <class T>
void LDU(HWND hDlg)
{
	matrixA<T> A;
	wchar_t path[255];
	GetDlgItemText(hDlg, IDC_EDIT2, path, sizeof(path));
	_bstr_t p(path);
	A.OpenFiles(TypeMatrix::Prof, p);
	matrixLDU<T> LDU;
	LDU.CountLDU(A.al, A.di, A.au, A.ia, TypeMatrix::Prof, A.Error);
	LDU.CountX(A.al, A.di, A.au, A.F, A.X, A.ia);
	if (A.Error)
		Error2Edit(A.Error, hDlg, IDC_EDIT3);
	else
	{
		vector2edit(A.X, hDlg, IDC_EDIT3);
		X<T>.X = A.X;
	}
}

template <class T>
void LDU(HWND hDlg, int k)
{
	matrixA<T> A;
	wchar_t path[255];
	GetDlgItemText(hDlg, IDC_EDIT2, path, sizeof(path));
	_bstr_t p(path);
	if (k == -1) {
		A.OpenFiles(TypeMatrix::Prof, p);
	}else
		A.OpenFiles(TypeMatrix::Prof, p, k);
	matrixLDU<T> LDU;
	LDU.CountLDU(A.al, A.di, A.au, A.ia, TypeMatrix::Prof, A.Error);
	LDU.CountX(A.al, A.di, A.au, A.F, A.X, A.ia);
	if (A.Error)
		Error2Edit(A.Error, hDlg, IDC_EDIT3);
	else
	{
		vector2edit(A.X, hDlg, IDC_EDIT3, k, p);
		X<T>.X = A.X;
	}
}

#pragma endregion


// Глобальные переменные:
HINSTANCE hInst;                                // текущий экземпляр

INT_PTR CALLBACK    DlgProc(HWND, UINT, WPARAM, LPARAM);

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

	CreateDialog(hInstance, MAKEINTRESOURCE(IDD_DIALOG1), NULL, DlgProc);
//	DialogBox(hInstance, MAKEINTRESOURCE(IDD_DIALOG1), NULL, DlgProc);

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_SLAU));

    MSG msg;

    // Цикл основного сообщения:
    while (GetMessage(&msg, nullptr, 0, 0))
    {
        if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    return (int) msg.wParam;
}

INT_PTR CALLBACK DlgProc(HWND hDlg, UINT message, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (message)
    {
    case WM_INITDIALOG:
	{
		HWND hCombo = GetDlgItem(hDlg, IDC_COMBO);
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"LDU");
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Метод Гаусса с выбором главного элемента по столбцам");

		hCombo = GetDlgItem(hDlg, IDC_COMBO2);
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Oдинарная точность (float)");
		SendMessage(hCombo, CB_ADDSTRING, 0, (LPARAM)L"Двойная точность (double)");
		break;
	}
    case WM_COMMAND:
		switch (wParam)
		{
		case IDC_BUTTON2:
		{
			wchar_t path[255];
			GetDlgItemText(hDlg, IDC_EDIT2, path, sizeof(path));
			_bstr_t p(path);
			if (SendDlgItemMessage(hDlg, IDC_COMBO, CB_GETCURSEL, 0, 0))
			{
				if (X<double>.X.size())	X<double>.SaveFile(p);
				else if (X<float>.X.size()) X<float>.SaveFile(p);
			}
			else
				if (X<float>.X.size()) X<float>.SaveFile(p);
				else if (X<double>.X.size()) X<double>.SaveFile(p);
					break;
		}
		case IDC_BUTTON3:
		{
			HWND hCombo = GetDlgItem(hDlg, IDC_COMBO);
			int i = SendMessage(hCombo, CB_GETCURSEL, 0, 0);
			int j = SendMessage(GetDlgItem(hDlg, IDC_COMBO2), CB_GETCURSEL, 0, 0);

			size_t t = GetWindowTextLength(GetDlgItem(hDlg, IDC_EDIT1));
			size_t k;

			if (t)
			{
				k = GetDlgItemInt(hDlg, IDC_EDIT1, 0, 1);
				if (i)
					if (j) Gaus<double>(hDlg, k);
					else   Gaus<float>(hDlg, k);
				else
					if (j) LDU<double>(hDlg, k);
					else   LDU<float>(hDlg, k);
			}
			else
				if (i)
					if (j) Gaus<double>(hDlg);
					else   Gaus<float>(hDlg);
				else
					if (j) LDU<double>(hDlg);
					else   LDU<float>(hDlg);
			break;
		}
		case IDCANCEL:
			EndDialog(hDlg, LOWORD(wParam));
			DestroyWindow(hDlg);
			return (INT_PTR)TRUE;
		default:
			break;
		}
        break;
	case WM_DESTROY:
		PostQuitMessage(0);
		break;
	default:
		return DefWindowProc(hDlg, message, wParam, lParam);
    }
    return (INT_PTR)FALSE;
}
