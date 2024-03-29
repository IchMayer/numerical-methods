#pragma once

#include "resource.h"
#include "stdafx.h"

template<class T>
struct endElem
{
	vector<T> X;
	bool SaveFile(char *path)
	{
		char Path[255];
		strcpy_s(Path, path);
		strcat_s(Path,"X.txt\0");
		ofstream oX(Path);
		if (!oX)
		{
			cout << "Error: File X";
			return false;
		}
		for (size_t i = 0; i < X.size(); i++)
			oX << X[i] << endl;
		return true;
	}
};

template<class T>
endElem<T> X;

template<class T>
void vector2edit(vector<T> X, HWND hDlg, size_t Item)
{
	SetDlgItemText(hDlg, Item, L"");
	for (size_t i = 0; i < X.size(); i++)
	{
		wchar_t w[20] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
		if (X[i] > 1000000000)
			swprintf_s(w, L"%e\r\n", X[i]);
		else
			swprintf_s(w, L"%lf\r\n", X[i]);
		SendDlgItemMessage(hDlg, Item, EM_REPLACESEL, TRUE, (LONG)w);
	}
}

void Error2Edit(size_t Error, HWND hDlg, size_t Item)
{
	SetDlgItemText(hDlg, Item, L"");
	switch (Error)
	{
	case 1:
		SendDlgItemMessage(hDlg, Item, EM_REPLACESEL, TRUE, (LPARAM)L"Не удалось\r\n открыть\r\n файлы\r\n");
		break;
	case 2:
		SendDlgItemMessage(hDlg, Item, EM_REPLACESEL, TRUE, (LPARAM)L"Не удалось\r\n сохранить\r\n в файл\r\n");
		break;
	case 3:
		SendDlgItemMessage(hDlg, Item, EM_REPLACESEL, TRUE, (LPARAM)L"无法计算矩阵");
		break;
	default:
		break;
	}
}

template<class T>
void vector2edit(vector<T> X, HWND hDlg, size_t Item, int k, char *path)
{
	char Path[255];
	strcpy_s(Path, path);
	strcat_s(Path, "Result.txt\0");
	ifstream ire(Path);
	if (!ire)
		return;

	SetDlgItemText(hDlg, Item, L"");
	SetDlgItemText(hDlg, IDC_EDIT10, L"");
	SetDlgItemText(hDlg, IDC_EDIT11, L"");
	for (size_t i = 0; i < X.size(); i++)
	{
		long double d;
		ire >> d;
		wchar_t w[40] = { '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0' };
		swprintf_s(w, is_same < T, double >::value ? L"%.15lf \r\n" : L"%.7lf \r\n", X[i]);
		SendDlgItemMessage(hDlg, Item, EM_REPLACESEL, TRUE, (LONG)w);

		for (size_t i = 0; i < 40; i++)
			w[i] = '\0';
		swprintf_s(w, is_same<T, double>::value ? L"%.15lf \r\n" : L"%.7lf \r\n", d);
		SendDlgItemMessage(hDlg, IDC_EDIT10, EM_REPLACESEL, TRUE, (LONG)w);

		for (size_t i = 0; i < 40; i++)
			w[i] = '\0';
		swprintf_s(w, L"%e \r\n", d - X[i]);
		SendDlgItemMessage(hDlg, IDC_EDIT11, EM_REPLACESEL, TRUE, (LONG)w);
	}

	ire.close();
}