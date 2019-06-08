// Project3.cpp : Определяет точку входа для приложения.
//

#include "stdafx.h"
#include "basewin.h"
#include "Lab4.h"

chm4 lab;

template <class T> void SafeRelease(T **ppT)
{
	if (*ppT)
	{
		(*ppT)->Release();
		*ppT = NULL;
	}
}

class MainWindow : public BaseWindow<MainWindow>
{
	ID2D1Factory            *pFactory;
	IDWriteFactory			*pDWriteFactory;
	ID2D1HwndRenderTarget   *pRenderTarget;
	ID2D1BitmapRenderTarget *pRenderBitmap;
	ID2D1SolidColorBrush    *pBrush;
	IDWriteTextFormat		*pTextFormat;
	ID2D1Bitmap				*pBitmap;


	HRESULT CreateGraphicsResources();
	void    DiscardGraphicsResources();
	void    OnPaint();

public:

	MainWindow() : pFactory(NULL), pRenderTarget(NULL), pBrush(NULL)
	{
	}
	
	PCWSTR  ClassName() const { return L"4m4 window"; }
	LRESULT HandleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam);
};

// Recalculate drawing layout when the size of the window changes.

HRESULT MainWindow::CreateGraphicsResources()
{
	HRESULT hr = S_OK;
	if (pRenderTarget == NULL)
	{
		RECT rc;
		GetClientRect(m_hwnd, &rc);

		D2D1_SIZE_U size = D2D1::SizeU(rc.right, rc.bottom);

		hr = pFactory->CreateHwndRenderTarget(
			D2D1::RenderTargetProperties(),
			D2D1::HwndRenderTargetProperties(m_hwnd, size),
			&pRenderTarget);


		if (SUCCEEDED(hr))
		{
			const D2D1_COLOR_F color = D2D1::ColorF(1.0f, 1.0f, 0.0f);
			hr = pRenderTarget->CreateSolidColorBrush(color, &pBrush);
		}


		static const WCHAR msc_fontName[] = L"Verdana";
		static const FLOAT msc_fontSize = 24;
		HRESULT hr = S_OK;

		if (SUCCEEDED(hr))
		{

			// Create a DirectWrite factory.
			hr = DWriteCreateFactory(
				DWRITE_FACTORY_TYPE_SHARED,
				__uuidof(pDWriteFactory),
				reinterpret_cast<IUnknown **>(&pDWriteFactory)
			);
		}
		if (SUCCEEDED(hr))
		{
			// Create a DirectWrite text format object.
			hr = pDWriteFactory->CreateTextFormat(
				msc_fontName,
				NULL,
				DWRITE_FONT_WEIGHT_NORMAL,
				DWRITE_FONT_STYLE_NORMAL,
				DWRITE_FONT_STRETCH_NORMAL,
				msc_fontSize,
				L"", //locale
				&pTextFormat
			);
		}

		hr = pRenderTarget->CreateCompatibleRenderTarget(D2D1::SizeF(1000.0f, 750.0f),&pRenderBitmap);
		hr = pRenderTarget->CreateBitmap(D2D1::SizeU(1000,750), D2D1::BitmapProperties(), &pBitmap);

	}
	return hr;
}

void MainWindow::DiscardGraphicsResources()
{
	SafeRelease(&pRenderTarget);
	SafeRelease(&pBrush);
}

void MainWindow::OnPaint()
{
	HRESULT hr = CreateGraphicsResources();
	if (SUCCEEDED(hr))
	{
		pRenderBitmap->BeginDraw();
		pRenderBitmap->Clear(D2D1::ColorF(D2D1::ColorF::White));
		float f;
		pBrush->SetColor(D2D1::ColorF(1.0f, 1.0f, 1.0f));
		for (double i = lab.Move[1].x, i1 = 0; i < lab.Move[1].y; i += 0.02, i1++)
		{
			for (double j = lab.Move[2].x, j1 = 0; j > lab.Move[2].y; j -= 0.02, j1++)
			{
				f = lab.GetNorm(i, j) / lab.FMax;
				if(f < 0)
				{
					pBrush->SetColor(D2D1::ColorF(0, 1, 0));
					pRenderBitmap->FillRectangle(D2D1::RectF(i1 - 1, j1 - 1, i1 + 1, j1 + 1), pBrush);
				}
				else
				{
		//			f *= 3;
					f = pow(f, 0.333);
					f = (int)(f * 20) / 20.0;
					pBrush->SetColor(D2D1::ColorF(f, f, f));
					pRenderBitmap->FillRectangle(D2D1::RectF(i1, j1, i1 + 1, j1 + 1), pBrush);
				}
			}
		}
		pBrush->SetColor(D2D1::ColorF(1, 0, 0));
		pRenderBitmap->FillRectangle(D2D1::RectF(0, lab.Move[0].y - 1, 1000, lab.Move[0].y + 1), pBrush);
		pRenderBitmap->FillRectangle(D2D1::RectF(lab.Move[0].x - 1, 0, lab.Move[0].x + 1, 750), pBrush);

		pBrush->SetColor(D2D1::ColorF(0, 0, 1));
		D2D1_POINT_2F o, o1;
		o = lab.Target[0];
		for (size_t i = 1; i < lab.Target.size(); i++)
		{
			o1 = o;
			o = lab.Target[i];
			pRenderBitmap->DrawLine(D2D1::Point2F(o1.x *50 + lab.Move[0].x, o1.y * -50 + lab.Move[0].y),
				D2D1::Point2F(o.x * 50 + lab.Move[0].x, o.y * -50 + lab.Move[0].y), pBrush, 2);
		}
		if (lab.dx.size() == 2 && lab.NormF / lab.NormF0 < lab.e1)
		{
			o.x = lab.dx[0];
			o.y = lab.dx[1];
			pBrush->SetColor(D2D1::ColorF(0, 0, 0.5));
			pRenderBitmap->DrawLine(D2D1::Point2F(o1.x * 50 + lab.Move[0].x, o1.y * -50 + lab.Move[0].y),
				D2D1::Point2F(o.x * 50 + lab.Move[0].x, o.y * -50 + lab.Move[0].y), pBrush, 1);
		}

		pBrush->SetColor(D2D1::ColorF(1, 0, 1));
		pRenderBitmap->FillEllipse(D2D1::Ellipse(D2D1::Point2F(lab.x[0] * 50 + lab.Move[0].x, lab.x[1] * -50 + lab.Move[0].y), 4, 4), pBrush);

		pBrush->SetColor(D2D1::ColorF(0, 0, 0));
		wchar_t Coor[100];
		swprintf_s(Coor, L"x: %lf\ny: %lf\nk: %i\n\0", lab.x[0], lab.x[1], lab.k - 1);
		pRenderBitmap->DrawTextW(Coor, wcslen(Coor), pTextFormat, D2D1::RectF(10, 10, 300, 500), pBrush);
		swprintf_s(Coor, L"F: %e\nB: %e\0", lab.NormF, lab.b);
		pRenderBitmap->DrawTextW(Coor, wcslen(Coor), pTextFormat, D2D1::RectF(10, 650, 300, 740), pBrush);

		//pRenderTarget->FillRectangle(rectangle, pBrush);
		pRenderBitmap->GetBitmap(&pBitmap);
		pRenderBitmap->EndDraw();


		PAINTSTRUCT ps;
		BeginPaint(m_hwnd, &ps);
		pRenderTarget->BeginDraw();

		pRenderTarget->DrawBitmap(pBitmap, D2D1::RectF(0, 0, 1000, 750), 1, D2D1_BITMAP_INTERPOLATION_MODE_NEAREST_NEIGHBOR, D2D1::RectF(0, 0, 1000, 750));

		hr = pRenderTarget->EndDraw();
		if (FAILED(hr) || hr == D2DERR_RECREATE_TARGET)
		{
			DiscardGraphicsResources();
		}
		EndPaint(m_hwnd, &ps);
	}
}

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE, PWSTR, int nCmdShow)
{
	MainWindow win;
	lab.Init();
	lab.TwoCicle();

	if (!win.Create(L"4m4", WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU ))
		return 0;
	ShowWindow(win.Window(), nCmdShow);

	// Run the message loop.

	MSG msg = { };
	while (GetMessage(&msg, NULL, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	return 0;
}

void MoveSkipY(int o, double l)
{
	lab.Move[0].y -= o;
	lab.Move[2].x -= l;
	lab.Move[2].y -= l;
}

void MoveSkipX(int o, double l)
{
	lab.Move[0].x -= o;
	lab.Move[1].x += l;
	lab.Move[1].y += l;
}

LRESULT MainWindow::HandleMessage(UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	switch (uMsg)
	{
	case WM_CREATE:
		if (FAILED(D2D1CreateFactory(
			D2D1_FACTORY_TYPE_SINGLE_THREADED, &pFactory)))
		{
			return -1;  // Fail CreateWindowEx.
		}
		return 0;
	case WM_DESTROY:
		DiscardGraphicsResources();
		SafeRelease(&pFactory);
		PostQuitMessage(0);
		return 0;

	case WM_PAINT:
		OnPaint();
		return 0;

	case WM_KEYDOWN:
		switch (wParam)
		{
		case VK_UP:
			MoveSkipY(-50, -1);
			InvalidateRect(m_hwnd, NULL, FALSE);
			break;
		case VK_DOWN:
			MoveSkipY(50, 1);
			InvalidateRect(m_hwnd, NULL, FALSE);
			break;
		case VK_LEFT:
			MoveSkipX(-50, -1);
			InvalidateRect(m_hwnd, NULL, FALSE);
			break;
		case VK_RIGHT:
			MoveSkipX(50, 1);
			InvalidateRect(m_hwnd, NULL, FALSE);
			break;
		default:
			break;
		}
		return 0;
	}
	return DefWindowProc(m_hwnd, uMsg, wParam, lParam);
}