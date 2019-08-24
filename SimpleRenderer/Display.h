#pragma once

#include "stdafx.h"
#include "Types.h"



class Display
{
	//origin point is at bottomleft of the screen


	//Window properties
	//window size with frame
	RECT wndSize;
	//actual draw area size
	RECT clientSize;
	HWND hWnd = nullptr;
	HDC hdcWnd = nullptr;
	POINT originPoint;

	//Frame buffer
	HDC hdcMem = nullptr;
	HBITMAP bmp = nullptr;
	BITMAPINFO bmpInfo;
	long bytesInLine;
	BYTE *pixelData = nullptr;

	//Size info
	long clientHeight = 0;
	long clientWidth = 0;
	long wndHeight = 0;
	long wndWidth = 0;


public:

	~Display() {
		//hWnd and hdcWnd are got from host window therefore no need to call DeleteXX() or DestroyXX()
		//BUT hdcWnd need to be destroyed!!

		//Assume that hWnd is not null
		if (hdcWnd != nullptr)
		{
			ReleaseDC(hWnd, hdcWnd);
		}

		//clean up framebuffer
		if (hdcMem != nullptr)
		{
			DeleteDC(hdcMem);
		}

		if (bmp != nullptr)
		{
			DeleteObject(bmp);
		}

		if (pixelData != nullptr)
			delete[] pixelData;
	}

	bool BindWindow(HWND hWindow)
	{
		if (hWindow == nullptr)
			return false;

		
		if (hdcWnd != nullptr)
		{
			ReleaseDC(hWnd, hdcWnd);
		}

		// Set window handle
		hWnd = hWindow;

		hdcWnd = GetWindowDC(hWindow);

		if (hdcWnd == nullptr)
			return false;

		wndSize = { -1, -1, -1, -1 };

		UpdateWnd();

		return true;
	}

	void UpdateWnd(bool forceUpdate = false)
	{
		RECT prevSize = clientSize;
		//Get window size
		GetWindowRect(hWnd, &wndSize);
		GetClientRect(hWnd, &clientSize);

		if (!forceUpdate && 
			prevSize.bottom - prevSize.top == clientSize.bottom - clientSize.top &&
			prevSize.right - prevSize.left == clientSize.right - clientSize.left)
			
			return;
		//Calculate client window origin point
		originPoint = { 0,0 };

		ClientToScreen(hWnd, &originPoint);

		originPoint.x -= wndSize.left;
		originPoint.y -= wndSize.top;

		clientHeight = clientSize.bottom - clientSize.top + 1;
		clientWidth = clientSize.right - clientSize.left + 1;
		wndHeight = wndSize.bottom - wndSize.top + 1;
		wndWidth = wndSize.right - wndSize.left + 1;

		// Set framebuffer
		if (hdcMem != nullptr)
		{
			DeleteDC(hdcMem);
		}
		hdcMem = CreateCompatibleDC(hdcWnd);

		if (bmp != nullptr)
		{
			DeleteObject(bmp);
		}
		bmp = CreateCompatibleBitmap(hdcWnd, clientSize.right - clientSize.left + 1, clientSize.bottom - clientSize.top + 1);
		//Bind bmp to framebuffer
		SelectObject(hdcMem, bmp);

		//BITMAPINFO bmpInfo;
		bmpInfo.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
		bmpInfo.bmiHeader.biWidth = clientSize.right - clientSize.left + 1;
		bmpInfo.bmiHeader.biHeight =  clientSize.bottom - clientSize.top + 1;
		bmpInfo.bmiHeader.biPlanes = 1;
		bmpInfo.bmiHeader.biBitCount = 24;
		bmpInfo.bmiHeader.biCompression = BI_RGB;
		bmpInfo.bmiHeader.biSizeImage = 0;
		bmpInfo.bmiHeader.biXPelsPerMeter = 3000;
		bmpInfo.bmiHeader.biYPelsPerMeter = 3000;
		bmpInfo.bmiHeader.biClrUsed = 0;
		bmpInfo.bmiHeader.biClrImportant = 0;

		//每行字节数，4字节对齐  
		bytesInLine = ((clientSize.right - clientSize.left + 1) * 3 + 3) / 4 * 4;

		//pixel data
		if (pixelData != nullptr)
			delete[] pixelData;
		pixelData = new BYTE[bytesInLine * (clientSize.bottom - clientSize.top + 1)];

		//Fill black
		FillBackgroundGreyscale(0);
	}

	void FillBackground(Color color)
	{
		for (int x = 0; x < GetWndWidth(); x++)
		{
			for (int y = 0; y < GetWndHeight(); y++)
			{
				//R
				pixelData[y * bytesInLine + x * 3 + 2] = 255 * color.r;
				//G
				pixelData[y * bytesInLine + x * 3 + 1] = 255 * color.g;
				//B
				pixelData[y * bytesInLine + x * 3] = 255 * color.b;
			}
		}
	}

	void FillBackgroundGreyscale(float g)
	{
		memset(pixelData, 255*g, bytesInLine * (clientSize.bottom - clientSize.top + 1));
	}
	//Get actual draw size
	LPRECT GetWndSize() { return &clientSize; }
	//Get actual draw height
	int GetWndHeight() { return clientHeight; }
	//Get actual draw width
	int GetWndWidth() { return clientWidth; }
	//Get total window height
	int GetWndHeightWithFrame() { return wndHeight; }
	//Get total window width
	int GetWndWidthWithFrame() { return wndWidth; }

	bool CheckPosition(int x, int y)
	{
		//if (pixelData == nullptr || bmp == nullptr || hdcMem == nullptr || hdcWnd == nullptr)
			//return false;
		if (x < 0 || x >= GetWndWidth() || y < 0 || y >= GetWndHeight())
			return false;
		return true;
	}

	void DrawPixel(int x, int y, Color color)
	{
		if (hdcWnd == nullptr)
			return;
		if (!CheckPosition(x, y))
			return;
		//R
		pixelData[y * bytesInLine + x * 3 + 2] = 255 * color.r;
		//G
		pixelData[y * bytesInLine + x * 3 + 1] = 255 * color.g;
		//B
		pixelData[y * bytesInLine + x * 3] = 255 * color.b;
	}

	void PushBuffer()
	{
		if (hdcWnd == nullptr)
			return;

		//Copy pixel data to framebuffer
		SetDIBits(hdcMem, bmp, 0, GetWndHeight(), pixelData, &bmpInfo, DIB_RGB_COLORS);

		//Push framebuffer to screen
		BitBlt(hdcWnd, originPoint.x, originPoint.y, GetWndWidth(), GetWndHeight(), hdcMem, 0, 0, SRCCOPY);		
	}
};