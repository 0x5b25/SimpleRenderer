//SimpleRenderer.cpp: 定义应用程序的入口点。
//

#include "stdafx.h"
#include "SimpleRenderer.h"
#include "Display.h"
#include "Render.h"
#include "cstdio"
#include "iostream"
#include "ctime"

#define MAX_LOADSTRING 100

// 全局变量: 
HINSTANCE hInst;                                // 当前实例
WCHAR szTitle[MAX_LOADSTRING];                  // 标题栏文本
WCHAR szWindowClass[MAX_LOADSTRING];            // 主窗口类名

static bool bIsRunning;//Run flag, program will quit if false
static HWND hMainWnd;//Main render window
static RECT MainWndRect;//Main window size

Scene scene;//Create scene
RenderCamera cam;//Create camera


// 此代码模块中包含的函数的前向声明: 
ATOM                MyRegisterClass(HINSTANCE hInstance);
BOOL                InitInstance(HINSTANCE, int);
LRESULT CALLBACK    WndProc(HWND, UINT, WPARAM, LPARAM);

void SetupScene(Scene* scene, RenderCamera* cam)
{
	//Triangles
	//Triangle3d t;

	RenderObject o;

	//Triangle 1
	//t.v0 = Vec3<float>(1, 0, 0);
	//t.v1 = Vec3<float>(0, 0, 1);
	//t.v2 = Vec3<float>(0, 0, 0);
	o.triangles.push_back(Triangle3d(Vec3<float>(1, 0, 0), Vec3<float>(0, 0, 1), Vec3<float>(0, 0, 0),
										Vec2<float>(1,0), Vec2<float>(0, 1), Vec2<float>(0, 0)));

	//Triangle 2
	//t.v0 = Vec3<float>(1, 0, 0);
	//t.v1 = Vec3<float>(0, 0, 0);
	//t.v2 = Vec3<float>(0, 1, 0);
	o.triangles.push_back(Triangle3d(Vec3<float>(1, 0, 0), Vec3<float>(0, 0, 0), Vec3<float>(0, 1, 0),
										Vec2<float>(1, 0), Vec2<float>(0, 0), Vec2<float>(0, 1)));

	//Triangle 3
	//t.v0 = Vec3<float>(0, 1, 0);
	//t.v1 = Vec3<float>(0, 0, 0);
	//t.v2 = Vec3<float>(0, 0, 1);
	o.triangles.push_back(Triangle3d(Vec3<float>(0, 1, 0), Vec3<float>(0, 0, 0), Vec3<float>(0, 0, 1),
										Vec2<float>(1, 0), Vec2<float>(0, 0), Vec2<float>(0, 1)));

	o.transform.Location = Vec3<float>(0, 1, 0);
	o.transform.Rotation = Vec3<float>(0, 0, 0);
	o.transform.Scale = Vec3<float>(1, 1, 1);

	scene->objects.push_back(o);

	//Camera
	cam->fieldOfView = 60;
	cam->isPerspective = true;
	cam->transform.Location = Vec3<float>(0, 0, 0.5);
	cam->transform.Rotation = Vec3<float>(DEG(0), 0, 0);
	cam->transform.Scale = Vec3<float>(1, 1, 1);
}

int APIENTRY wWinMain(_In_ HINSTANCE hInstance,
                     _In_opt_ HINSTANCE hPrevInstance,
                     _In_ LPWSTR    lpCmdLine,
                     _In_ int       nCmdShow)
{
    UNREFERENCED_PARAMETER(hPrevInstance);
    UNREFERENCED_PARAMETER(lpCmdLine);

    // TODO: 在此放置代码。

    // 初始化全局字符串
    LoadStringW(hInstance, IDS_APP_TITLE, szTitle, MAX_LOADSTRING);
    LoadStringW(hInstance, IDC_SIMPLERENDERER, szWindowClass, MAX_LOADSTRING);
    MyRegisterClass(hInstance);

    // 执行应用程序初始化: 
    if (!InitInstance (hInstance, nCmdShow))
    {
        return FALSE;
    }

    HACCEL hAccelTable = LoadAccelerators(hInstance, MAKEINTRESOURCE(IDC_SIMPLERENDERER));

    MSG msg;

	//Create display
	Display disp;
	Renderer rnd;

	

	//Bind window
	disp.BindWindow(hMainWnd);
	rnd.BindDisplay(&disp);

	AllocConsole();
	freopen("conout$", "w", stdout);
	printf("hello hplonline!-_-\n");
	std::cout << "i'm cout" << std::endl;
	freopen("conout$", "w", stderr);
	std::cerr << "i'm cerr" << std::endl;

	SetupScene(&scene, &cam);

	//Mark run flag
	bIsRunning = true;

	DWORD prevTime = GetTickCount();
	int delayCycle = 0;
    //Main loop
	while (bIsRunning)
	{
		// 主消息循环: 
		if (PeekMessage(&msg, nullptr, 0, 0, true))
		{
			//Got a pending message
			if (!TranslateAccelerator(msg.hwnd, hAccelTable, &msg))
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
		}
		
		//Update frame info
		disp.UpdateWnd();

		//Clear display buffer
		//disp.FillBackgroundGreyscale(0);
		disp.FillBackground(Color(0, 0, 0));
		//Vec2<int> st(5,7), en(180,900);
		//rnd.DrawLine(en, st, Color(1, 0.5, 0.7));
		rnd.RenderSceneToDisplay(&scene, &cam);
		
		//Display
		disp.PushBuffer();
		
		if (delayCycle <= 0) {
			delayCycle = 10;
			DWORD currTime = GetTickCount();
			float frameTime = (currTime - prevTime)/10;
			float fps = 1000.0 / frameTime;
			prevTime = currTime;
			printf("fps:%f frametime:%f\n", fps, frameTime);
		}
		else {
			delayCycle--;
		}
	}

    return (int) msg.wParam;
}

ATOM MyRegisterClass(HINSTANCE hInstance)
{
    WNDCLASSEXW wcex;

    wcex.cbSize = sizeof(WNDCLASSEX);

    wcex.style          = CS_HREDRAW | CS_VREDRAW;
    wcex.lpfnWndProc    = WndProc;
    wcex.cbClsExtra     = 0;
    wcex.cbWndExtra     = 0;
    wcex.hInstance      = hInstance;
    wcex.hIcon          = LoadIcon(hInstance, MAKEINTRESOURCE(IDI_SIMPLERENDERER));
    wcex.hCursor        = LoadCursor(nullptr, IDC_ARROW);
    wcex.hbrBackground  = (HBRUSH)(COLOR_WINDOW+1);
    wcex.lpszMenuName   = NULL;
    wcex.lpszClassName  = szWindowClass;
    wcex.hIconSm        = LoadIcon(wcex.hInstance, MAKEINTRESOURCE(IDI_SMALL));

    return RegisterClassExW(&wcex);
}

BOOL InitInstance(HINSTANCE hInstance, int nCmdShow)
{
   hInst = hInstance; // 将实例句柄存储在全局变量中

   hMainWnd = CreateWindowW(szWindowClass, szTitle, WS_OVERLAPPEDWINDOW,
      CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, nullptr, nullptr, hInstance, nullptr);

   if (!hMainWnd)
   {
      return FALSE;
   }

   ShowWindow(hMainWnd, nCmdShow);
   UpdateWindow(hMainWnd);

   return TRUE;
}

LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
    switch (message)
    {
	case WM_KEYDOWN:
	{
		switch (wParam)
		{
		case 'A':
		{
			cam.transform.Location.x+=0.01;
		}break;
		case 'D':
		{
			cam.transform.Location.x -= 0.01;
		}break;
		case 'Q':
		{
			cam.transform.Rotation.x -= 0.01;
		}break;
		case 'E':
		{
			cam.transform.Rotation.x += 0.01;
		}break;
		case 'W':
		{
			scene.objects.front().transform.Rotation.z += 0.1;
		}break;
		case 'S':
		{
			scene.objects.front().transform.Rotation.z -= 0.1;
		}break;
		}
	}
	break;
    case WM_DESTROY:
		{
			PostQuitMessage(0);
			bIsRunning = false;
			DestroyWindow(hWnd);
		}
        break;
    default:
        return DefWindowProc(hWnd, message, wParam, lParam);
    }
    return 0;
}