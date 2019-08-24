#pragma once
#include <list>
#include <vector>
#include "Display.h"
#include "math.h"
#include "Types.h"

class FragmentShader;
class VertexShader;
class Renderer;

class RenderTarget
{
	uint32_t xRes, yRes;
	float *pixelData = nullptr;//[RGBA]
public:
	//Ctor
	RenderTarget(int xRes, int yRes)
	{
		this->xRes = xRes;
		this->yRes = yRes;

		//allocate pixel buffer
		pixelData = new float[xRes * yRes * 4];
	}

	Color GetPixel(int xPos, int yPos)
	{
		if (pixelData == nullptr)
			return Color();
		if (xPos < 0 || xPos >= xRes || yPos < 0 || yPos >= yRes)
			return Color();

		Color retc;
		retc.r = pixelData[(xPos * yRes + yPos) * 4];
		retc.g = pixelData[(xPos * yRes + yPos) * 4 + 1];
		retc.b = pixelData[(xPos * yRes + yPos) * 4 + 2];
		retc.a = pixelData[(xPos * yRes + yPos) * 4 + 3];
		return retc;
	}

	void SetPixel(int xPos, int yPos, Color color)
	{
		if (pixelData == nullptr)
			return ;
		if (xPos < 0 || xPos >= xRes || yPos < 0 || yPos >= yRes)
			return ;

		pixelData[(xPos * yRes + yPos) * 4] = color.r;
		pixelData[(xPos * yRes + yPos) * 4 + 1] = color.g;
		pixelData[(xPos * yRes + yPos) * 4 + 2] = color.b;
		pixelData[(xPos * yRes + yPos) * 4 + 3] = color.a;
	}
};

class FragmentShader {
public:
	Vec3<float> coord_c;
	Vec3<float> normal_w, normal_l;
	Vec2<float> uv;

	virtual Color Run() {
		if(uv.x < 0.5 && uv.y < 0.5)
			return Color(0,0,0);
		if (uv.x < 0.5 && uv.y > 0.5)
			return Color(1,1,1);
		if (uv.x > 0.5 && uv.y < 0.5)
			return Color(1,1,1);
		//if (uv.x > 0.5 && uv.y > 0.5)
			return Color(0,0,0);
	}
};

class VertexShader {

};

class RenderObject
{
public:
	//Transform, vertices, triangles
	RenderObject *parent = nullptr;

	Transform transform;

	FragmentShader fs;

	VertexShader vs;

	std::vector<Triangle3d> triangles;
};

class RenderCamera
{
public:
	Transform transform;
	
	float fieldOfView;//FOV in degree

	bool isPerspective;

	float screenPlaneDist;

	float CalcPlaneDist(Vec2<int> screenRes) {
		return screenPlaneDist = max(screenRes.x, screenRes.y) / tan(DEG(fieldOfView));
	}

	Vec2<float> ProjectPointWorldToScreen(Vec3<float> inPoint)
	{
		inPoint = transform.TransformPointWorldToLocal(inPoint);

		Vec2<float> t;

		t.x = (fabs(inPoint.x) < 1e-6) ? 0 : (screenPlaneDist * inPoint.x / inPoint.y);
		t.y = (fabs(inPoint.z) < 1e-6) ? 0 : (screenPlaneDist * inPoint.z / inPoint.y);

		return t;
	}
};

class Scene
{
public:
	//Contain objects and lightsources

	std::list<RenderObject> objects;
};



class Renderer
{
	Display *targetDisp = nullptr;

public:
	void BindDisplay(Display *display)
	{
		targetDisp = display;
	}

	//Bresenham method draw line
	void DrawLine(Vec2<int> start, Vec2<int> end, Color color)
	{
		if (targetDisp == nullptr)
			return;

		int x, y, dx, dy, dx2, dy2, xstep, ystep, error, index;
		x = start.x;
		y = start.y;
		dx = end.x - start.x;
		dy = end.y - start.y;

		if (dx >= 0) // 从左往右画  
		{
			xstep = 1; // x步进正1  
		}
		else // 从右往左画  
		{
			xstep = -1; // x步进负1  
			dx = -dx; // 取绝对值  
		}

		if (dy >= 0) // 从上往下画  
		{
			ystep = 1; // y步进正1  
		}
		else // 从下往上画  
		{
			ystep = -1; // y步进负1  
			dy = -dy; // 取绝对值  
		}

		dx2 = dx << 1; // 2 * dx  
		dy2 = dy << 1; // 2 * dy  

		if (dx > dy) // 近X轴直线  
		{
			error = dy2 - dx;
			for (index = 0; index <= dx; ++index)
			{
				targetDisp->DrawPixel(x, y, color);
				if (error >= 0)
				{
					error -= dx2;
					y += ystep;
				}
				error += dy2;
				x += xstep;
			}
		}
		else // 近Y轴直线  
		{
			error = dx2 - dy;
			for (index = 0; index <= dy; ++index)
			{
				targetDisp->DrawPixel(x, y, color);
				if (error >= 0)
				{
					error -= dy2;
					x += xstep;
				}
				error += dx2;
				y += ystep;
			}
		}

		return;
	}

	void DrawLine(Vec2<float> start, Vec2<float> end, Color color)
	{
		Vec2<int> st(floor(start.x), floor(start.y)), en(floor(end.x), floor(end.y));
		DrawLine(st, en, color);
	}

	void DrawTriangle(Triangle2d *t, Color color)
	{
		if (t == nullptr || targetDisp == nullptr)
			return;
		//LPRECT rect = targetDisp->GetWndSize();
		Vec2<float> r0(0,0), r1(targetDisp->GetWndWidth(), targetDisp->GetWndHeight());
		Line2d l0(t->v0, t->v1), l1(t->v1, t->v2), l2(t->v2, t->v0);
		if(l0.ClipByRect(r0, r1))
			DrawLine(l0.GetP0(), l0.GetP1(), color);
		if(l1.ClipByRect(r0, r1))
			DrawLine(l1.GetP0(), l1.GetP1(), color);
		if(l2.ClipByRect(r0, r1))
			DrawLine(l2.GetP0(), l2.GetP1(), color);
	}

	void RenderTriangle(Triangle3d source, RenderCamera &cam, FragmentShader &fs /*TODO: Material input*/) {

		Vec2<int> res(targetDisp->GetWndWidth() - 1, targetDisp->GetWndHeight() - 1);
		Vec2<float> offset(targetDisp->GetWndWidth() / 2, targetDisp->GetWndHeight() / 2);
		
		//Rasterize vertices
		Vec2<float> v0 = cam.ProjectPointWorldToScreen(source.V0()) + offset;
		Vec2<float> v1 = cam.ProjectPointWorldToScreen(source.V1()) + offset;
		Vec2<float> v2 = cam.ProjectPointWorldToScreen(source.V2()) + offset;

		//Split triangle into upper and lower halves
		bool upresent = true, lpresent = true;

		if (v0.y < v1.y) {
			//Highest point is v1
			//swap v1 and v0
			Vec2<float> t = v0;
			v0 = v1;
			v1 = t;
		}
		if (v0.y < v2.y) {
			//Highest point is v2
			//swap v2 and v0
			Vec2<float> t = v0;
			v0 = v2;
			v2 = t;
		}
		if (v1.y < v2.y) {
			//Second highest point is v2
			//swap v2 and v0
			Vec2<float> t = v1;
			v1 = v2;
			v2 = t;
		}
		
		source.TransformToLocal(cam.transform);
		source.UpdateProjectionParams();

		if (abs(v0.y - v1.y) < 1e-6)		upresent = false;
		if (abs(v2.y - v1.y) < 1e-6)	lpresent = false;
		bool v1OnLeft = (v2.x > v1.x);
		float dMain = (v2.x - v0.x) / (v2.y - v0.y);
		//Render upper half
		if (upresent) {
			float d01 = (v1.x - v0.x) / (v1.y - v0.y);
			for (int y = floor(min(v0.y, res.y));y > floor(max(v1.y,0)); --y) {
				float xleft, xright;
				if (v1OnLeft) {
					xleft = v0.x + d01*(y - v0.y);
					xright = v0.x + dMain*(y - v0.y);
				}
				else {
					xright = v0.x + d01*(y - v0.y);
					xleft = v0.x + dMain*(y - v0.y);
				}
				for (int x = floor(max(0, xleft)); x < floor(min(res.x, xright)); ++x) {
					//draw pixel according to uv
					fs.coord_c = source.PerspProject(Vec2<float>(x, y) - offset, cam.screenPlaneDist);
					fs.uv = source.GetUV(fs.coord_c);
					targetDisp->DrawPixel(x, y, fs.Run());
					//targetDisp->DrawPixel(x, y, Color(1, 1, 1));
				}
			}
		}

		//Render lower half
		if (lpresent) {
			float d12 = (v1.x - v2.x) / (v1.y - v2.y);
			for (int y = floor(max(v2.y, 0)); y <= floor(min(v1.y, res.y)); ++y) {
				float xleft, xright;
				if (v1OnLeft) {
					xleft = v2.x + d12*(y - v2.y);
					xright = v2.x + dMain*(y - v2.y);
				}
				else {
					xright = v2.x + d12*(y - v2.y);
					xleft = v2.x + dMain*(y - v2.y);
				}
				for (int x = floor(max(0, xleft)); x < floor(min(res.x, xright)); ++x) {
					//draw pixel according to uv
					fs.coord_c = source.PerspProject(Vec2<float>(x, y) - offset, cam.screenPlaneDist);
					fs.uv = source.GetUV(fs.coord_c);
					targetDisp->DrawPixel(x, y, fs.Run());
					//targetDisp->DrawPixel(x, y, Color(1,1,1));
				}
			}
		}

	}

	void RenderSceneToDisplay(Scene *scene, RenderCamera *cam)
	{
		if (scene == nullptr || cam == nullptr)
			return;

		//Calculate offset
		Vec2<float> offset(
			targetDisp->GetWndWidth() / 2,
			targetDisp->GetWndHeight() / 2
		);

		Vec2<int> ScreenRes(targetDisp->GetWndWidth(), targetDisp->GetWndHeight());

		cam->CalcPlaneDist(ScreenRes);

		for (RenderObject &obj : scene->objects)
		{
			for (Triangle3d t : obj.triangles)
			{
				obj.fs.normal_l = t.Normal();
				Triangle2d t2d;
				t.TransformToWorld(obj.transform);
				obj.fs.normal_w = t.Normal();
				//Project points
				t2d.v0 = cam->ProjectPointWorldToScreen(t.V0()) + offset;
				t2d.v1 = cam->ProjectPointWorldToScreen(t.V1()) + offset;
				t2d.v2 = cam->ProjectPointWorldToScreen(t.V2()) + offset;
				DrawTriangle(&t2d, Color(1, 0, 1));
				RenderTriangle(t, *cam, obj.fs);
			}
		}
	}


};

