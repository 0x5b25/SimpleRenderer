#pragma once
#include "stdafx.h"
#include "cstdint"
#include "cmath"

#define DEG(x) x/180.0*3.14159
inline float FInf() { int z = 0; return 1.0 / z; }
inline float FNInf() { int z = 0; return -1.0 / z; }
inline float FNan() { return sqrt(-1); }

class Color
{
public:
	float r, g, b, a;
	//Ctor
	Color()
	{
		this->r = 1;
		this->g = 1;
		this->b = 1;
		this->a = 1;
	}
	Color(float r, float g, float b, float a)
	{
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = a;
	}
	Color(float r, float g, float b)
	{
		this->r = r;
		this->g = g;
		this->b = b;
		this->a = 1;
	}
};

template <typename T>
class Vec2
{
public: 
	T x, y;
	//Ctors
	Vec2()
	{
		x = 0, y = 0;
	}

	Vec2(T x, T y)
	{
		this->x = x, this->y = y;
	}


	T CrossProduct(Vec2<T> v) {
		return x*v.y - v.x * y;
	}

	//Operators
	Vec2 operator- (Vec2 b)
	{
		Vec2 t;
		t.x = x - b.x;
		t.y = y - b.y;
		return t;
	}

	Vec2 operator+ (Vec2 b)
	{
		Vec2 t;
		t.x = x + b.x;
		t.y = y + b.y;
		return t;
	}

	T operator* (Vec2 b)
	{

		return x*b.x + y*b.y;
	}

	Vec2 operator* (int b)
	{
		Vec2 t;
		t.x = x * b;
		t.y = y * b;
		return t;

	}

	Vec2 operator* (float b)
	{
		Vec2 t;
		t.x = x * b;
		t.y = y * b;
		return t;

	}

	
};

template <typename T>
class Vec3
{
public:
	T x, y, z;

	Vec3()
	{
		x = 0, y = 0, z = 0;
	}

	Vec3(T x, T y, T z)
	{
		this->x = x, this->y = y, this->z = z;
	}
	//Functions
	void Normalize() {
		float l = Length();
		x /= l;
		y /= l;
		z /= l;
	}
	float Length() { return sqrt(x*x + y*y + z*z); }
	T LengthSquared() { return x*x + y*y + z*z; }
	Vec3<T> CrossProduct(Vec3<T> v) {
		Vec3<T> t;
		t.x = (y*v.z - v.y*z);
		t.y = (z*v.x - v.z*x);
		t.z = (x*v.y - v.x*y);
		return t;
	}
	//Operators
	Vec3<T> operator- (Vec3<T> b)
	{
		Vec3<T> t;
		t.x = x - b.x;
		t.y = y - b.y;
		t.z = z - b.z;
		return t;
	}

	Vec3<T> operator+ (Vec3<T> b)
	{
		Vec3<T> t;
		t.x = x + b.x;
		t.y = y + b.y;
		t.z = z + b.z;
		return t;
	}

	T operator* (Vec3<T> b)
	{
		return x*b.x + y*b.y + z*b.z;
	}

	Vec3<T> operator* (int b)
	{
		Vec3<T> t;
		t.x = x * b;
		t.y = y * b;
		t.z = z * b;
		return t;

	}

	Vec3<T> operator* (float b)
	{
		Vec3<T> t;
		t.x = x * b;
		t.y = y * b;
		t.z = z * b;
		return t;
	}

	Vec3<T> operator/ (float b)
	{
		Vec3<T> t;
		t.x = x / b;
		t.y = y / b;
		t.z = z / b;
		return t;
	}
};

template <typename T>
class Vec4
{
public:
	T x, y, z, w;

	Vec4()
	{
		x = 0, y = 0, z = 0, w = 0;
	}

	Vec4(T x, T y, T z, T w)
	{
		this->x = x, this->y = y, this->z = z, this->w = w;
	}
	//Functions
	float Length() { return sqrt<T>(x*x + y*y + z*z + w*w); }
	T LengthSquared() { return x*x + y*y + z*z + w*w; }
	//Operators
	Vec4 operator- (Vec4 b)
	{
		Vec4 t;
		t.x = x - b.x;
		t.y = y - b.y;
		t.z = z - b.z;
		t.w = w - b.w;
		return t;
	}

	Vec4 operator+ (Vec4 b)
	{
		Vec4 t;
		t.x = x + b.x;
		t.y = y + b.y;
		t.z = z + b.z;
		t.w = w + b.w;
		return t;
	}

	T operator* (Vec4 b)
	{
		return x*b.x + y*b.y + z*b.z + w*b.w;
	}

	Vec3 operator* (int b)
	{
		Vec3 t;
		t.x = x * b;
		t.y = y * b;
		t.z = z * b;
		t.w = w * b;
		return t;

	}

	Vec3 operator* (float b)
	{
		Vec3 t;
		t.x = x * b;
		t.y = y * b;
		t.z = z * b;
		t.w = w * b;
		return t;

	}
};

class Mat3
{
public:
	float a1, a2, a3;
	float b1, b2, b3;
	float c1, c2, c3;

	Mat3()
	{
		this->a1 = 1, this->a2 = 0, this->a3 = 0;
		this->b1 = 0, this->b2 = 1, this->b3 = 0;
		this->c1 = 0, this->c2 = 0, this->c3 = 1;
	}

	Mat3(float a1, float a2, float a3,
		float b1, float b2, float b3,
		float c1, float c2, float c3)
	{
		this->a1 = a1, this->a2 = a2, this->a3 = a3;
		this->b1 = b1, this->b2 = b2, this->b3 = b3;
		this->c1 = c1, this->c2 = c2, this->c3 = c3;
	}

	Vec3<float> operator* (Vec3<float> b)
	{
		Vec3<float> t;
		t.x = a1 * b.x + a2 * b.y + a3 * b.z;
		t.y = b1 * b.x + b2 * b.y + b3 * b.z;
		t.z = c1 * b.x + c2 * b.y + c3 * b.z;
		return t;
	}

	Mat3 operator* (Mat3 b)
	{
		Mat3 t;
		t.a1 = a1 * b.a1 + a2 * b.b1 + a3 * b.c1;
		t.a2 = a1 * b.a2 + a2 * b.b2 + a3 * b.c2;
		t.a3 = a1 * b.a3 + a2 * b.b3 + a3 * b.c3;
		t.b1 = b1 * b.a1 + b2 * b.b1 + b3 * b.c1;
		t.b2 = b1 * b.a2 + b2 * b.b2 + b3 * b.c2;
		t.b3 = b1 * b.a3 + b2 * b.b3 + b3 * b.c3;
		t.c1 = c1 * b.a1 + c2 * b.b1 + c3 * b.c1;
		t.c2 = c1 * b.a2 + c2 * b.b2 + c3 * b.c2;
		t.c3 = c1 * b.a3 + c2 * b.b3 + c3 * b.c3;
		return t;
	}
};

class Mat4
{

};

// Coordinate system:
// Z
// |       Y    
// |      /
// |     /
// |    /
// |   /
// |  /
// | /
// |/_________________ X
class Transform
{
public:
	//Rotation order: x, y, z
	Vec3<float> Location, Rotation, Scale;

	Vec3<float> TransformPointLocalToWorld(Vec3<float> inPoint)
	{
		//Scale
		inPoint.x *= Scale.x;
		inPoint.y *= Scale.y;
		inPoint.z *= Scale.z;

		//Rotate
		//x
		Mat3 rx(
			1,	0,					0,
			0,	cos(Rotation.x),	-sin(Rotation.x),
			0,	sin(Rotation.x),	cos(Rotation.x)
		);
		//y
		Mat3 ry(
			cos(Rotation.y),	0,	sin(Rotation.y),
			0,					1,	0,
			-sin(Rotation.y),	0,	cos(Rotation.y)
		);
		//z
		Mat3 rz(
			cos(Rotation.z),	-sin(Rotation.z),	0,
			sin(Rotation.z),	cos(Rotation.z),	0,
			0,					0,					1
		);
		//Sum
		inPoint = rz*ry*rx*inPoint;

		//Translate
		inPoint = inPoint + Location;

		return inPoint;
	}

	Vec3<float> TransformPointWorldToLocal(Vec3<float> inPoint)
	{
		//Translate
		inPoint = inPoint - Location;

		//Rotate
		//x
		Mat3 rx(
			1, 0, 0,
			0, cos(Rotation.x), sin(Rotation.x),
			0, -sin(Rotation.x), cos(Rotation.x)
		);
		//y
		Mat3 ry(
			cos(Rotation.y), 0, -sin(Rotation.y),
			0, 1, 0,
			sin(Rotation.y), 0, cos(Rotation.y)
		);
		//z
		Mat3 rz(
			cos(Rotation.z), sin(Rotation.z), 0,
			-sin(Rotation.z), cos(Rotation.z), 0,
			0, 0, 1
		);
		//Sum
		inPoint = rx*ry*rz*inPoint;


		//Scale
		inPoint.x /= Scale.x;
		inPoint.y /= Scale.y;
		inPoint.z /= Scale.z;

		return inPoint;
	}
};

class Line2d {

	Vec2<float> p0, p1;
	float k;

	void CalcK() {

		Vec2<float> pt = p0 - p1;
		if ((pt.x == FInf() && pt.y == FInf())||
			(pt.x == FNInf() && pt.y == FNInf())){
			k = 1;
		}
		else if ((pt.x == FNInf() && pt.y == FInf()) ||
			(pt.x == FInf() && pt.y == FNInf())) {
			k = -1;
		}
		else if (fabs(pt.x) <= 1e-6) {
			k = FInf();//Infinite
		}
		else if (fabs(pt.y) <= 1e-6) {
			k = 0;
		}
		else {
			k = pt.y / pt.x;
		}
	}

public:
	Line2d(Vec2<float> &p0, Vec2<float> &p1){
		this->p0 = p0, this->p1 = p1;
		CalcK();
	}

	Line2d(Line2d &&l) {
		this->p0 = l.GetP0(), this->p1 = l.GetP1(), this->k = l.GetK();
	}

	Vec2<float> GetP0() { return p0; }

	void SetP0(Vec2<float> p) { p0 = p; CalcK(); }

	Vec2<float> GetP1() {return p1;}

	void SetP1(Vec2<float> p) { p1 = p; CalcK(); }

	float GetK() { return k; }

	float GetXfromY(float y, bool inSection = false) {

		if ((y < fmin(p0.y, p1.y) || y > fmax(p0.y, p1.y)) && inSection ||
			fabs(k) <= 1e-6) {
			return FNan();
		}

		if (k == FInf()) {
			return p0.x;
		}

		
		if (fabs(p0.y) != FInf()) {
			return p0.x + (y - p0.y) * (1 / k);
		}
		else if (fabs(p1.y) != FInf()) {
			return p1.x + (y - p1.y) * k;
		}
		else return FNan();
	}

	float GetYfromX(float x, bool inSection = false) {
		if ((x < fmin(p0.x, p1.x) || x > fmax(p0.x, p1.x)) && inSection ||
			k == FInf()) {
			return FNan();
		}

		if (fabs(k) <= 1e-6) {
			return p0.y;
		}

		if (fabs(p0.x) != FInf()) {
			return p0.y + (x - p0.x) * k;
		}
		else if (fabs(p1.x) != FInf()) {
			return p1.y + (x - p1.x) * k;
		}
		else return FNan();
	}

private:
	// Cohen–Sutherland algorithm
	// Compute the bit code for a point (x, y) using the clip rectangle
	// bounded diagonally by (xmin, ymin), and (xmax, ymax)
	//	INSIDE = 0; // 0000
	// 	LEFT = 1;   // 0001
	// 	RIGHT = 2;  // 0010
	// 	BOTTOM = 4; // 0100
	// 	TOP = 8;    // 1000

	int ComputeOutcode(Vec2<float> p, float xmin, float xmax, float ymin, float ymax) {
		int code;
		code = 0;          // initialised as being inside of clip window

		if (p.x < xmin)           // to the left of clip window
			code |= 1;
		else if (p.x > xmax)      // to the right of clip window
			code |= 2;
		if (p.y < ymin)           // below the clip window
			code |= 4;
		else if (p.y > ymax)      // above the clip window
			code |= 8;

		return code;
	}

public:
	// Cohen–Sutherland clipping algorithm clips a line from
	// P0 = (x0, y0) to P1 = (x1, y1) against a rectangle with 
	// diagonal from (xmin, ymin) to (xmax, ymax).
	// return true if line is visible
	bool ClipByRect(Vec2<float> rectBL, Vec2<float> rectTR) {

		float xmin = fmin(rectBL.x, rectTR.x);
		float xmax = fmax(rectBL.x, rectTR.x);
		float ymin = fmin(rectBL.y, rectTR.y);
		float ymax = fmax(rectBL.y, rectTR.y);
		
		if (p0.x == FInf()) { p0.x = xmax; } 
		else if (p0.x == FNInf()) { p0.x = xmin; }

		if (p0.y == FInf()) { p0.y = ymax; }
		else if (p0.y == FNInf()) { p0.y = ymin; }

		if (p1.x == FInf()) { p1.x = xmax; }
		else if (p1.x == FNInf()) { p1.x = xmin; }

		if (p1.y == FInf()) { p1.y = ymax; }
		else if (p1.y == FNInf()) { p1.y = ymin; }

		CalcK();

		// compute outcodes for P0, P1, and whatever point lies outside the clip rectangle
		int outcode0 = ComputeOutcode(p0, xmin, xmax, ymin, ymax);
		int outcode1 = ComputeOutcode(p1, xmin, xmax, ymin, ymax);
		bool accept = false;

		while (true) {
			if (!(outcode0 | outcode1)) { //相或为0，接受并且退出循环
				accept = true;
				
				break;
			}
			else if (outcode0 & outcode1) { // 相与为1，拒绝且退出循环
				break;
			}
			else {
				// failed both tests, so calculate the line segment to clip
					// from an outside point to an intersection with clip edge
				float x, y;

				//找出在界外的点
				int outcodeOut = outcode0 ? outcode0 : outcode1;

				// 找出和边界相交的点
				// 使用点斜式 y = y0 + slope * (x - x0), x = x0 + (1 / slope) * (y - y0)
				//	INSIDE = 0; // 0000
				// 	LEFT = 1;   // 0001
				// 	RIGHT = 2;  // 0010
				// 	BOTTOM = 4; // 0100
				// 	TOP = 8;    // 1000
				if (outcodeOut & 8) {           // point is above the clip rectangle
					//x = p0.x + (p1.x - p0.x) * (ymax - p0.y) / (p1.y - p0.y);
					x = GetXfromY(ymax);
					y = ymax;
				}
				else if (outcodeOut & 4) { // point is below the clip rectangle
					//x = p0.x + (p1.x - p0.x) * (ymin - p0.y) / (p1.y - p0.y);
					x = GetXfromY(ymin);
					y = ymin;
				}
				else if (outcodeOut & 2) {  // point is to the right of clip rectangle
					//y = p0.y + (p1.y - p0.y) * (xmax - p0.x) / (p1.x - p0.x);
					y = GetYfromX(xmax);
					x = xmax;
				}
				else if (outcodeOut & 1) {   // point is to the left of clip rectangle
					//y = p0.y + (p1.y - p0.y) * (xmin - p0.x) / (p1.x - p0.x);
					y = GetYfromX(xmin);
					x = xmin;
				}

				// Now we move outside point to intersection point to clip
				// 为什么继续循环，两个端点都有可能在外面
				if (outcodeOut == outcode0) {
					p0.x = x;
					p0.y = y;
					CalcK();
					outcode0 = ComputeOutcode(p0, xmin, xmax, ymin, ymax);
				}
				else {
					p1.x = x;
					p1.y = y;
					CalcK();
					outcode1 = ComputeOutcode(p1, xmin, xmax, ymin, ymax);
				}
			}
		}
		return accept;
	}
};


class Triangle2d
{
public:
	//Vertices
	Vec2<float> v0, v1, v2, v3;

	//UVs
	Vec2<float> uv0, uv1, uv2, uv3;



	Triangle2d() {}

	Triangle2d(	Vec2<float> &v0, Vec2<float> &v1, Vec2<float> &v2,
		Vec2<float> &uv0, Vec2<float> &uv1, Vec2<float> &uv2) {

		this->v0 = v0; this->v1 = v1; this->v2 = v2;
		this->uv0 = uv0; this->uv1 = uv1; this->uv2 = uv2;
	}
	Vec2<float> V0() { return v0; }
	Vec2<float> V1() { return v1; }
	Vec2<float> V2() { return v2; }
	Vec2<float> UV0() { return uv0; }
	Vec2<float> UV1() { return uv1; }
	Vec2<float> UV2() { return uv2; }
	
	Vec2<float> GetUV(Vec2<float> p) {
		//           v0   s: total space of triangle
		//         / /|   a: space of triangle v1-p-v2
		//       c  / |   b: space of triangle v0-p-v2
		//     /   /  |   c: space of triangle v0-p-v1
		//   v1---p   b
		//     \   \  |
		//       a  \ |
		//         \ \|
		//           v2
		Vec2<float> p0(v0 - p), p1(v1 - p), p2(v2 - p);
		float s = ((v1 - v0).CrossProduct(v2 - v0));
		float a = ((p1).CrossProduct(p2));
		float b = ((p2).CrossProduct(p0));
		float c = ((p0).CrossProduct(p1));

		return uv0 * (a / s) + uv1 * (b / s) + uv2 * (c / s);
	}

};

// Normal direction and vertex order
// Normal vec = crossproduct( v1-v0, v2-v0 )
//           v0
//         /  |
//       /    |
//     /      |
//   v1   UP  |
//     \      |
//       \    |
//         \  |
//           v2
class Triangle3d{


	//Vertices
	Vec3<float> v0, v1, v2;

	Vec3<float> normal;

	Vec3<float> plane;
	//UVs
	Vec2<float> uv0, uv1, uv2;

	//Projection parameters
	Triangle2d t_projected;
	//u = a1x + b1y + c1
	//v = a2x + b2y + c2
	Vec3<float> uparam, vparam;
	//plane that project 3d point to
	enum {
		plane_xy,
		plane_xz,
		plane_yz
	} projPlane;

	void CalcNorm() {
		Vec3<float> vec0(v1-v0), vec1(v2-v0);
		normal = (v1 - v0).CrossProduct(v2 - v0);
		//normal.Normalize();
	}
public:
	void UpdateProjectionParams() {
		CalcNorm();
		float d = normal * v0;
		plane = normal / d;

		//Calculate plane according to normal vector
		if (fabs(normal.y) > fabs(normal.x)) {
			if (fabs(normal.z) > fabs(normal.y)) {
				projPlane = plane_xy;
			}
			else {
				projPlane = plane_xz;
			}
		}
		else {
			projPlane = plane_yz;
		}

		//project

		switch (projPlane)
		{
		case plane_xy: {
			t_projected.v0.x = v0.x, t_projected.v0.y = v0.y;
			t_projected.v1.x = v1.x, t_projected.v1.y = v1.y;
			t_projected.v2.x = v2.x, t_projected.v2.y = v2.y;
		}break;
		case plane_xz: {
			t_projected.v0.x = v0.x, t_projected.v0.y = v0.z;
			t_projected.v1.x = v1.x, t_projected.v1.y = v1.z;
			t_projected.v2.x = v2.x, t_projected.v2.y = v2.z;
		}break;
		case plane_yz: {
			t_projected.v0.x = v0.y, t_projected.v0.y = v0.z;
			t_projected.v1.x = v1.y, t_projected.v1.y = v1.z;
			t_projected.v2.x = v2.y, t_projected.v2.y = v2.z;
		}break;
		}
		Vec2<float> duv1(uv0 - uv1), duv2(uv0-uv2);
		Vec2<float> dxy1(t_projected.v0 - t_projected.v1);
		Vec2<float> dxy2(t_projected.v0 - t_projected.v2);
		float div = dxy2.x * dxy1.y - dxy1.x * dxy2.y;
		uparam.x = (duv2.x*dxy1.y - duv1.x*dxy2.y) / div;
		vparam.x = (duv2.y*dxy1.y - duv1.y*dxy2.y) / div;
		uparam.y = (duv1.x*dxy2.x - duv2.x*dxy1.x) / div;
		vparam.y = (duv1.y*dxy2.x - duv2.y*dxy1.x) / div;
		uparam.z = uv0.x - uparam.x * t_projected.v0.x - uparam.y * t_projected.v0.y;
		vparam.z = uv0.y - vparam.x * t_projected.v0.x - vparam.y*t_projected.v0.y;
	}

	Triangle3d(Vec3<float> &v0, Vec3<float> &v1, Vec3<float> &v2) {
		this->v0 = v0; this->v1 = v1; this->v2 = v2;
		CalcNorm();
		float d = normal * v0;
		plane = normal / d;
	}
	Triangle3d(	Vec3<float> &v0, Vec3<float> &v1, Vec3<float> &v2, 
				Vec2<float> &uv0, Vec2<float> &uv1, Vec2<float> &uv2) {
		this->v0 = v0; this->v1 = v1; this->v2 = v2;
		this->uv0 = uv0; this->uv1 = uv1; this->uv2 = uv2;
		CalcNorm();
		float d = normal * v0;
		plane = normal / d;
	}
	Vec3<float> V0() { return v0; }
	Vec3<float> V1() { return v1; }
	Vec3<float> V2() { return v2; }
	Vec2<float> UV0() { return uv0; }
	Vec2<float> UV1() { return uv1; }
	Vec2<float> UV2() { return uv2; }
	Vec3<float> Normal() { return normal; }

	//Ax + By + Cz = 1
	//return Vec3{A, B, C}
	Vec3<float> GetPlaneParams() { return plane; }

	void TransformToWorld(Transform &t) {
		v0 = t.TransformPointLocalToWorld(v0);
		v1 = t.TransformPointLocalToWorld(v1);
		v2 = t.TransformPointLocalToWorld(v2);
	}

	void TransformToLocal(Transform &t) {
		v0 = t.TransformPointWorldToLocal(v0);
		v1 = t.TransformPointWorldToLocal(v1);
		v2 = t.TransformPointWorldToLocal(v2);
	}

	//Not really in z axis, only means depth
	//depth is reperesented in y axis in 3d coords
	float OneOverZ(Vec2<float> p, float d) {
		return (plane.x / d)*p.x + (plane.z / d)*p.y + plane.y;
	}

	Vec3<float> PerspProject(Vec2<float> p, float d) {
		Vec3<float> p3d;
		p3d.y = 1 / OneOverZ(p, d);
		p3d.x = p.x*p3d.y / d;
		p3d.z = p.y*p3d.y / d;
		return p3d;
	}

	Vec2<float> GetUV(Vec3<float> p3d) {
/*
		//Project triangle to plane which it has maximum surface area on
		enum {
			plane_xy,
			plane_xz,
			plane_yz
		} projPlane;

		//Calculate plane according to normal vector
		if (fabs(normal.y) > fabs(normal.x)) {
			if (fabs(normal.z) > fabs(normal.y)) {
				projPlane = plane_xy;
			}
			else {
				projPlane = plane_xz;
			}
		}
		else {
			projPlane = plane_yz;
		}

		//project
		Triangle2d t; Vec2<float> p;
		t.uv0 = uv0, t.uv1 = uv1, t.uv2 = uv2;
		switch (projPlane)
		{
		case plane_xy: {
			t.v0.x = v0.x, t.v0.y = v0.y;
			t.v1.x = v1.x, t.v1.y = v1.y;
			t.v2.x = v2.x, t.v2.y = v2.y;
			p.x = p3d.x, p.y = p3d.y;
		}break;
		case plane_xz: {
			t.v0.x = v0.x, t.v0.y = v0.z;
			t.v1.x = v1.x, t.v1.y = v1.z;
			t.v2.x = v2.x, t.v2.y = v2.z;
			p.x = p3d.x, p.y = p3d.z;
		}break;
		case plane_yz: {
			t.v0.x = v0.y, t.v0.y = v0.z;
			t.v1.x = v1.y, t.v1.y = v1.z;
			t.v2.x = v2.y, t.v2.y = v2.z;
			p.x = p3d.y, p.y = p3d.z;
		}break;
		}

		return t.GetUV(p);*/
		switch (projPlane)
		{
		case plane_xy: {
			//t_projected.v0.x = v0.x, t_projected.v0.y = v0.y;
			//t_projected.v1.x = v1.x, t_projected.v1.y = v1.y;
			//t_projected.v2.x = v2.x, t_projected.v2.y = v2.y;
			return Vec2<float>(	uparam.x*p3d.x + uparam.y*p3d.y + uparam.z,
								vparam.x*p3d.x + vparam.y*p3d.y + vparam.z);
		}break;
		case plane_xz: {
			//t_projected.v0.x = v0.x, t_projected.v0.y = v0.z;
			//t_projected.v1.x = v1.x, t_projected.v1.y = v1.z;
			//t_projected.v2.x = v2.x, t_projected.v2.y = v2.z;
			return Vec2<float>(	uparam.x*p3d.x + uparam.y*p3d.z + uparam.z,
								vparam.x*p3d.x + vparam.y*p3d.z + vparam.z);
		}break;
		case plane_yz: {
			//t_projected.v0.x = v0.y, t_projected.v0.y = v0.z;
			//t_projected.v1.x = v1.y, t_projected.v1.y = v1.z;
			//t_projected.v2.x = v2.y, t_projected.v2.y = v2.z;
			return Vec2<float>(	uparam.x*p3d.y + uparam.y*p3d.z + uparam.z,
								vparam.x*p3d.y + vparam.y*p3d.z + vparam.z);
		}break;
		}
		
	}
};



class Image2d
{
	uint32_t xRes, yRes;
	float *pixelData = nullptr;//[RGBA]
public:
	//Ctor
	Image2d(int xRes, int yRes)
	{
		this->xRes = xRes;
		this->yRes = yRes;

		//allocate pixel buffer
		pixelData = new float[xRes * yRes * 4];

		for (int xPos = 0; xPos < xRes; xPos++)
		{
			for (int yPos = 0; yPos < yRes; yPos++)
			{
				pixelData[(xPos * yRes + yPos) * 4] = 1;
				pixelData[(xPos * yRes + yPos) * 4 + 1] = 1;
				pixelData[(xPos * yRes + yPos) * 4 + 2] = 1;
				pixelData[(xPos * yRes + yPos) * 4 + 3] = 1;
			}
		}
	}

	Image2d(int xRes, int yRes, Color fillColor)
	{
		this->xRes = xRes;
		this->yRes = yRes;

		//allocate pixel buffer
		pixelData = new float[xRes * yRes * 4];

		for (int xPos = 0; xPos < xRes; xPos++)
		{
			for (int yPos = 0; yPos < yRes; yPos++)
			{
				pixelData[(xPos * yRes + yPos) * 4] = fillColor.r;
				pixelData[(xPos * yRes + yPos) * 4 + 1] = fillColor.g;
				pixelData[(xPos * yRes + yPos) * 4 + 2] = fillColor.b;
				pixelData[(xPos * yRes + yPos) * 4 + 3] = fillColor.a;
			}
		}
	}

	~Image2d() {
		delete[] pixelData;
	}

	Vec2<uint32_t> GetImageRes() { return Vec2<uint32_t>(xRes, yRes); }

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

	Color GetPixelFromUV(float u, float v)
	{
		if (pixelData == nullptr)
			return Color();

		return GetPixel((u - floor(u))*xRes, (v - floor(v))*yRes);
	}

	void SetPixel(int xPos, int yPos, Color color)
	{
		if (pixelData == nullptr)
			return;
		if (xPos < 0 || xPos >= xRes || yPos < 0 || yPos >= yRes)
			return;

		pixelData[(xPos * yRes + yPos) * 4] = color.r;
		pixelData[(xPos * yRes + yPos) * 4 + 1] = color.g;
		pixelData[(xPos * yRes + yPos) * 4 + 2] = color.b;
		pixelData[(xPos * yRes + yPos) * 4 + 3] = color.a;
	}

	//Operators
	float* operator&() { return pixelData; }
};