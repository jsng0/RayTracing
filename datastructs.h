// Name: Jordon Ng 860977164
// Quarter, Year: Fall 2014
// Project: 2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#ifndef DATA_STRUCTS
#define DATA_STRUCTS

#include <cmath>
#include <algorithm> //for min and max

// This is an overly simplified container
// for an RGB Color
struct Color3d
{
	double r;
	double g;
	double b;
	double a;

	Color3d()
		: r(0.0), g(0.0), b(0.0), a(0.0)
	{}
	Color3d(double r, double g, double b)
		: r(r), g(g), b(b), a(0.0)
	{}

	// Changes color based on the rotation by
   // the color wheel, ...
	// using steps RGB -> HSL -> transform -> RGB
	void rotateHue(double degrees)
	{
		// hue, saturation, and luminosity
		double h = 0.0;
		double s = 0.0;
		double l = 0.0;
		double cmin = std::min(r, std::min(g, b));
		double cmax = std::max(r, std::max(g, b));

		// Convert to HSV
		l = (cmax - cmin) / 2.0;
		if (cmax - cmin != 0.0)
		{
			if (cmax == r)
			{
				h = (g - b) / (cmax - cmin);
				h = h < 0.0 ? h + 6.0 : h;
			}
			else if (cmax == g)
				h = (b - r) / (cmax - cmin) + 2.0;
			else
				h = (r - g) / (cmax - cmin) + 4.0;
			h *= 60;
			s = (cmax - cmin) / (1 - abs(2 * l - 1));
		}

		// Rotate Hue
		h += degrees;
		h = h > 360.0 ? h - 360.0 : h;
		h = h < 0.0 ? h + 360.0 : h;

		// Convert back to RGB and reassign
		double c = 1 - std::abs(2.0 * l - 1.0) * s;
		double x = c * (1 - std::abs(std::fmod(h / 60.0, 2.0) - 1.0));
		double m = l - c / 2.0;;

		r = 0.0;
		g = 0.0;
		b = 0.0;
		if (h < 60.0) {
			r = c;
			g = x;
		} else if (h < 120.0) {
			r = x;
			g = c;
		} else if (h < 180.0) {
			g = c;
			b = x;
		} else if (h < 240.0) {
			g = x;
			b = c;
		} else if (h < 300.0) {
			b = c;
			r = x;
		} else {
			b = x;
			r = c;
		}
		r += m;
		g += m;
		b += m;
	}
};

struct my_vertices //holds vertices to connect each triangle
{
	int v1;
	int v2;
	int v3;
	my_vertices(int x, int y, int z)
		:v1(x), v2(y), v3(z) {}
};

// A simple wrapper for store 3D vectors
struct Vector3
{
	float x;
	float y;
	float z;

	Vector3()
		:x(0.0), y(0.0), z(0.0) {}

	Vector3(float x, float y, float z)
		:x(x), y(y), z(z) {}

	Vector3(const Vector3 & v)
		:x(v.x), y(v.y), z(v.z) {}

	Vector3 operator+(const Vector3 & rhs) const
	{ return Vector3(x + rhs.x, y + rhs.y, z + rhs.z); }
	Vector3 operator-(const Vector3 & rhs) const
	{ return Vector3(x - rhs.x, y - rhs.y, z - rhs.z); }
	Vector3 operator*(const Vector3 & rhs) const
	{ return Vector3(x * rhs.x, y * rhs.y, z * rhs.z); }
	Vector3 operator*(float rhs) const
	{ return Vector3(x * rhs, y * rhs, z * rhs); }
	Vector3 operator/(float rhs) const
	{ return Vector3(x / rhs, y / rhs, z / rhs); }
	Vector3 operator+=(const Vector3 & rhs)
	{ x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	Vector3 operator-=(const Vector3 & rhs)
	{ x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
	Vector3 operator*=(float rhs)
	{ x *= rhs; y *= rhs; z *= rhs; return *this; }
	Vector3 operator/=(float rhs)
	{ x /= rhs; y /= rhs; z /= rhs; return *this; }

	float magnitude() const
	{ return sqrt(x * x + y * y + z * z); }

	void normalize() { *this /= magnitude(); }

	Vector3 normalized() const	{ return *this / magnitude(); } //returns L_2 norm

	float dot(const Vector3 & rhs) const
	{
		return x * rhs.x + y * rhs.y + z * rhs.z;
	}

	Vector3 cross(const Vector3 & rhs) const
	{
		return Vector3(y * rhs.z - z * rhs.y,
					   z * rhs.x - x * rhs.z,
		               x * rhs.y - y * rhs.x);
	}

	float vec_max(const Vector3 & rhs) const
	{
		return (this->magnitude() > rhs.magnitude()) ? this->magnitude() : rhs.magnitude();
	}
};

struct Triangle
{
	Vector3 v1;
	Vector3 v2;
	Vector3 v3;
	Vector3 v_norm;
	Color3d rgb_color;

	//Vector3 material_ambient;
	//Vector3 material_diffuse;


	float shininess_coeff;
	Triangle(Vector3 one, Vector3 two, Vector3 three, Color3d c)
	:v1(one), v2(two), v3(three), rgb_color(c), shininess_coeff(100)
	{
		Vector3 ba = two-one;
		Vector3 ca = three-one;
		v_norm = ba.cross(ca);
		//v_norm.normalize();
	}
};

struct Ray
{
	Vector3 origin;
	Vector3 direction;
	Ray() {Vector3 d(0.0, 0.0, 1.0);  direction = d;}
	Ray(const Vector3& o, const Vector3& dir)
	{
		origin = o;
		Vector3 d(0.0, 0.0, 1.0);
		float mag = dir.magnitude();
		if (mag > 0.0) {d = dir;}
		direction = d;
	}
};

struct Lighting
{
	Vector3 diffuse;
	Vector3 specular;
	Vector3 total_illu;

	Lighting()
		:diffuse(Vector3(0.0,0.0,0.0)), specular(Vector3(0.0,0.0,0.0)) {}

	Lighting(Vector3 d, Vector3 s)
		:diffuse(d), specular(s) {}
};

struct Light
{
	Vector3 pos;
	Vector3 diffuse_col;
	float diffuse_power;
	Vector3 specular_col;
	float specular_power;

	Light(Vector3 p, Vector3 dc, float dp, Vector3 sc, float sp)
		:pos(p),diffuse_col(dc),diffuse_power(dp),specular_col(sc),specular_power(sp) {}
};

#endif
