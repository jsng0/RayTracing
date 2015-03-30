// Name: Jordon Ng 860977164
// Quarter, Year: Fall 2014
// Project: 2
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#include <cmath> //for sin and cos
#include <vector>
#include <GL/glut.h>
#include <algorithm>
#include <iostream>
#include <fstream> //for ifstream
#include "datastructs.h"


using namespace std;

//--------------------------------GLOBAL DATA-----------------------------------
const int WINDOW_WIDTH = 1000;
const int WINDOW_HEIGHT = 1000;
const float VIEW_LEFT = 0.0;
const float VIEW_RIGHT = WINDOW_WIDTH;
const float VIEW_BOTTOM = 0.0;
const float VIEW_TOP = WINDOW_HEIGHT;
const float VIEW_NEAR = -400;
const float VIEW_FAR = 400;

Vector3 viewer(400,400,0);
Vector3 light_loc(400,800,0);

Vector3 light_diffuse(0.2,0.2,0.2);
Vector3 light_specular(0.2,0.2,0.2);
Vector3 light_ambient(0.2,0.2,0.2);

Vector3 material_spec(0.2,0.2,0.2);
Vector3 material_diff(0.4,0.2,0.6);
Vector3 material_ambient(0.4,0.2,0.6);
Vector3 surface_col(0.0,0.0,1.0);

float diffuse_power = 0.1;
float specular_power = 25;
Light L(light_loc,light_diffuse,diffuse_power,light_specular,specular_power);

bool render_toggle = false;
//--------------------------------GLOBAL DATA-----------------------------------


//--------------------------------GLOBAL DATA CONTAINERS------------------------
vector<Vector3> poly_pts; //holds the points
vector<my_vertices> poly_vertices; //holds which 3 points to connect
vector<Triangle> triangles; //holds number of triangles each poly has
vector<Vector3> poly_normals; //holds normals of the vertices of each triangle
//--------------------------------GLOBAL DATA CONTAINERS------------------------

// Renders a quad at cell (x, y) with dimensions CELL_LENGTH
void renderPixel(const int& x, const int& y, Color3d& color, float sz = 1.0)
{
	glPointSize(sz);
	glColor4f(color.r, color.g, color.b, color.a);
	glBegin(GL_POINTS);
	glVertex2i(x, y);
	glEnd();
}

/* @param: references to corresponding 3x3 matrix indices
 * @return: the cross product of the 3x3 matrix
*/
float get_det(const float& topx, const float& topy, const float& topz,
			  const float& midx, const float& midy, const float& midz,
			  const float& botx, const float& boty, const float& botz)
{
	return topx * (midy*botz - boty*midz)
		  -topy * (midx*botz - botx*midz)
		  +topz * (midx*boty - botx*midy);
}

/* @param: current ray, which triangle, and pointer to t value
 * @return: whether ray intersects the given triangle using barycentric coord
*/
bool rayIntersectsTriangle(const Ray& ray, const Triangle& tri, float* t = NULL)
{
	if(t == NULL) return false;

	Vector3 a = tri.v1; Vector3 b = tri.v2; Vector3 c = tri.v3;
	Vector3 d = ray.direction; Vector3 e = ray.origin;
	float abx = a.x-b.x; float aby = a.y-b.y; float abz = a.z-b.z;
	float acx = a.x-c.x; float acy = a.y-c.y; float acz = a.z-c.z;
	float aex = a.x-e.x; float aey = a.y-e.y; float aez = a.z-e.z;

	float det_a = get_det(abx, acx, d.x,
						  aby, acy, d.y,
						  abz, acz, d.z);

	float beta = get_det(aex, acx, d.x,
	                     aey, acy, d.y,
						 aez, acz, d.z) / det_a;

	float gamma = get_det(abx, aex, d.x,
	                      aby, aey, d.y,
						  abz, aez, d.z) / det_a;

	float nt = get_det(abx, acx, aex,
	                   aby, acy, aey,
				       abz, acz, aez) / det_a;

	if(nt > 0)
	{
		if(beta >= 0 && gamma >= 0)
		{
			if(beta + gamma < 1)
			{
				*t = nt; //update t
				return true;
			}
		}
	}
	return false;
}

/* @param: none
 * @return: sets all triangles in polygon back to og colors
*/
void reset_colors()
{
	for(int i = 0; i < triangles.size(); ++i)
	{
		triangles[i].rgb_color = Color3d(surface_col.x,surface_col.y,surface_col.z);
	}
}

/* @param: ray to check intersection, number of triangles in the polygon
 * @return: index of nearest triangle hit, -1 if no triangle hit
*/
int Compute_Closest_Intersect(const Ray& S, const int& num_triangles)
{
	float t_min = -1;
	float min = 1000;
	int min_tri_index = -1;
	for (int i = 0; i < num_triangles; ++i)
	{
		if(rayIntersectsTriangle(S, triangles[i], &t_min))
		{
			if(t_min <  min)
			{
				min = t_min;
				if(min_tri_index == -1) min_tri_index = i;
				min_tri_index = (i < min_tri_index) ? i : min_tri_index; //replace min_i if i is lower
			}
		}
	}
	return min_tri_index;
}

/* @param: vector of intersection, p, and min tri that intersected the triangle
 * @return: interpolated normal
*/
Vector3 Get_Interpolate(const Vector3& intersect, const int& tri_index)
{
	float d0 = (triangles[tri_index].v1 - intersect).magnitude(); //get distance of each side
	float d1 = (triangles[tri_index].v2 - intersect).magnitude();
	float d2 = (triangles[tri_index].v3 - intersect).magnitude();

	Vector3 a, b, c, Na, Nb, Nc;

	if(d0 > d1 && d1 > d2) //reorder vertices so a is furthest away
	{
		a = triangles[tri_index].v1;
		b = triangles[tri_index].v2;
		c = triangles[tri_index].v3;
		Na = poly_normals[poly_vertices[tri_index].v1];
		Nb = poly_normals[poly_vertices[tri_index].v2];
		Nc = poly_normals[poly_vertices[tri_index].v3];
	}
	else if(d1 > d0 && d1 > d2)
	{
		a = triangles[tri_index].v2;
		b = triangles[tri_index].v1;
		c = triangles[tri_index].v3;
		Na = poly_normals[poly_vertices[tri_index].v2];
		Nb = poly_normals[poly_vertices[tri_index].v1];
		Nc = poly_normals[poly_vertices[tri_index].v3];
	}
	else
	{
		a = triangles[tri_index].v3;
		b = triangles[tri_index].v1;
		c = triangles[tri_index].v2;
		Na = poly_normals[poly_vertices[tri_index].v3];
		Nb = poly_normals[poly_vertices[tri_index].v1];
		Nc = poly_normals[poly_vertices[tri_index].v2];
	}
	//compute new normal
	Vector3 edge_origin(c);
	Vector3 edge_normal((c-b).cross(triangles[tri_index].v_norm));

	Ray R(a, (intersect - a).normalized());

	float t = (edge_normal.dot(edge_origin-R.origin)) / edge_normal.dot(R.direction);
	Vector3 q = a + R.direction * t; //R_o + R_d * t from lecture3 slide 11

	float distanceB_Q = abs(Vector3(b-q).magnitude());
	float distanceB_C = abs(Vector3(b-c).magnitude());
	Vector3 Nq = Nb + (Nc - Nb) * (distanceB_Q / distanceB_C);

	float distanceQ_I = abs(Vector3(q - intersect).magnitude());
	float distanceQ_A = abs(Vector3(q - a).magnitude());
	Vector3 Ni = Nq + (Na - Nq) * (distanceQ_I / distanceQ_A);

	return Ni;
}

/* @param: variables from pixel on for ray tracing
 * @return: updates color of the closest triangle based on phong illumination
*/
void smooth_shade(const Ray& eye,const int& tri_index, const int& min, const int& num_tri)
{
	Vector3 intersect(eye.origin.x + eye.direction.x*min,
					  eye.origin.y + eye.direction.y*min,
					  eye.origin.z + eye.direction.z*min);

	Vector3 illum;
	Color3d tri_col = triangles[tri_index].rgb_color;

	Vector3 light_dir = intersect - light_loc;
	Ray shadow_dragon(light_loc, light_dir);

	int shadow_tri = Compute_Closest_Intersect(shadow_dragon,num_tri); //Check shadow

	if(shadow_tri >= 0 && shadow_tri != tri_index)
	{
		illum.x = (tri_col.r*0.2);
		illum.y = (tri_col.g*0.2);
		illum.z = (tri_col.b*0.2);
		triangles[tri_index].rgb_color = Color3d(illum.x, illum.y, illum.z);
		return;
	}

	Vector3 normal = triangles[tri_index].v_norm;

	if(tri_index < triangles.size() - 2) //don't include ground triangles
		normal = Get_Interpolate(intersect,tri_index); //Get new interpolated triangle
	else
		normal = triangles[tri_index].v_norm; //otherwise use the triangle normal

	if(normal.dot(eye.direction) > 0) //check direction of normal
	{
		normal.x = -normal.x;
		normal.y = -normal.y;
		normal.z = -normal.z;
	}
	normal.normalize();

	Vector3 L = (light_loc - intersect).normalized();
	Vector3 V = (viewer - intersect).normalized();

	float NdotL = normal.dot(L);
	Vector3 R_L = (normal * (2 * normal.dot(L))) - L;
	//Vector3 H = (V + L).normalized(); //replace below calc with N.H for blinn phong
	if(NdotL < 0) R_L = Vector3(0.0,0.0,0.0);

	illum.x = (L.dot(normal)*0.3) + (pow(V.dot(R_L),15)*0.3) + (tri_col.r*0.3); //do phong
	illum.y = (L.dot(normal)*0.3) + (pow(V.dot(R_L),15)*0.3) + (tri_col.g*0.3);
	illum.z = (L.dot(normal)*0.3) + (pow(V.dot(R_L),15)*0.3) + (tri_col.b*0.3);

	triangles[tri_index].rgb_color = Color3d(illum.x, illum.y, illum.z);
}

/* @param: variables from pixel on for ray tracing
 * @return: updates color of the closest triangle based on phong illumination
*/
void flat_shade(const Ray& eye,const int& tri_index, const int& min, const int& num_tri)
{
	Vector3 intersect(eye.origin.x + eye.direction.x*min,
					  eye.origin.y + eye.direction.y*min,
					  eye.origin.z + eye.direction.z*min);

	Vector3 illum;
	Color3d tri_col = triangles[tri_index].rgb_color;

	Vector3 light_dir = intersect - light_loc;
	Ray shadow_dragon(light_loc, light_dir);

	int shadow_tri = Compute_Closest_Intersect(shadow_dragon,num_tri); //check if shadow hits

	if(shadow_tri >= 0 && shadow_tri != tri_index) //if shadow is valid, add ambient only
	{
		illum.x = (tri_col.r*0.2);
		illum.y = (tri_col.g*0.2);
		illum.z = (tri_col.b*0.2);
		triangles[tri_index].rgb_color = Color3d(illum.x, illum.y, illum.z);
		return;
	}

	Vector3 normal = triangles[tri_index].v_norm;
	if(normal.dot(eye.direction) > 0)
	{
		normal.x = -normal.x;
		normal.y = -normal.y;
		normal.z = -normal.z;
	}
	normal.normalize();

	Vector3 L = (light_loc - intersect).normalized();
	Vector3 V = (viewer - intersect).normalized();
	float NdotL = normal.dot(L);
	Vector3 R = (normal * (2 * normal.dot(L))) - L;
	//Vector3 H = (V + L).normalized();
	if(NdotL < 0) R = Vector3(0.0,0.0,0.0);

	illum.x = (L.dot(normal)*0.3) + (pow(V.dot(R),15)*0.3) + (tri_col.r*0.2); //do phong illumination for diff,spec,ambient
	illum.y = (L.dot(normal)*0.3) + (pow(V.dot(R),15)*0.3) + (tri_col.g*0.2);
	illum.z = (L.dot(normal)*0.3) + (pow(V.dot(R),15)*0.3) + (tri_col.b*0.2);

	triangles[tri_index].rgb_color = Color3d(illum.x, illum.y, illum.z);
}

//flat shading
int pixelOnFlat(const float& pixel_x, const float& pixel_y)
{
	Ray current(Vector3(pixel_x, pixel_y, 0.0), Vector3(0.0, 0.0, 1.0));
	int num_triangles = triangles.size();
	float t_min = -1;
	float min = 1000;
	int min_tri_index = -1;
	for (int i = 0; i < num_triangles; ++i)
	{
		if(rayIntersectsTriangle(current, triangles[i], &t_min))
		{
			if(t_min <  min)
			{
				min = t_min;
				if(min_tri_index == -1) min_tri_index = i;
				min_tri_index = (i < min_tri_index) ? i : min_tri_index; //replace min_i if i is lower
			}
		}
	}

	//int min_tri_index = Compute_Closest_Intersect(current,num_triangles);
	if(min_tri_index >= 0) //we have a valid index to triangle vector
	{
		reset_colors();
		flat_shade(current,min_tri_index,min,num_triangles);
	}
	return min_tri_index;
}

//Smooth shading
int pixelOnSmooth(const float& pixel_x, const float& pixel_y)
{
	Ray current(Vector3(pixel_x, pixel_y, 0.0), Vector3(0.0, 0.0, 1.0));
	int num_triangles = triangles.size();
	float t_min = -1;
	float min = 1000;
	int min_tri_index = -1;
	for (int i = 0; i < num_triangles; ++i)
	{
		if(rayIntersectsTriangle(current, triangles[i], &t_min))
		{
			if(t_min <  min)
			{
				min = t_min;
				if(min_tri_index == -1) min_tri_index = i;
				min_tri_index = (i < min_tri_index) ? i : min_tri_index; //replace min_i if i is lower
			}

		}
	}

	if(min_tri_index >= 0) //we have a valid index to triangle vector
	{
		reset_colors();
		smooth_shade(current,min_tri_index,min,num_triangles);
	}
	return min_tri_index;
}

void GLKeyboardPress(unsigned char key, int x, int y)
{
	if(key == '0')
	{
		if(!render_toggle)
		{
			cout << "Rendering Smooth Shading\n";
			render_toggle = true;
		}
		else
		{
			cout << "Rendering Flat Shading\n";
			render_toggle = false;
		}
		glutPostRedisplay();
	}
}

void GLrender()
{
	glClear(GL_COLOR_BUFFER_BIT);

	cout << "Starting Rendering...\n";
	for (int i = 0; i < WINDOW_WIDTH; i++)
	{
		for (int j = 0; j < WINDOW_HEIGHT; j++)
		{
			int which_triangle = (render_toggle) ? pixelOnSmooth(i,800-j) : pixelOnFlat(i, 800-j); //if hit, update color
			if(which_triangle >= 0) //if hit, render pixel
			{
				renderPixel(i, 800-j, triangles[which_triangle].rgb_color);
			}
			else //if not hit, set to background color
			{
				Color3d bg_col(1.0,1.0,1.0);
				renderPixel(i, 800-j, bg_col);
			}
		}
	}
	cout << "Finished Rendering - Press '0' To Change Render\n";
	glFlush();
	glutSwapBuffers();
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
	glutInit(argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutCreateWindow("CS 130 - Project 2: Flat & Smooth Shading - Jordon Ng");

	glutKeyboardFunc(GLKeyboardPress);
	glutDisplayFunc(GLrender);

	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
	glOrtho(VIEW_LEFT, VIEW_RIGHT, VIEW_BOTTOM, VIEW_TOP, VIEW_NEAR, VIEW_FAR);
}

/* @param: none (global data structures for points and vertices)
 * @return: populate a vector of normals of each vertex in the polygon
*/
void get_vertex_normals()
{
	int num_vertices = poly_vertices.size();
	int num_points = poly_pts.size();
	for(int i = 0; i < num_points; ++i)
	{
		Vector3 normal_numer;
		float normal_denom;
		for(int j = 0;  j < num_vertices; ++j)
		{
			if(poly_vertices[j].v1 == i || poly_vertices[j].v2 == i || poly_vertices[j].v3 == i)
			{
				normal_numer += triangles[j].v_norm;
				normal_denom += triangles[j].v_norm.magnitude();
			}
		}
		poly_normals.push_back(normal_numer / normal_denom);
		//cout << "tri: " << i  << " " << poly_normals[i].x  << " "<< poly_normals[i].y << " " << poly_normals[i].z << endl;
	}
}

/* @param: nothing (takes global triangle vector)
*  @return: populates a ground plane into the triangles vector
*/
void Make_Ground()
{
	triangles.push_back(Triangle(Vector3(0.0,0.0,0.0), Vector3(1000,0.0,0.0), Vector3(0.0,150,1000), Color3d(0.0,1.0,0.0)));
	triangles.push_back(Triangle(Vector3(1000,150,1000), Vector3(1000,0.0,0.0), Vector3(0.0,150,1000), Color3d(0.0,1.0,0.0)));
}

/* @param: vector of which vertices to make a tri with, and vector of vertices
*  @return: populates the tri vector
*/
void Pts_To_Tri(const vector<my_vertices>& V, const vector<Vector3>& P)
{
	int num_ver = V.size();
	for(unsigned i = 0; i < num_ver; ++i)
	{
		triangles.push_back(Triangle(P[V[i].v1], P[V[i].v2], P[V[i].v3], Color3d(0.0,0.0,1.0)));
	}
}

/* @param: pass in a filename to be read
*  @return: fills in the elements of three vectors: 1. poly_pts 2. poly_vertices 3. zbuffer array
*/
void ProcessFileInput(string fname)
{
	ifstream inbuf;
	inbuf.open(fname.c_str());
	int num_pts, num_tri;
	inbuf >> num_pts >> num_tri; //get input bounds

	for(; num_pts > 0; --num_pts)
	{
		int x_in, y_in, z_in;
		inbuf >> x_in >> y_in >> z_in;
		poly_pts.push_back(Vector3(x_in, y_in, z_in)); //put into vertex vector
	}
	for(; num_tri > 0; --num_tri)
	{
		int x_in, y_in, z_in;
		inbuf >> x_in >> y_in >> z_in;
		poly_vertices.push_back(my_vertices(x_in, y_in, z_in));
	}
	inbuf.close();

	Pts_To_Tri(poly_vertices,poly_pts); //populate the triangle vector
	Make_Ground(); //add triangles to make ground plane
	get_vertex_normals();
}

int main(int argc, char** argv)
{
	if(argc != 2) {cout << "Specify Input File!\n"; return 0;}
	string file_name = argv[1];
	ProcessFileInput(file_name);

	GLInit(&argc, argv);

	glutMainLoop();

	return 0;
}
