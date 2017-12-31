/**
 * minigl.cpp
 * -------------------------------
 * Implement miniGL here.
 *
 * You may include minigl.h and any of the standard C++ libraries.
 * No other includes are permitted.  Other preprocessing directives
 * are also not permitted.  These requirements are strictly
 * enforced.  Be sure to run a test grading to make sure your file
 * passes the sanity tests.
 *
 * The behavior of the routines your are implenting is documented here:
 * https://www.opengl.org/sdk/docs/man2/
 * Note that you will only be implementing a subset of this.  In particular,
 * you only need to implement enough to pass the tests in the suite.
 */

#include "minigl.h"
#include "vec.h"
#include "mat.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <vector>
#include <cstdio>

using namespace std;

/**
 * Useful data types
 */
typedef mat<MGLfloat,4,4> mat4; //data structure storing a 4x4 matrix, see mat.h
typedef mat<MGLfloat,3,3> mat3; //data structure storing a 3x3 matrix, see mat.h
typedef vec<MGLfloat,4> vec4;   //data structure storing a 4 dimensional vector, see vec.h
typedef vec<MGLfloat,3> vec3;   //data structure storing a 3 dimensional vector, see vec.h
typedef vec<MGLfloat,2> vec2;   //data structure storing a 2 dimensional vector, see vec.h

/*--------- GLOBALS ----------*/
MGLpoly_mode polyMode;
MGLmatrix_mode matrixMode;
struct Vertex;
struct Pixel;

vector< Vertex > vertices;
vector< Pixel > pixelVertices;
vector< vector<MGLfloat> > zBuffer;
vector< vec<MGLpixel, 3> > baryVertices;

MGLfloat FAR = 20;
MGLfloat NEAR = -20;

vector< vector< Vertex > > triangles;
vector< vector< Pixel > > pixelTriangles;


vector< mat<MGLfloat,4,4> > projection; 
vector< mat<MGLfloat,4,4> > modelview; 
mat<MGLfloat,4,4> currentModel;
mat<MGLfloat,4,4> currentProj;

MGLpixel RGB[3];

MGLfloat TRI_AREA;

vector< vector<MGLfloat> > colorStack;
int cCount1 = 0;
int cCount2 = 0;
/**
 * Standard macro to report errors
*/

struct Vertex
{
	MGLpixel vColor[3];
	vec<MGLfloat, 4> vertices;
	
	Vertex() 
	{
		vertices[0] = 0;
		vertices[1] = 0;
		vertices[2] = 0;
		vertices[3] = 0;
		
		vColor[0] = RGB[0];
		vColor[1] = RGB[1];
		vColor[2] = RGB[2];
	}
	
	Vertex(MGLfloat x, MGLfloat y, MGLfloat z)
	{
		vertices[0] = x;
		vertices[1] = y;
		vertices[2] = z;
		vertices[3] = 1;
		
		vColor[0] = RGB[0];
		vColor[1] = RGB[1];
		vColor[2] = RGB[2];
	}
	
	const MGLfloat operator[] (int i) const
    {return vertices[i];}
	
	MGLfloat operator[] (int i)
    {return vertices[i];}
	

};


struct Pixel
{
	MGLfloat pColor[3];
	vec<MGLpixel, 2> pixels;
	
	Pixel() {}
	
	Pixel(MGLfloat x, MGLfloat y)
	{
		pixels[0] = static_cast<MGLpixel>(x);
		pixels[1] = static_cast<MGLpixel>(y);
	}
	
	const MGLpixel operator[] (int i) const
    {return pixels[i];}
	
	MGLpixel operator[] (int i)
    {return pixels[i];}
	

};



inline void MGL_ERROR(const char* description) {
    printf("%s\n", description);
    exit(1);
}

MGLfloat getArea(Pixel a, Pixel b, Pixel c)
{
	//do something
	MGLfloat area = 0.5*(a[0]*b[1] - a[1]*b[0] + b[0]*c[1] - b[1]*c[0] + c[0]*a[1] - c[1]*a[0]);
	if(area == TRI_AREA)
	{
		return 0;
	}
	else if(area == 0)
	{
		return TRI_AREA;
	}
	else
	{
		return area;
	}
}

void multiply(Vertex &v, mat<MGLfloat,4,4> m)
{
	int k = 0;
	int ind = 0;
	MGLfloat sum = 0;
	Vertex r;
	
	for(int i = 0; i < 16; i += 4)
	{
		for(int j = 0; j < 4; j++)
		{
			sum += v.vertices[j] * m.values[i+k];
			k++;
		}
		r.vertices[ind] = sum;
		ind++;
		k = 0;
		sum = 0;
	}
	
	for(int i = 0; i < 4; i++)
	{
		v.vertices[i] = r.vertices[i];
	}
		
	return;
}

void multiply(mat<MGLfloat,4,4> ortho, mat<MGLfloat,4,4> &m)
{
	int h = 0;
	int ind = 0;
	MGLfloat sum = 0;
	mat<MGLfloat,4,4> r;
	r.make_zero();

	for(int i = 0; i < 4; i++)
	{
		for(int j = 0; j < 16; j += 4)
		{
			for(int k = 0; k < 4; k++)
			{
				sum += m.values[k+j] * ortho.values[i+h];
				h += 4;
			}
			r.values[j+i] = sum;
			ind++;
			h = 0;
			sum = 0;
			
		}
	}
	m = r;
		
	return;
}





void setColor(MGLfloat alpha, MGLfloat beta, MGLfloat gamma, Pixel &pixel, int t)
{
	MGLfloat w_a = triangles[t][0][3];
	MGLfloat w_b = triangles[t][1][3];
	MGLfloat w_c = triangles[t][2][3];

	MGLfloat realAlpha = (alpha / w_a) / ((alpha / w_a) + (beta / w_b) + (gamma / w_c));
	MGLfloat realBeta = (beta / w_b) / ((alpha / w_a) + (beta / w_b) + (gamma / w_c));
	MGLfloat realGamma = (gamma / w_c) / ((alpha / w_a) + (beta / w_b) + (gamma / w_c));

	MGLfloat red = realAlpha * triangles[t][0].vColor[0] + realBeta*triangles[t][1].vColor[0] + realGamma*triangles[t][2].vColor[0];
	MGLfloat green = realAlpha * triangles[t][0].vColor[1] + realBeta*triangles[t][1].vColor[1] + realGamma*triangles[t][2].vColor[1];
	MGLfloat blue = realAlpha * triangles[t][0].vColor[2] + realBeta*triangles[t][1].vColor[2] + realGamma*triangles[t][2].vColor[2];
	
	
	pixel.pColor[0] = red;
	pixel.pColor[1] = green;
	pixel.pColor[2] = blue;
	
	//cout << "cols: " << red << " " << green << " " << blue << endl;
	
	return;
}

void barycentric(Pixel pixel, MGLfloat &alpha, MGLfloat &beta, MGLfloat &gamma, int i, MGLfloat &zeta)
{
	Pixel a, b, c;
	

	a.pixels[0] = pixelTriangles[i][0][0];
	a.pixels[1] = pixelTriangles[i][0][1];
	
	b.pixels[0] = pixelTriangles[i][1][0];
	b.pixels[1] = pixelTriangles[i][1][1];
	
	c.pixels[0] = pixelTriangles[i][2][0];
	c.pixels[1] = pixelTriangles[i][2][1];
	
	MGLfloat z_a = triangles[i][0][2] / triangles[i][0][3];
	MGLfloat z_b = triangles[i][1][2] / triangles[i][1][3];
	MGLfloat z_c = triangles[i][2][2] / triangles[i][2][3];
	
	
	TRI_AREA = getArea(a, b, c);

	alpha = getArea(pixel, b, c) / TRI_AREA;

	beta = getArea(a, pixel, c) / TRI_AREA;

	gamma = getArea(a, b, pixel) / TRI_AREA;

	
	vec<MGLpixel,3> v(alpha, beta, gamma);
	baryVertices.push_back(v);
	
	zeta = (alpha*z_a + beta*z_b + gamma*z_c);
	
	return;
}

void rasterizeTriangle(MGLsize width, MGLsize height, MGLpixel *data, int triangleNum)
{
	// add pixels between 2 points
	//for(int i = pixelVertices[0][0]; i < pixelVertices[1][0]; i++)
	//{
	//	vec<MGLpixel, 2> v(i,150);
	//	pixelVertices.push_back(v);
	//}
	MGLfloat alpha = 0.0;
	MGLfloat beta = 0.0;
	MGLfloat gamma = 0.0;
	
	MGLpixel x_min = 10000, y_min = 10000, x_max = 0, y_max = 0;
	
	for(int i = 0; i < pixelTriangles[triangleNum].size(); i++)
	{
		if(pixelTriangles[triangleNum][i][0] < x_min)
		{
			x_min = pixelTriangles[triangleNum][i][0];
		}
		if(pixelTriangles[triangleNum][i][0] > x_max)
		{
			x_max = pixelTriangles[triangleNum][i][0];
		}
		if(pixelTriangles[triangleNum][i][1] < y_min)
		{
			y_min = pixelTriangles[triangleNum][i][1];
		}
		if(pixelTriangles[triangleNum][i][1] > y_max)
		{
			y_max = pixelTriangles[triangleNum][i][1];
		}
	}
	
	for(MGLfloat y = y_min; y < y_max; y += 0.5)
	{
		for(MGLfloat x = x_min; x < x_max; x += 0.5)
		{
			MGLfloat zeta = 0;
			Pixel pixel(x, y);
			barycentric(pixel, alpha, beta, gamma, triangleNum, zeta);
			
			if( (alpha >= 0 && alpha <= 1) && (beta >= 0 && beta <= 1) && (gamma >= 0 && gamma <= 1) )
			{
				//it is in the triangle, give it a placeholder color
				pixel.pColor[0] = triangles[triangleNum][0].vColor[0];
				pixel.pColor[1] = triangles[triangleNum][0].vColor[1];
				pixel.pColor[2] = triangles[triangleNum][0].vColor[2];
			
				if(zeta < zBuffer[y][x] && zeta <= FAR && zeta >= NEAR) //may need to switch x and y
				{
					//MGLpixel pColor[3];
					setColor(alpha, beta, gamma, pixel, triangleNum);
					//cout << "colss: " << pixel.pColor[0] << " " << pixel.pColor[1] << " " << pixel.pColor[2] << endl;
					
					*(data + pixel[0] + pixel[1]*width) = Make_Pixel(pixel.pColor[0]*255, pixel.pColor[1]*255, pixel.pColor[2]*255);
					
					//pixelVertices.push_back(p);
					zBuffer[y][x] = zeta;
				}
			}
		}
	}

			// color white
	//MGLfloat R = colorStack.at(cCount2).at(0);
	//MGLfloat G = colorStack.at(cCount2).at(1);
	//MGLfloat B = colorStack.at(cCount2).at(2);

	//for(int i = 0; i < pixelVertices.size(); i++)
	//{


		//*(data + pixelVertices[i][0] + pixelVertices[i][1] * width) = Make_Pixel(R*255,G*255,B*255);
		
	//}
	
	//cCount2++;
}


/**
 * Read pixel data starting with the pixel at coordinates
 * (0, 0), up to (width,  height), into the array
 * pointed to by data.  The boundaries are lower-inclusive,
 * that is, a call with width = height = 1 would just read
 * the pixel at (0, 0).
 *
 * Rasterization and z-buffering should be performed when
 * this function is called, so that the data array is filled
 * with the actual pixel values that should be displayed on
 * the two-dimensional screen.
 */
void mglReadPixels(MGLsize width,
                   MGLsize height,
                   MGLpixel *data)
{
	
	zBuffer.resize(height);
	for(int u = 0; u < height; u++)
	{
		zBuffer.at(u).resize(width, 2);
	}
	//maybe
	


	// convert our 3 vertices into 3 different pixel coordinates
	for(int k = 0; k < triangles.size(); k++)
	{
		for(int t = 0; t < triangles[k].size(); t++)
		{
			//normalize
			triangles[k][t].vertices[0] = triangles[k][t][0] / triangles[k][t][3];
			triangles[k][t].vertices[1] = triangles[k][t][1] / triangles[k][t][3];
			triangles[k][t].vertices[2] = triangles[k][t][2] / triangles[k][t][3];
			
			//convert to pixel coord
			MGLpixel i = (triangles[k][t][0] + 1) * width / 2;
			MGLpixel j = (triangles[k][t][1] + 1) * height / 2;
			
			//create pixel vector
			//vec<MGLpixel, 2> v(i,j);
			Pixel v(i,j);
			
			//push to list
			pixelVertices.push_back(v);
		}
		pixelTriangles.push_back(pixelVertices);
		pixelVertices.clear();
	}
	
	//MGLpixel pixel[2] = { *data, *(data+1) };
	//barycentric(pixel);
	
	
	// change color of our 3 pixels
	for(int i = 0; i < triangles.size(); i++)
	{
		rasterizeTriangle(width, height, data, i);
		
		pixelVertices.clear();
	}
	
	pixelTriangles.clear();
	triangles.clear();
	baryVertices.clear();
	zBuffer.clear();
	vertices.clear();

}

/**
 * Start specifying the vertices for a group of primitives,
 * whose type is specified by the given mode.
 */
void mglBegin(MGLpoly_mode mode)
{
	if(mode == MGL_TRIANGLES)
	{
		polyMode = MGL_TRIANGLES;
	}
	else if(mode == MGL_QUADS)
	{
		polyMode = MGL_QUADS;
	}
	else
	{
		MGL_ERROR("Invalid polygon.");
		exit(1);
	}
}


/**
 * Stop specifying the vertices for a group of primitives.
 */
void mglEnd()
{
	if(polyMode == MGL_TRIANGLES)
	{
		int count = 0;
		vector< Vertex > temp;
		
		for(int i = 0; i < vertices.size(); i++)
		{
			temp.push_back(vertices[i]);
			count++;
			if(count == 3)
			{
				triangles.push_back(temp);
				count = 0;
				temp.clear();
			}
		}
			
		vertices.clear();
	}
	else if(polyMode == MGL_QUADS)
	{
		
		vector< Vertex > A {vertices[0], vertices[1], vertices[2]};
		vector< Vertex > B {vertices[0], vertices[2], vertices[3]};
	
		triangles.push_back(A);
		triangles.push_back(B);
		
		vertices.clear();
	}
}

/**
 * Specify a two-dimensional vertex; the x- and y-coordinates
 * are explicitly specified, while the z-coordinate is assumed
 * to be zero.  Must appear between calls to mglBegin() and
 * mglEnd().
 */
void mglVertex2(MGLfloat x,
                MGLfloat y)
{
	mglVertex3(x, y, 0.00);
}

/**
 * Specify a three-dimensional vertex.  Must appear between
 * calls to mglBegin() and mglEnd().
 */
void mglVertex3(MGLfloat x,
                MGLfloat y,
                MGLfloat z)
{

	Vertex vec(x, y, z);
	
	
	multiply(vec, currentModel);
	multiply(vec, currentProj);
	
	
	vertices.push_back(vec);
	
	
}

/**
 * Set the current matrix mode (modelview or projection).
 */
void mglMatrixMode(MGLmatrix_mode mode)
{
	matrixMode = mode;
	return;
}

/**
 * Push a copy of the current matrix onto the stack for the
 * current matrix mode.
 */
void mglPushMatrix()
{
	if(matrixMode == MGL_PROJECTION)
	{
		cout << "currentProj: " << currentProj << endl;
		projection.push_back(currentProj);
	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		cout << "currentModel: " << currentModel << endl;
		modelview.push_back(currentModel);
	}
}

/**
 * Pop the top matrix from the stack for the current matrix
 * mode.
 */
void mglPopMatrix()
{
	if(matrixMode == MGL_PROJECTION)
	{
		currentProj = projection.back();
		projection.pop_back();
	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		currentModel = modelview.back();
		modelview.pop_back();
	}
}

/**
 * Replace the current matrix with the identity.
 */
void mglLoadIdentity()
{
	if(matrixMode == MGL_PROJECTION)
	{
		currentProj.make_zero();
		for(int i = 0; i < 16; i += 5)
		{
			currentProj.values[i] = 1;
		}
	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		currentModel.make_zero();
		for(int i = 0; i < 16; i += 5)
		{
			currentModel.values[i] = 1;
		}
	}
	else
	{
		MGL_ERROR("Invalid matrix type");
		exit(1);
	}
	return;
}

/**
 * Replace the current matrix with an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglLoadMatrix(const MGLfloat *matrix)
{
	if(matrixMode == MGL_MODELVIEW)
	{
		int j = 0;
		for(int i = 0; i < 4; i++)
		{
			for(int k = 0; k < 16; k += 4)
			{
				currentModel.values[i+k] = *(matrix + j);
				j++;
			}
		}
	}
	else if(matrixMode == MGL_PROJECTION)
	{
		int j = 0;
		for(int i = 0; i < 4; i++)
		{
			for(int k = 0; k < 16; k += 4)
			{
				currentProj.values[i+k] = *(matrix + j);
				j++;
			}
		}
	}
	else
	{
		MGL_ERROR("error");
		exit(1);
	}
}

/**
 * Multiply the current matrix by an arbitrary 4x4 matrix,
 * specified in column-major order.  That is, the matrix
 * is stored as:
 *
 *   ( a0  a4  a8  a12 )
 *   ( a1  a5  a9  a13 )
 *   ( a2  a6  a10 a14 )
 *   ( a3  a7  a11 a15 )
 *
 * where ai is the i'th entry of the array.
 */
void mglMultMatrix(const MGLfloat *matrix)
{
	mat4 dummy;
	dummy.make_zero();
	
	dummy.values[0] = *(matrix + 0);
	dummy.values[4] = *(matrix + 1);
	dummy.values[8] = *(matrix + 2);
	dummy.values[12] = *(matrix + 3);
	
	dummy.values[1] = *(matrix + 4);
	dummy.values[5] = *(matrix + 5);
	dummy.values[9] = *(matrix + 6);
	dummy.values[13] = *(matrix + 7);
	
	dummy.values[2] = *(matrix + 8);
	dummy.values[6] = *(matrix + 9);
	dummy.values[10] = *(matrix + 10);
	dummy.values[14] = *(matrix + 11);
	
	dummy.values[3] = *(matrix + 12);
	dummy.values[7] = *(matrix + 13);
	dummy.values[11] = *(matrix + 14);
	dummy.values[15] = *(matrix + 15);
	
	
	if(matrixMode == MGL_MODELVIEW)
	{
		multiply(currentModel, dummy);
	}
	else if(matrixMode == MGL_PROJECTION)
	{
		//cout << "currentProj: " << currentProj << endl;
		multiply(dummy, currentProj);
		//cout << "currentProj: " << currentProj << endl;
	}
	else
	{
		MGL_ERROR("Invalid");
		exit(1);
	}
}

/**
 * Multiply the current matrix by the translation matrix
 * for the translation vector given by (x, y, z).
 */
void mglTranslate(MGLfloat x,
                  MGLfloat y,
                  MGLfloat z)
{
	mat<MGLfloat, 4, 4> translate;
	translate.make_zero();
	
	translate.values[0] = 1;
	translate.values[5] = 1;
	translate.values[10] = 1;
	translate.values[3] = x;
	translate.values[7] = y;
	translate.values[11] = z;
	translate.values[15] = 1;
	
	
	if(matrixMode == MGL_PROJECTION)
	{
	
		multiply(translate, currentProj);

	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		multiply(translate, currentModel);
	}
	else
	{
		MGL_ERROR("Invalid matrix type");
		exit(1);
	}
}

/**
 * Multiply the current matrix by the rotation matrix
 * for a rotation of (angle) degrees about the vector
 * from the origin to the point (x, y, z).
 */
void mglRotate(MGLfloat angle,
               MGLfloat x,
               MGLfloat y,
               MGLfloat z)
{
	mat4 rotate;
	rotate.make_zero();
	
	MGLfloat norm = sqrt(x*x + y*y + z*z);
	x = x/norm;
	y = y/norm;
	z = z/norm;
	
	MGLfloat cosine = cos((angle*3.14159) / 180);
	MGLfloat sine = sin((angle*3.14159) / 180);
	
	rotate.values[0] = ((x*x) * (1 - cosine)) + cosine;
	rotate.values[1] = ((x*y) * (1 - cosine)) - z*sine;
	rotate.values[2] = ((x*z) * (1 - cosine)) + y*sine;
	
	rotate.values[4] = ((y*x) * (1 - cosine)) + z*sine;
	rotate.values[5] = ((y*y) * (1 - cosine)) + cosine;
	rotate.values[6] = ((y*z) * (1 - cosine)) - x*sine;
	
	rotate.values[8] = ((z*x) * (1 - cosine)) - y*sine;
	rotate.values[9] = ((z*y) * (1 - cosine)) + x*sine;
	rotate.values[10] = ((z*z) * (1 - cosine)) + cosine;
	
	rotate.values[15] = 1;
	
	if(matrixMode == MGL_PROJECTION)
	{
	
		multiply(rotate, currentProj);

	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		multiply(rotate, currentModel);
	}
	else
	{
		MGL_ERROR("Invalid matrix type");
		exit(1);
	}
}

/**
 * Multiply the current matrix by the scale matrix
 * for the given scale factors.
 */
void mglScale(MGLfloat x,
              MGLfloat y,
              MGLfloat z)
{
	mat<MGLfloat, 4, 4> scalar;
	scalar.make_zero();
	
	scalar.values[0] = x;
	scalar.values[5] = y;
	scalar.values[10] = z;
	scalar.values[15] = 1;
	
	
	if(matrixMode == MGL_PROJECTION)
	{
	
		multiply(scalar, currentProj);

	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		multiply(scalar, currentModel);
	}
	else
	{
		MGL_ERROR("Invalid matrix type");
		exit(1);
	}
}

/**
 * Multiply the current matrix by the perspective matrix
 * with the given clipping plane coordinates.
 */
void mglFrustum(MGLfloat left,
                MGLfloat right,
                MGLfloat bottom,
                MGLfloat top,
                MGLfloat near,
                MGLfloat far)
{
	mat<MGLfloat, 4, 4> frustum;
	frustum.make_zero();
	
	
	frustum.values[0] = (2 * near) / (right - left);
	frustum.values[2] = (right + left) / (right - left);
	
	frustum.values[5] = (2 * near) / (top - bottom);
	frustum.values[6] = (top + bottom) / (top - bottom);
	
	frustum.values[10] = -((far + near) / (far - near));
	frustum.values[11] = (-2 * far * near) / (far - near);
	
	frustum.values[14] = -1;
	
	if(matrixMode == MGL_PROJECTION)
	{
	
		multiply(frustum, currentProj);

	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		multiply(frustum, currentModel);
	}
	else
	{
		MGL_ERROR("Invalid matrix type");
		exit(1);
	}
	
	//cout << "frust: " << frustum << endl;
}

/**
 * Multiply the current matrix by the orthographic matrix
 * with the given clipping plane coordinates.
 */
void mglOrtho(MGLfloat left,
              MGLfloat right,
              MGLfloat bottom,
              MGLfloat top,
              MGLfloat near,
              MGLfloat far)
{
	// ortho init
	mat<MGLfloat, 4, 4> ortho;
	ortho.make_zero();
	ortho.values[0] = 2 / (right - left);
	ortho.values[5] = 2 / (top - bottom);
	ortho.values[10] = -2 / (far - near);
	ortho.values[15] = 1;
	
	ortho.values[3] = -(right + left)/(right - left);
	ortho.values[7] = -(top + bottom)/(top - bottom);
	ortho.values[11] = -(far + near)/(far - near);
	
	mat<MGLfloat, 4, 4> result;
	result.make_zero();
	
	FAR = far;
	NEAR = near;
	
	if(matrixMode == MGL_PROJECTION)
	{
		multiply(ortho, currentProj);
	}
	else if(matrixMode == MGL_MODELVIEW)
	{
		multiply(ortho, currentModel);
	}
	else
	{
		MGL_ERROR("Invalid matrix type");
		exit(1);
	}
	return;
	
	
	
}

/**
 * Set the current color for drawn shapes.
 */
void mglColor(MGLfloat red,
              MGLfloat green,
              MGLfloat blue)
{
	RGB[0] = red;
	RGB[1] = green;
	RGB[2] = blue;
}
