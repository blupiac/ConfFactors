// --------------------------------------------------------------------------
// Copyright(C) 2009-2016
// Tamy Boubekeur
// 
// Permission granted to use this code only for teaching projects and 
// private practice.
//
// Do not distribute this code outside the teaching assignements.                                                                           
// All rights reserved.                                                       
// --------------------------------------------------------------------------

#include <GL/glew.h>
#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>

#define _USE_MATH_DEFINES
#include <math.h>

#include "Vec3.h"
#include "Camera.h"
#include "Mesh.h"
#include "GLProgram.h"
#include "Exception.h"

using namespace std;

//#define USE_SHADER
//#define USE_BVH
//#define USE_NPR


//#define GRADIENT
//#define DEBUG_LAP
//#define DEBUG_GAUSS

static const unsigned int DEFAULT_SCREENWIDTH = 1024;
static const unsigned int DEFAULT_SCREENHEIGHT = 768;
static const string DEFAULT_MESH_FILE ("models/pierre.off");

// Rayons envoyes en AO, et portee maximale pour une intersection
static const unsigned int AO_SAMPLES = 10;
static const float AO_RADIUS = 0.7;

static const string appTitle ("Informatique Graphique & Realite Virtuelle - Travaux Pratiques - Algorithmes de Rendu");
static const string myName ("Bernard Lupiac");
static GLint window;
static unsigned int FPS = 0;
static bool fullScreen = false;

// Rend la BBox du BVH visible
#ifdef USE_BVH
static bool bboxVisible = false;
#endif

// Mode de rendu en CPU
static unsigned int mode = 4;

// Numero de lumieres affichees pour rendu CPU
static unsigned int nlights = 1;

// Definition du sol, 0 = rien, 1 = 20*20 points, 2 = 200*200 points
static unsigned int floorType = 0;

static Camera camera;
static Mesh mesh;
GLProgram * glProgram;
GLuint defaultShader;

Vec3f red = Vec3f(0xFF, 0x41, 0x36);
Vec3f yellow = Vec3f(0xFF, 0xDC, 0x00);
Vec3f green = Vec3f(0x2E, 0xCC, 0x40);
Vec3f blue = Vec3f(0x00, 0x00, 0xFF);
Vec3f indigo = Vec3f(0x00, 0x1F, 0x3F);
Vec3f purple = Vec3f(0xB1, 0x0D, 0xC9);

//-----------------------------------------------------------------------------
//------------------------------LightSource------------------------------------
//-----------------------------------------------------------------------------

// Utilise pour indexes qui representent des coordonees de position
enum Dim { x, y, z };
// Utilise pour indexes qui representent des coordonees de couleur
enum Col { r, g, b };

//-----------------------------------------------------------------------------
//-----------------------------------Ray---------------------------------------
//-----------------------------------------------------------------------------

// Represente un rayon, avec point de depart et direction
class Ray {
	public:
	Vec3f origin;
	Vec3f direction;
	
	inline Ray (Vec3f o, Vec3f dir)
	{origin = o; direction = normalize(dir);}
	
	// Retourne si le rayon intersecte le triangle donne par les points
	bool intersectsPoly(const Triangle tri) {
		Vec3f e0 = mesh.positions()[tri[1]] - mesh.positions()[tri[0]];
		Vec3f e1 = mesh.positions()[tri[2]] - mesh.positions()[tri[0]];
		Vec3f n = cross(e0, e1) / length(cross(e0, e1));
		Vec3f q = cross(direction, e1);
		float a = dot(e0, q);
		
		if(fabs(a) < 0.0001)
			return false;
			
		Vec3f s = (origin - mesh.positions()[tri[0]]) / a;
		Vec3f r = cross(s, e0);
		
		float b0 = dot(s, q);
		float b1 = dot(r, direction);
		float b2 = 1 - b0 - b1;
		
		if(b0 < 0 || b1 < 0 || b2 < 0)
			return false;
			
		float t = dot(e1, r);	
		
		if(t >= 0)
			return true;
		
		return false;
	}
	
};

//-----------------------------------------------------------------------------
//-----------------------------------Vec4--------------------------------------
//-----------------------------------------------------------------------------

template <typename T>
class Vec4
{

public:
	inline Vec4 (void)  { m_p[x] = m_p[y] = m_p[z] = m_p[3] = 0.0; }

	inline Vec4 (T p0, T p1, T p2, T p3) {
		m_p[0] = p0;
		m_p[1] = p1;
		m_p[2] = p2;
		m_p[3] = p3;
	};
	
	inline Vec4 & init (T x, T y, T z, T w) {
		m_p[0] = x;
		m_p[1] = y;
		m_p[2] = z;
		m_p[3] = w;
		return (*this);
  };

	inline Vec4 (const Vec4 & v) {
		init (v[0], v[1], v[2], v[3]);
	}
	
  inline T& operator[] (int Index) {
    return (m_p[Index]);
  };

  inline const T& operator[] (int Index) const {
    return (m_p[Index]);
  };

  inline Vec4& operator= (const Vec4 & p) {
    m_p[0] = p[0];
    m_p[1] = p[1];
    m_p[2] = p[2];
    m_p[3] = p[3];
    return (*this);
  };

protected:
  T m_p[4];

};
typedef Vec4<float> Vec4f;

//-----------------------------------------------------------------------------
//--------------------------AxisAlignedBoundingBox-----------------------------
//-----------------------------------------------------------------------------

// Represente bounding box, avec points min et max
class AxisAlignedBoundingBox {
	
public:
	inline AxisAlignedBoundingBox ()	{}
	
	// Donne des triangles, trouve leurs bounding boxes
	inline AxisAlignedBoundingBox(vector<Vec3f>& pos, vector<Triangle>& tri)
	{
		float max_x, max_y, max_z, min_x, min_y, min_z;
		
		max_x = max_y = max_z = pow(10,-10);
		min_x = min_y = min_z = pow(10,10);
		
		for(unsigned int i = 0; i < tri.size(); i++) {
			for(unsigned int j = 0; j < 3; j++) {
				if(pos[tri[i][j]][x] > max_x)
					max_x = pos[tri[i][j]][x];
				if(pos[tri[i][j]][x] < min_x)
					min_x = pos[tri[i][j]][x];

				if(pos[tri[i][j]][y] > max_y)
					max_y = pos[tri[i][j]][y];
				if(pos[tri[i][j]][y] < min_y)
					min_y = pos[tri[i][j]][y];
	
				if(pos[tri[i][j]][z] > max_z)
					max_z = pos[tri[i][j]][z];
				if(pos[tri[i][j]][z] < min_z)
					min_z = pos[tri[i][j]][z];
			}
		}
		
		uppLefFro = Vec3f(min_x, max_y, min_z);
		uppLefBac = Vec3f(min_x, max_y, max_z);
		uppRigFro = Vec3f(max_x, max_y, min_z);
		uppRigBac = Vec3f(max_x, max_y, max_z);
		lowLefFro = Vec3f(min_x, min_y, min_z);
		lowLefBac = Vec3f(min_x, min_y, max_z);
		lowRigFro = Vec3f(max_x, min_y, min_z);
		lowRigBac = Vec3f(max_x, min_y, max_z);	
	}
	
	// Retourne si un point est dans la bounding box
	inline bool inside (Vec3f& point)
	{
		if(	point[x] <= lowRigBac[x] && point[x] >= uppLefFro[x] &&
			point[y] <= uppLefFro[y] && point[y] >= lowRigBac[y] &&
			point[z] <= lowRigBac[z] && point[z] >= lowRigFro[z])
				return true;
				
		return false;
	}
	
	// Algorithme pour voir si un rayon intersecte la bounding box
	// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
	inline bool intersects(const Ray &r, float radius)
	{
		float tmin = (lowLefFro[x] - r.origin[x]) / r.direction[x]; 
		float tmax = (uppRigBac[x] - r.origin[x]) / r.direction[x]; 
	 
		if (tmin > tmax) swap(tmin, tmax); 
	 
		float tymin = (lowLefFro[y] - r.origin[y]) / r.direction[y]; 
		float tymax = (uppRigBac[y] - r.origin[y]) / r.direction[y]; 
	 
		if (tymin > tymax) swap(tymin, tymax); 
	 
		if ((tmin > tymax) || (tymin > tmax)) 
		    return false; 
	 
		if (tymin > tmin) 
		    tmin = tymin; 
	 
		if (tymax < tmax) 
		    tmax = tymax; 
	 
		float tzmin = (lowLefFro[z] - r.origin[z]) / r.direction[z]; 
		float tzmax = (uppRigBac[z] - r.origin[z]) / r.direction[z]; 
	 
		if (tzmin > tzmax) swap(tzmin, tzmax); 
	 
		if ((tmin > tzmax) || (tzmin > tmax)) 
		    return false; 
	 
		if (tzmin > tmin) 
		    tmin = tzmin; 
	 
		if (tzmax < tmax) 
		    tmax = tzmax; 
	 	
	 	if(tmin < radius && tmax > 0.0)
			return true; 
		
		return false;
	}
	
	// Affiche la BBox sur le programe
	inline void drawBBox()
	{
	
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
		
		// BACK
		glBegin(GL_POLYGON);
		glColor3f( 1.0, 0.0, 0.0 );
		glVertex3f( lowRigBac[x], lowRigBac[y], lowRigBac[z] );
		glVertex3f( uppRigBac[x], uppRigBac[y], uppRigBac[z] );
		glVertex3f( uppLefBac[x], uppLefBac[y], uppLefBac[z] );
		glVertex3f( lowLefBac[x], lowLefBac[y], lowLefBac[z] );
		glEnd();
		 
		// RIGHT
		glBegin(GL_POLYGON);
		glColor3f( 1.0, 0.0, 0.0 );
		glVertex3f( lowRigFro[x], lowRigFro[y], lowRigFro[z] );
		glVertex3f( uppRigFro[x], uppRigFro[y], uppRigFro[z] );
		glVertex3f( uppRigBac[x], uppRigBac[y], uppRigBac[z] );
		glVertex3f( lowRigBac[x], lowRigBac[y], lowRigBac[z] );
		glEnd();
		 
		// LEFT
		glBegin(GL_POLYGON);
		glColor3f( 1.0, 0.0, 0.0 );
		glVertex3f( lowLefBac[x], lowLefBac[y], lowLefBac[z] );
		glVertex3f( uppLefBac[x], uppLefBac[y], uppLefBac[z] );
		glVertex3f( uppLefFro[x], uppLefFro[y], uppLefFro[z] );
		glVertex3f( lowLefFro[x], lowLefFro[y], lowLefFro[z] );
		glEnd();
		 
		// TOP
		glBegin(GL_POLYGON);
		glColor3f( 1.0, 0.0, 0.0 );
		glVertex3f( uppRigBac[x], uppRigBac[y], uppRigBac[z] );
		glVertex3f( uppRigFro[x], uppRigFro[y], uppRigFro[z] );
		glVertex3f( uppLefFro[x], uppLefFro[y], uppLefFro[z] );
		glVertex3f( uppLefBac[x], uppLefBac[y], uppLefBac[z] );
		glEnd();
		 
		// BOTTOM
		glBegin(GL_POLYGON);
		glColor3f( 1.0, 0.0, 0.0 );
		glVertex3f( lowRigFro[x], lowRigFro[y], lowRigFro[z] );
		glVertex3f( lowRigBac[x], lowRigBac[y], lowRigBac[z] );
		glVertex3f( lowLefBac[x], lowLefBac[y], lowLefBac[z] );
		glVertex3f( lowLefFro[x], lowLefFro[y], lowLefFro[z] );
		glEnd();
		
		glFlush();
		//glutSwapBuffers();
		
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	
	}

	Vec3f uppLefFro, uppLefBac, uppRigFro, uppRigBac, 
			lowLefFro, lowLefBac, lowRigFro, lowRigBac;	
	
};

// Retourne centroide du triangle
Vec3f getCentroid(const vector<Vec3f>& pos, const Triangle& tri)
{		
	return Vec3f(	(pos.at(tri[0])[x]+pos.at(tri[1])[x]+pos.at(tri[2])[x])/3.0f,
					(pos.at(tri[0])[y]+pos.at(tri[1])[y]+pos.at(tri[2])[y])/3.0f,
					(pos.at(tri[0])[z]+pos.at(tri[1])[z]+pos.at(tri[2])[z])/3.0f);
}

// Fonction pour ordonner les triangles selon la position de leur centroide
// par rapport  la dimention specifiee
bool SortFunction(const Triangle& i, const Triangle& j, const vector<Vec3f>& pos, const Dim dim) {
	Vec3f centri, centrj;
	
	centri = getCentroid(pos, i);
	centrj = getCentroid(pos, j);
	
    return centri[dim] < centrj[dim];
}

// Construit une SortFonction pour donner a une fonction Sort
class Comparator {
public:
    Comparator(const vector<Vec3f>& pos, const Dim dim) { this->pos = pos; this->dim = dim;}

    bool operator () (const Triangle& i, const Triangle& j) {
    	return SortFunction(i, j, pos, dim);
    }

    vector<Vec3f> pos;
    Dim dim;
};
				
//-----------------------------------------------------------------------------
//-----------------------------------BVH---------------------------------------
//-----------------------------------------------------------------------------

// BVH pour accelerer appels d'AO et d'ombre
class BVH {

public:

	inline BVH ()	{}
	
	// Cree une BVH recursivement 
	inline BVH (vector<Vec3f>& pos, vector<Triangle> tri, Dim dim, bool sorted)
	{
		// Initialisation des variables
		feuille = drawn = false;
		this->bbox = AxisAlignedBoundingBox(pos, tri);
		
		// Cree sous noeuds
		if(tri.size() > 1)
		{			
			vector<Triangle> tri1, tri2;
			distributeSort(pos, tri, tri1, tri2, dim, sorted);
			tri.clear();
		
			leftChild = new BVH(pos, tri1, Dim((dim+1)%3), true);
			rightChild = new BVH(pos, tri2, Dim((dim+1)%3), true);	
		} // Condition pour que le noeud soit une feuille
		else if(tri.size() == 1)
		{
			feuille = true;
			leftChild = rightChild = NULL;
			triangle = tri[0];
		}
	}
	
	inline ~BVH ()
	{
		delete leftChild;
		delete rightChild;
	}
	
	inline bool isFeuille()
	{
		return feuille;
	}

	// Affiche la BBox du BVH courrant
	inline void drawBVH()
	{
		if(!drawn)
		{
			bbox.drawBBox();
		}
		else if(!feuille)
		{
			if(leftChild != NULL)
				leftChild->drawBVH();
			if(rightChild != NULL)
				rightChild->drawBVH();
		}
	}

	// Affice la BBox, un niveau de profondeur de plus a chaque appel
	inline void drawNextBVH()
	{
		if(!drawn)
		{
			drawn = true;
		}
		else if(!feuille)
		{
			if(leftChild != NULL)
				leftChild->drawNextBVH();
			if(rightChild != NULL)
				rightChild->drawNextBVH();
		}
	}
	
	// Retourne s'il y a intersection entre feuille ou sous-feuilles et rayon
	inline bool searchIntersection(vector<Vec3f>& pos, Ray& r, Triangle& resp, float radius)
	{
		// Si feuille, donne directement s'il y a intersection ou pas
		if(feuille)
		{
			if(r.intersectsPoly(triangle))
			{
				resp = triangle;
				return true;	
			}
			return false;
		} // Sinon, appel recursif sur le fils du noeud courrant qui intersecte le rayon
		else 
		{
			if(leftChild->getBBox().intersects(r, radius))
				leftChild->searchIntersection(pos, r, resp, radius);
			
			if(rightChild->getBBox().intersects(r, radius))
				rightChild->searchIntersection(pos, r, resp, radius);
				
			return false;
		}
	}
	
	inline AxisAlignedBoundingBox getBBox()
	{
		return bbox;
	}

private:

	Triangle triangle;
	bool feuille, drawn;
	AxisAlignedBoundingBox bbox;
	BVH* leftChild;
	BVH* rightChild;
	bool sorted = false;
	
	// Ordonne les triangles selon la position de leur centroide sur la dimention donnee
	inline void distributeSort(vector<Vec3f>& positions, vector<Triangle>& triangles,
							vector<Triangle>& tri1, vector<Triangle>& tri2, Dim dim, bool sorted)
	{
		if(!sorted)
		{
			std::sort(triangles.begin(), triangles.end(), Comparator(positions, dim));
			sorted = true;
		}

		tri1.assign(triangles.begin(), triangles.begin() + triangles.size()/2);
		tri2.assign(triangles.begin() + triangles.size()/2, triangles.end());
	}
	
};

static std::vector<Vec4f> colorResponses; // Cached per-vertex color response, updated at each frame
static std::vector<bool> shadowResponses; 	// vecteur de booleens qui disent si un point est dans l'ombre
static std::vector<float> aoResponses; // vecteur de floats qui disent la quantite d'AO sur un point
#ifdef USE_BVH
static BVH* bvh; // BVH utilise dans le programme
#endif

//-----------------------------------------------------------------------------
//----------------------------------Usage--------------------------------------
//-----------------------------------------------------------------------------

void printUsage () {
	std::cerr << std::endl 

		 << appTitle << std::endl

         << "Author: " << myName << std::endl << std::endl

         << "Usage: ./main [<" << DEFAULT_MESH_FILE << ">]" << std::endl

         << "Commands:" << std::endl 

         << "------------------" << std::endl

         << " ?: Print help" << std::endl

		 << " w: Toggle wireframe mode" << std::endl 
		 << " n: reload mesh" << std::endl
		 << " y: put floor with 20x20 resolution (press again to have 200x200, and a 3rd time to disable)" << std::endl
		 << " r: recompute normals" << std::endl
		 
		 << std::endl
		 << " ====== TP1: Synthese d’Image : Eclairage et BRDF ======" << std::endl
		 << " 0: switch to 3 sources of color (press again to go back to 1 source)" << std::endl
		 << " 4: Lambert BRDF" << std::endl
		 << " 5: Blinn-Phong BRDF" << std::endl
		 << " 6: Cook-Torrance BRDF" << std::endl
		 << " 7: GGX BRDF" << std::endl
		 << " to use the shader, uncomment \"define SHADER\" in code" << std::endl
		 
		 << std::endl
		 << " ====== TP2: Synthese d’Image : Ombrage ======" << std::endl
		 << " s: show shadow" << std::endl
		 << " o: show ambient occlusion" << std::endl
		 << " to use BVH, uncomment \"define BVH\" in code" << std::endl
		 << " d: show bounding box of BVH (pressing again show bounding boxes of sons recursively)" << std::endl
		 
		 << std::endl
		 << " ====== TP3: Synthese d’Image : NPR ======" << std::endl
		 << " to use NPR shader, uncomment \"define NPR\" in code" << std::endl
		 
		 << std::endl
		 << " ====== TP4: Modélisation Geometrique : Filtrage ======" << std::endl
		 << " l: topologic laplacian filter with alpha 0.1" << std::endl
		 << " 1: geometric laplacian filter with alpha 0.1" << std::endl
		 << " 2: geometric laplacian filter with alpha 0.5" << std::endl
		 << " 3: geometric laplacian filter with alpha 1.0" << std::endl
		 
		 << std::endl
		 << " ====== TP5: Modélisation Geometrique : Simplification ======" << std::endl
		 << " 8: simplify with grid using 16 cubes" << std::endl
		 << " 9: simplify with octree using at most 14 points per node" << std::endl
		 << " obs: ne marche pas tres bien avec Monkey car tres simple deja" << std::endl
		 
		 << std::endl
		 << " ====== TP6: Modelisation Geometrique : Subdivision ======" << std::endl
		 << " z: subdivide surface" << std::endl

		 << std::endl
         << " <drag>+<left button>: rotate model" << std::endl 
         << " <drag>+<right button>: move model" << std::endl

         << " <drag>+<middle button>: zoom" << std::endl

         << " q, <esc>: Quit" << std::endl << std::endl; 

}

//-----------------------------------------------------------------------------
//----------------------------------init---------------------------------------
//-----------------------------------------------------------------------------

void init (const char * modelFilename) {

	srand (time(NULL));

    glewExperimental = GL_TRUE;

    glewInit (); // init glew, which takes in charges the modern OpenGL calls (v>1.2, shaders, etc)

    glCullFace (GL_BACK);     // Specifies the faces to cull (here the ones pointing away from the camera)

    glEnable (GL_CULL_FACE); // Enables face culling (based on the orientation defined by the CW/CCW enumeration).

    glDepthFunc (GL_LESS); // Specify the depth test for the z-buffer

    glEnable (GL_DEPTH_TEST); // Enable the z-buffer in the rasterization

    glEnableClientState (GL_VERTEX_ARRAY);

    glEnableClientState (GL_NORMAL_ARRAY);

    glEnableClientState (GL_COLOR_ARRAY);

    glEnable (GL_NORMALIZE);

	glLineWidth (2.0); // Set the width of edges in GL_LINE polygon mode

    //glClearColor (0.0f, 0.0f, 0.0f, 1.0f); // Background color
    glClearColor (0.0f, 0.5f, 0.8f, 1.0f); // Background color    

	mesh.loadOFF (modelFilename);

    colorResponses.resize(mesh.positions().size());

    shadowResponses.resize(mesh.positions().size());

	for(unsigned int i = 0; i < shadowResponses.size(); i++ )
		shadowResponses[i] = false;

    aoResponses.resize (mesh.positions ().size ());

	for(unsigned int i = 0; i < aoResponses.size(); i++ )
		aoResponses[i] = 1.0;

    camera.resize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    
    #ifdef USE_BVH
    bvh = new BVH(mesh.positions(), mesh.triangles(), x, false);
    #endif
    
    #ifdef USE_SHADER
    try {
        #ifdef USE_NPR
        glProgram = GLProgram::genVFProgram ("Simple GL Program", "shader.vert", "toon.frag");
        #else
        glProgram = GLProgram::genVFProgram ("Simple GL Program", "shader.vert", "shader.frag");
        #endif
		glProgram->use (); // Activate the shader program
    } catch (Exception & e) {

        cerr << e.msg () << endl;

    }
    
   	glProgram->setUniform3f("lightPos", lights[0].pos()[x] , lights[0].pos()[y], lights[0].pos()[z]);
   	glProgram->setUniform3f("lightCol", lights[0].color()[r] , lights[0].color()[g], lights[0].color()[b]);
   	glProgram->setUniform1f("lightIntens", lights[0].intensity());
   	#ifndef USE_NPR
   	glProgram->setUniform1i("mode", mode);
   	#endif

   	#endif

}

// Calcule attenuation de la lumiere
float attLum (float dist, float intens) {
	float ac = 0.011;
	float al = 0.008;
	float aq = 0.002;
	
    return intens / (ac + al * dist + aq * pow(dist, 2));
}

float phong (Vec3f n, Vec3f wh) {
	float shinyness = 15.0f;
	
	if(dot(n, wh) < 0)
		return 0;
	
	return pow(dot(n, wh), shinyness);
}

float cooktorrance (Vec3f n, Vec3f wh, Vec3f wi, Vec3f w0) {
	
	float alpha = 0.8;
	float refrac = 0.3;
	
	float D = ( exp( (pow(dot(n,wh), 2) - 1) / (pow(alpha, 2) * pow(dot(n,wh), 2) )) ) / ( M_PI * pow(alpha,2) * pow( dot(n, wh) , 4) );
	float F = refrac + (1 - refrac) * pow(1 - max(0.0f, dot(wi, wh)), 5);
	float G = min(1.0f, min((2*dot(n,wh)*dot(n,wi))/(dot(w0,wh)), (2*dot(n,wh)*dot(n,w0))/(dot(w0,wh))));
	
	if(D < 0 || F < 0 || G < 0)
		return 0;
	
	return D*F*G;

}

float gsmith(Vec3f n, Vec3f w, float alpha)
{
	return (2*dot(n,w))/(dot(n,w) + sqrt(pow(alpha,2) + (1-pow(alpha, 2) * pow(dot(n, w),2))));
}

float gschlick(Vec3f n, Vec3f w, float alpha)
{
	float k = alpha * sqrt(2/M_PI);
	return dot(n,w)/(dot(n,w)*(1-k)+k);
}

float ddx (Vec3f n, Vec3f wh, Vec3f wi, Vec3f w0) {
	
	float alpha = 0.9;
	float refrac = 0.9;
	
	float D = pow(alpha,2) / ( M_PI * pow((1.0 + (pow(alpha,2)-1) * pow(dot(n, wh), 2)),2));
	float F = refrac + (1 - refrac) * pow(1 - max(0.0f, dot(wi, wh)), 5);
	float G = gsmith(n, wi, alpha) * gsmith(n, w0, alpha);
	
	if(D < 0 || F < 0 || G < 0)
		return 0;

	return D*F*G;

}

//-----------------------------------------------------------------------------
//---------------------updatePerVertexColorResponse----------------------------
//-----------------------------------------------------------------------------

// https://stackoverflow.com/questions/21835739/smooth-color-transition-algorithm

Vec3f colorGradient(Vec3f col1, Vec3f col2, float val)
{
	Vec3f result;

	result[r] = (col1[r] * (1 - val)) + (col2[r] * val);
	result[g] = (col1[g] * (1 - val)) + (col2[g] * val);
	result[b] = (col1[b] * (1 - val)) + (col2[b] * val);

	return result;
}

void updatePerVertexColorResponse () {

    for (unsigned int i = 0; i < colorResponses.size (); i++)
    {
    	#ifdef DEBUG_LAP	
    	//std::cout << mesh.normalizeLaplacian(i) << std::endl;
    	float response = mesh.normalizeLaplacian(i);
    	#elif defined(DEBUG_GAUSS)
    	//std::cout << mesh.normalizeConf(i) << std::endl;
    	float response = mesh.normalizeGausscurv(i);
    	#else
    	//std::cout << mesh.normalizeConf(i) << std::endl;
		float response = mesh.normalizeConf(i);
		#endif

		#ifdef GRADIENT
		Vec3f resCol;

		if(response < 0.2)
		{
			resCol = colorGradient(red, yellow, (response - 0.0) * 5);
			colorResponses[i] = Vec4f(resCol[r] , resCol[g] , resCol[b] , 0.0);
		}
		else if(response < 0.4)
		{
			resCol = colorGradient(yellow, green, (response - 0.2) * 5);
			colorResponses[i] = Vec4f(resCol[r] , resCol[g] , resCol[b] , 0.0);
		}
		else if(response < 0.6)
		{
			resCol = colorGradient(green, blue, (response - 0.4) * 5);
			colorResponses[i] = Vec4f(resCol[r] , resCol[g] , resCol[b] , 0.0);
		}
		else if(response < 0.8)
		{
			resCol = colorGradient(blue, indigo, (response - 0.6) * 5);
			colorResponses[i] = Vec4f(resCol[r] , resCol[g] , resCol[b] , 0.0);
		}
		else
		{
			resCol = colorGradient(indigo, purple, (response - 0.8) * 5);
			colorResponses[i] = Vec4f(resCol[r] , resCol[g] , resCol[b] , 0.0);
		}
		std::cout << "response: " << response <<  " color given: " << resCol << std::endl;
		#else
		colorResponses[i] = Vec4f(response * aoResponses[i],
								response * aoResponses[i],
								response * aoResponses[i],0.0);
		#endif
    }
}


//-----------------------------------------------------------------------------
//-----------------------------computePerVertexAO------------------------------
//-----------------------------------------------------------------------------

// Calcule AO sur chaque point en creant des rayons avec directions aleatoires
// et en voyant combient d'intersections il y a avec la geometrie

void computePerVertexAO (unsigned int samples, unsigned int radius) {

	for(unsigned int i = 0; i < mesh.positions().size() ; i++)
	{
		Vec3f origin = mesh.positions()[i] + mesh.normals()[i] * 0.001;	
		
		Vec3f orth1, orth2;
		mesh.normals()[i].getTwoOrthogonals(orth1, orth2);
		
		aoResponses[i] = 0;
		
		for(unsigned int sample = 0; sample < samples ; sample++)
		{
			double rnd1 = ((double)rand()/(double)RAND_MAX)*2 - 1;
			double rnd2 = ((double)rand()/(double)RAND_MAX)*2 - 1;
			
			Vec3f offset = orth1*rnd1 + orth2*rnd2;

			Ray raycast = Ray(origin, mesh.normals()[i] + offset);
	
			#ifdef USE_BVH
			Triangle inter = Triangle(0, 0, 0);
			bvh->searchIntersection(mesh.positions(), raycast, inter, radius);
			if(inter[0] == 0 && inter[1] == 0 && inter[2] == 0)
				aoResponses[i] += 1.0f/(float)samples;
			#else
			bool intersection = false;	
			for(unsigned int j = 0; j < mesh.triangles().size() ; j++)
			{		
				if( (length(mesh.positions()[mesh.triangles()[j][x]] - mesh.positions()[i]) < 1.0) &&
					raycast.intersectsPoly(mesh.triangles()[j]) )
				{
					intersection = true;
					break;
				}
			}
			
			if(!intersection)
			{
				aoResponses[i] += 1.0f/(float)samples;
			}
			#endif
		}
	}
	
	cout << "ambient occlusion done \n";
}

//-----------------------------------------------------------------------------
//---------------------------------renderScene---------------------------------
//-----------------------------------------------------------------------------

void renderScene () {
	
	#ifdef USE_BVH
	if(bboxVisible)
		bvh->drawBVH();
	#endif

    updatePerVertexColorResponse ();
    glVertexPointer (3, GL_FLOAT, sizeof (Vec3f), (GLvoid*)(&(mesh.positions()[0])));

    glNormalPointer (GL_FLOAT, 3*sizeof (float), (GLvoid*)&(mesh.normals()[0]));
	glColorPointer (4, GL_FLOAT, sizeof (Vec4f), (GLvoid*)(&(colorResponses[0])));
    glDrawElements (GL_TRIANGLES, 3*mesh.triangles().size(), GL_UNSIGNED_INT, (GLvoid*)((&mesh.triangles()[0])));
}

void reshape(int w, int h) {
    camera.resize (w, h);
}

void display () {
    glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    camera.apply (); 
    renderScene ();
    glFlush ();
    glutSwapBuffers (); 
}

void key (unsigned char keyPressed, int x, int y) {

    switch (keyPressed) {
    case 'f':
        if (fullScreen) {
            glutReshapeWindow (camera.getScreenWidth (), camera.getScreenHeight ());
            fullScreen = false;
        } else {
            glutFullScreen ();
            fullScreen = true;
        }      
        break;
    case '0':
        if(nlights == 1)
        	nlights = 3;
        else
        	nlights = 1;
    	break;
    case '1':
        mesh.laplacianFilterGeom(0.1);
        break;
    case '2':
        mesh.laplacianFilterGeom(0.5);
        break;
    case '3':
        mesh.laplacianFilterGeom(1.0);
        break;
    case '4':
    	// lambert
        mode = 1;
        #ifdef USE_SHADER
        glProgram->setUniform1i("mode", mode);
        #endif
        break;
    case '5':
    	// phong
        mode = 2;
        #ifdef USE_SHADER
        glProgram->setUniform1i("mode", mode);
        #endif
        break;
    case '6':
    	// cook torrance
        mode = 3;
        #ifdef USE_SHADER
        glProgram->setUniform1i("mode", mode);
        #endif
        break;
    case '7':
    	// ggx
        mode = 4;
        #ifdef USE_SHADER
        glProgram->setUniform1i("mode", mode);
        #endif
        break;
    case '8':
        mesh.simplify(16);
        break;
    case '9':
        mesh.simplifyAdaptive(14);  
        break;
    case 'z':
        mesh.subdivide ();
        colorResponses.resize(mesh.positions().size());
        
        for(unsigned int i = 0; i < shadowResponses.size(); i++ )
			colorResponses[i] = Vec4f(0.0,0.0,0.0,0.0);

		shadowResponses.resize(mesh.positions().size());

		for(unsigned int i = 0; i < shadowResponses.size(); i++ )
			shadowResponses[i] = false;

		aoResponses.resize (mesh.positions ().size ());

		for(unsigned int i = 0; i < aoResponses.size(); i++ )
			aoResponses[i] = 1.0;
			
        break;
    case 'o':
        computePerVertexAO (AO_SAMPLES, AO_RADIUS); 
        break;
    case 'd':
    	#ifdef USE_BVH
    	if(!bboxVisible)
    		bboxVisible = true;
    	else
	    	bvh->drawNextBVH();
	    #endif
    	break;
    case 'l':
    	mesh.laplacianFilter(0.1);
    	break;
    case 'n':
    	mesh.loadOFF (DEFAULT_MESH_FILE);
		colorResponses.resize(mesh.positions().size());
		shadowResponses.resize(mesh.positions().size());

		for(unsigned int i = 0; i < shadowResponses.size(); i++ )
			shadowResponses[i] = false;

		aoResponses.resize (mesh.positions ().size ());

		for(unsigned int i = 0; i < aoResponses.size(); i++ )
			aoResponses[i] = 1.0;
			
    	break;
    case 'r':
    	mesh.recomputeNormals();
    	break;
    case 'y':
    	if(floorType == 0)	{    	
			mesh.addFloor(1);
			floorType = 1;
		}
		else if(floorType == 1)	{
			mesh.removeFloor(1);
			mesh.addFloor(2);
			floorType = 2;
		}
		else {
			mesh.removeFloor(2);
			floorType = 0;
		}
		
		colorResponses.resize(mesh.positions().size());
		shadowResponses.resize(mesh.positions().size());

		for(unsigned int i = 0; i < shadowResponses.size(); i++ )
			shadowResponses[i] = false;

		aoResponses.resize (mesh.positions ().size ());

		for(unsigned int i = 0; i < aoResponses.size(); i++ )
			aoResponses[i] = 1.0;
			
    	break;
    case 'q':
    case 27:
        exit (0);
        break;
    case 'w':
        GLint mode[2];
		glGetIntegerv (GL_POLYGON_MODE, mode);
		glPolygonMode (GL_FRONT_AND_BACK, mode[1] ==  GL_FILL ? GL_LINE : GL_FILL);
        break;
        break;
    default:
        printUsage ();
        break;
    }
}

void mouse (int button, int state, int x, int y) {
    camera.handleMouseClickEvent (button, state, x, y);
}

void motion (int x, int y) {
    camera.handleMouseMoveEvent (x, y);
}

void idle () {
    static float lastTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    static unsigned int counter = 0;
    counter++;
    float currentTime = glutGet ((GLenum)GLUT_ELAPSED_TIME);
    if (currentTime - lastTime >= 1000.0f) {
        FPS = counter;
        counter = 0;
        static char winTitle [128];
        unsigned int numOfTriangles = mesh.triangles ().size ();
        sprintf (winTitle, "Number Of Triangles: %d - FPS: %d", numOfTriangles, FPS);
        string title = appTitle + " - By " + myName  + " - " + winTitle;
        glutSetWindowTitle (title.c_str ());
        lastTime = currentTime;
    }
    glutPostRedisplay (); 
}

int main (int argc, char ** argv) {
    if (argc > 2) {
        printUsage ();
        exit (1);
    }
    glutInit (&argc, argv);
    glutInitDisplayMode (GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutInitWindowSize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    window = glutCreateWindow (appTitle.c_str ());
    init (argc == 2 ? argv[1] : DEFAULT_MESH_FILE.c_str ());
    glutIdleFunc (idle);
    glutReshapeFunc (reshape);
    glutDisplayFunc (display);
    glutKeyboardFunc (key);
    glutMotionFunc (motion);
    glutMouseFunc (mouse);
    printUsage ();  
    glutMainLoop ();
    return 0;
}
