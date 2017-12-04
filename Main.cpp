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


#define GRADIENT
//#define DEBUG_LAP
//#define DEBUG_GAUSS

static const unsigned int DEFAULT_SCREENWIDTH = 1024;
static const unsigned int DEFAULT_SCREENHEIGHT = 768;
static const std::string DEFAULT_MESH_FILE ("models/armadillo1.off");

// Rayons envoyes en AO, et portee maximale pour une intersection
static const unsigned int AO_SAMPLES = 10;
static const float AO_RADIUS = 0.7;

static const std::string appTitle ("Informatique Graphique & Realite Virtuelle - Travaux Pratiques - Algorithmes de Rendu");
static const std::string myName ("Bernard Lupiac");
static GLint window;
static unsigned int FPS = 0;
static bool fullScreen = false;


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

static std::vector<Vec4f> colorResponses; // Cached per-vertex color response, updated at each frame

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

    camera.resize (DEFAULT_SCREENWIDTH, DEFAULT_SCREENHEIGHT);
    
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
	float F = refrac + (1 - refrac) * pow(1 - std::max(0.0f, dot(wi, wh)), 5);
	float G = std::min(1.0f, std::min((2*dot(n,wh)*dot(n,wi))/(dot(w0,wh)), (2*dot(n,wh)*dot(n,w0))/(dot(w0,wh))));
	
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
	float F = refrac + (1 - refrac) * pow(1 - std::max(0.0f, dot(wi, wh)), 5);
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
			colorResponses[i] = Vec4f(resCol[r]/255.0 , resCol[g]/255.0 , resCol[b]/255.0 , 0.0);
		}
		else if(response < 0.4)
		{
			resCol = colorGradient(yellow, green, (response - 0.2) * 5);
			colorResponses[i] = Vec4f(resCol[r]/255.0 , resCol[g]/255.0 , resCol[b]/255.0 , 0.0);
		}
		else if(response < 0.6)
		{
			resCol = colorGradient(green, blue, (response - 0.4) * 5);
			colorResponses[i] = Vec4f(resCol[r]/255.0 , resCol[g]/255.0 , resCol[b]/255.0 , 0.0);
		}
		else if(response < 0.8)
		{
			resCol = colorGradient(blue, indigo, (response - 0.6) * 5);
			colorResponses[i] = Vec4f(resCol[r]/255.0 , resCol[g]/255.0 , resCol[b]/255.0 , 0.0);
		}
		else
		{
			resCol = colorGradient(indigo, purple, (response - 0.8) * 5);
			colorResponses[i] = Vec4f(resCol[r]/255.0 , resCol[g]/255.0 , resCol[b]/255.0 , 0.0);
		}
		//std::cout << "response: " << response <<  " color given: " << resCol << std::endl;
		#else
		colorResponses[i] = Vec4f(response * aoResponses[i],
								response * aoResponses[i],
								response * aoResponses[i],0.0);
		#endif
    }
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
    case '8':
        mesh.simplify(16);
        break;
    case '9':
        mesh.simplifyAdaptive(14);  
        break;
    case 'z':
        mesh.subdivide ();
        colorResponses.resize(mesh.positions().size());
        
        for(unsigned int i = 0; i < colorResponses.size(); i++ )
			colorResponses[i] = Vec4f(0.0,0.0,0.0,0.0);
        break;
    case 'n':
    	mesh.loadOFF (DEFAULT_MESH_FILE);
		colorResponses.resize(mesh.positions().size());			
    	break;
    case 'r':
    	mesh.recomputeNormals();
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
        std::string title = appTitle + " - By " + myName  + " - " + winTitle;
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
