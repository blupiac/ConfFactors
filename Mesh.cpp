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

#include "Mesh.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cmath>

#define MED_X min_x + (max_x - min_x)/2
#define MED_Y min_y + (max_y - min_y)/2
#define MED_Z min_z + (max_z - min_z)/2

using namespace std;

struct Cube {
	float max_x;
	float min_x;
	float max_y;
	float min_y;
	float max_z;
	float min_z;

	unsigned int element;
	
} cube;

void Mesh::clear () {
    m_positions.clear ();
    m_normals.clear ();
    m_triangles.clear ();
    m_nneighbours.clear();
    middle_points.clear();
}

// Utilise pour indexes qui representent des coordonees de position
enum Dim { x, y, z };

//-----------------------------------------------------------------------------
//--------------------------------loadOFF--------------------------------------
//-----------------------------------------------------------------------------

void Mesh::loadOFF (const std::string & filename) {
    clear ();
	ifstream in (filename.c_str ());
    if (!in) 
    {
        exit (1);
    }
	string offString;
    unsigned int sizeV, sizeT, tmp;
    in >> offString >> sizeV >> sizeT >> tmp;
    m_positions.resize (sizeV);
    m_triangles.resize (sizeT);
    m_confFact.resize (sizeV);
    m_nneighbours.resize(sizeV);
    totalArea = 0; totalCurv = 0;
    for (unsigned int i = 0; i < sizeV; i++)
        in >> m_positions[i];
    int s;
    for (unsigned int i = 0; i < sizeT; i++) {
        in >> s;
        for (unsigned int j = 0; j < 3; j++)
        {
            in >> m_triangles[i][j];
        }
        // fills neighbour vector
        for (unsigned int j = 0; j < 3; j++)
        {
            m_nneighbours[m_triangles[i][j]].push_back(m_triangles[i]);
        }

        // should fill area vector
        // https://www.opengl.org/discussion_boards/showthread.php/159771-How-can-I-find-the-area-of-a-3D-triangle
        // uses a sqrt, might look for smth else later

        //BA vector
		Vec3f v1 = Vec3f(m_positions[m_triangles[i][1]][x] - m_positions[m_triangles[i][0]][x],
					m_positions[m_triangles[i][1]][y] - m_positions[m_triangles[i][0]][y],
					m_positions[m_triangles[i][1]][z] - m_positions[m_triangles[i][0]][z]);

		//BC vector
		Vec3f v2 = Vec3f(m_positions[m_triangles[i][2]][x] - m_positions[m_triangles[i][0]][x],
					m_positions[m_triangles[i][2]][y] - m_positions[m_triangles[i][0]][y],
					m_positions[m_triangles[i][2]][z] - m_positions[m_triangles[i][0]][z]);

		//cross prod
		Vec3f v3 = Vec3f(v1[y] * v2[z] - v1[z] * v2[y],
						v1[z] * v2[x] - v1[x] * v2[z],
						v1[x] * v2[y] - v1[y] * v2[x]);

		totalArea += 0.5 * sqrt(v3[x]*v3[x] + v3[y]*v3[y] + v3[z]*v3[z]);
    }
    in.close ();
    centerAndScaleToUnit ();
    recomputeNormals ();
}

//-----------------------------------------------------------------------------
//--------------------------calculateConfFact----------------------------------
//-----------------------------------------------------------------------------

void Mesh::calculateConfFact () 
{
	for (unsigned int i = 0; i < m_positions.size (); i++) 
	{
		float gaussCurv = getGaussCurv(i);
		float targCurv = getTargetCurv(i);

		float laplacian = getLaplacian(i);

		m_confFact[i] = (targCurv - gaussCurv) / laplacian;
	}
}

float Mesh::getGaussCurv(unsigned int pointIdx)
{
	std::vector<Triangle> neighbours = m_nneighbours[pointIdx];
	float curv;
	float totalAngle = 0;

	for(unsigned int i = 0; i < neighbours.size(); i++)
	{
		totalAngle += getAngle(neighbours[i], pointIdx);
	}

	// interior
	if(neighbours.size() == 6)
	{
		curv = 2 * M_PI - totalAngle;
	}
	// edge
	else if(neighbours.size() == 3)
	{
		curv = M_PI - totalAngle;
	}
	// how did this point get here?
	else
	{
		std::cout << "Point " << pointIdx << " has an unsupported number of neighbours." << std::endl;
		return 0;
	}

	totalCurv += curv;
	return curv;
}

// get angle between 3 points in a 3d space
float Mesh::getAngle(Triangle tri, int pointIdx)
{
	// https://stackoverflow.com/questions/19729831/angle-between-3-points-in-3d-space

	Vec3f v1, v2;

	// finding out which is the current point
	if(m_positions[tri[0]] == m_positions[pointIdx])
	{
		//BA vector
		v1 = Vec3f(m_positions[tri[1]][x] - m_positions[tri[0]][x],
					m_positions[tri[1]][y] - m_positions[tri[0]][y],
					m_positions[tri[1]][z] - m_positions[tri[0]][z]);

		//BC vector
		v2 = Vec3f(m_positions[tri[2]][x] - m_positions[tri[0]][x],
					m_positions[tri[2]][y] - m_positions[tri[0]][y],
					m_positions[tri[2]][z] - m_positions[tri[0]][z]);
	}
	else if(m_positions[tri[1]] == m_positions[pointIdx])
	{
		//BA vector
		v1 = Vec3f(m_positions[tri[0]][x] - m_positions[tri[1]][x],
					m_positions[tri[0]][y] - m_positions[tri[1]][y],
					m_positions[tri[0]][z] - m_positions[tri[1]][z]);

		//BC vector
		v2 = Vec3f(m_positions[tri[2]][x] - m_positions[tri[1]][x],
					m_positions[tri[2]][y] - m_positions[tri[1]][y],
					m_positions[tri[2]][z] - m_positions[tri[1]][z]);
	}
	else if(m_positions[tri[2]] == m_positions[pointIdx])
	{
		//BA vector
		v1 = Vec3f(m_positions[tri[0]][x] - m_positions[tri[2]][x],
					m_positions[tri[0]][y] - m_positions[tri[2]][y],
					m_positions[tri[0]][z] - m_positions[tri[2]][z]);

		//BC vector
		v2 = Vec3f(m_positions[tri[1]][x] - m_positions[tri[2]][x],
					m_positions[tri[1]][y] - m_positions[tri[2]][y],
					m_positions[tri[1]][z] - m_positions[tri[2]][z]);
	}
	else
	{
		std::cout << "Point " << pointIdx << " has an unsupported number of neighbours." << std::endl;
	}

	float v1mag = sqrt(v1[x] * v1[x] + v1[y] * v1[y] + v1[z] * v1[z]);
	Vec3f v1norm = Vec3f(v1[x] / v1mag, v1[y] / v1mag, v1[z] / v1mag);

	float v2mag = sqrt(v2[x] * v2[x] + v2[y] * v2[y] + v2[z] * v2[z]);
	Vec3f v2norm = Vec3f(v2[x] / v2mag, v2[y] / v2mag, v2[z] / v2mag);

	float res = v1norm[x] * v2norm[x] + v1norm[y] * v2norm[y] + v1norm[z] * v2norm[z];

	return  acos(res);
}

float Mesh::getTargetCurv(unsigned int i)
{
	std::vector<Triangle> tris = m_nneighbours[i];
	float area = 0;

	for(unsigned int i = 0; i < tris.size(); i++)
	{
		area += getArea(tris[i]);
	}

	return totalCurv * area / totalArea;
}

float Mesh::getArea(unsigned int i)
{
    // https://www.opengl.org/discussion_boards/showthread.php/159771-How-can-I-find-the-area-of-a-3D-triangle
    // uses a sqrt, might look for smth else later

    //BA vector
	Vec3f v1 = Vec3f(m_positions[m_triangles[i][1]][x] - m_positions[m_triangles[i][0]][x],
				m_positions[m_triangles[i][1]][y] - m_positions[m_triangles[i][0]][y],
				m_positions[m_triangles[i][1]][z] - m_positions[m_triangles[i][0]][z]);

	//BC vector
	Vec3f v2 = Vec3f(m_positions[m_triangles[i][2]][x] - m_positions[m_triangles[i][0]][x],
				m_positions[m_triangles[i][2]][y] - m_positions[m_triangles[i][0]][y],
				m_positions[m_triangles[i][2]][z] - m_positions[m_triangles[i][0]][z]);

	//cross prod
	Vec3f v3 = Vec3f(v1[y] * v2[z] - v1[z] * v2[y],
					v1[z] * v2[x] - v1[x] * v2[z],
					v1[x] * v2[y] - v1[y] * v2[x]);

	return 0.5 * sqrt(v3[x]*v3[x] + v3[y]*v3[y] + v3[z]*v3[z]);
}

float Mesh::getArea(Triangle tri)
{
    // https://www.opengl.org/discussion_boards/showthread.php/159771-How-can-I-find-the-area-of-a-3D-triangle
    // uses a sqrt, might look for smth else later

    //BA vector
	Vec3f v1 = Vec3f(m_positions[tri[1]][x] - m_positions[tri[0]][x],
				m_positions[tri[1]][y] - m_positions[tri[0]][y],
				m_positions[tri[1]][z] - m_positions[tri[0]][z]);

	//BC vector
	Vec3f v2 = Vec3f(m_positions[tri[2]][x] - m_positions[tri[0]][x],
				m_positions[tri[2]][y] - m_positions[tri[0]][y],
				m_positions[tri[2]][z] - m_positions[tri[0]][z]);

	//cross prod
	Vec3f v3 = Vec3f(v1[y] * v2[z] - v1[z] * v2[y],
					v1[z] * v2[x] - v1[x] * v2[z],
					v1[x] * v2[y] - v1[y] * v2[x]);

	return 0.5 * sqrt(v3[x]*v3[x] + v3[y]*v3[y] + v3[z]*v3[z]);
}

float Mesh::getLaplacian(unsigned int i)
{
	return 0;
}

//-----------------------------------------------------------------------------
//-------------------------calculateSignature----------------------------------
//-----------------------------------------------------------------------------

// TODO: memory effective function, orders triangles indexes by area size,
// calculates probability
void Mesh::calculateSignature () 
{
	initializeSignature(-99, 100);

	for (unsigned int i = 0; i < 5 * m_positions.size (); i++) 
	{
		int triangleIdx = getRandTri();
		Vec3f randPoint = getRandPoint(m_triangles[triangleIdx]);
		float randConfFactor = getConfFactor(randPoint, triangleIdx);

		incrSignature(randConfFactor);
	}

	for (unsigned int i = 0; i < 5 * m_positions.size (); i++) 
	{
		incrSignature(m_confFact[i]);
	}
}

void Mesh::initializeSignature(int min, int max)
{
	signature.resize (max - min);

	for(unsigned int i = 0; i < signature.size(); i ++)
	{
		Bin b;
		b.occur = 0;
		b.idx = min + i;
		signature[i] = b;
	}
}

int Mesh::getRandTri()
{
	return 0;
}

Vec3f Mesh::getRandPoint(Triangle tri)
{
	return Vec3f(0,0,0);
}

float Mesh::getConfFactor(Vec3f point, unsigned int triIdx)
{
	return 0;
}

void Mesh::incrSignature(float confFact)
{

}

//-----------------------------------------------------------------------------
//---------------------------recomputeNormals----------------------------------
//-----------------------------------------------------------------------------

void Mesh::recomputeNormals () {
    m_normals.clear ();
    m_normals.resize (m_positions.size (), Vec3f (0.f, 0.f, 0.f));
    for (unsigned int i = 0; i < m_triangles.size (); i++) {
        Vec3f e01 = m_positions[m_triangles[i][1]] -  m_positions[m_triangles[i][0]];
        Vec3f e02 = m_positions[m_triangles[i][2]] -  m_positions[m_triangles[i][0]];
        Vec3f n = cross (e01, e02);
        Vec3f cross = n;
        n.normalize ();
        for (unsigned int j = 0; j < 3; j++)
        {
        	// avec poids
            m_normals[m_triangles[i][j]] += n * std::asin(length(cross));
        }
    }
    for (unsigned int i = 0; i < m_normals.size (); i++)
        m_normals[i].normalize ();
}

//-----------------------------------------------------------------------------
//---------------------------laplacianFilters----------------------------------
//-----------------------------------------------------------------------------

void Mesh::laplacianFilter(float alpha) {
	
	for(unsigned int i = 0; i < m_nneighbours.size (); i++)
	{
		// Prends tous les points voisins au point
		std::set<unsigned int> voisins = getVoisins(m_nneighbours[i], i);

		// Calcule barycentre de ces points
		Vec3f bar = calculerBarycentre(voisins);

		// Repositionne point
		m_positions[i] = m_positions[i]*(1-alpha) + bar*alpha;
    }
	
	recomputeNormals ();
}

// Idem aue laplacien normal, sauf pour poids lors du calcul du barycentre
void Mesh::laplacianFilterGeom(float alpha) {
	
	for(unsigned int i = 0; i < m_nneighbours.size (); i++)

	{
		std::set<unsigned int> voisins = getVoisins(m_nneighbours[i], i);

		Vec3f bar = calculerBarycentreGeom(i, voisins, m_nneighbours[i]);

		m_positions[i] = m_positions[i]*(1-alpha) + bar*alpha;

    }
	recomputeNormals ();
}

// Prends min et max d'une mesh
void Mesh::getMaxMin(float& max_x, float& max_y, float& max_z, float& min_x,
				float& min_y, float& min_z)
{
		
	max_x = max_y = max_z = pow(10,-10);
	min_x = min_y = min_z = pow(10,10);
	
	for(unsigned int i = 0; i < m_positions.size(); i++) {
		if(m_positions[i][0] > max_x)
			max_x = m_positions[i][0];
		if(m_positions[i][0] < min_x)
			min_x = m_positions[i][0];

		if(m_positions[i][1] > max_y)
			max_y = m_positions[i][1];
		if(m_positions[i][1] < min_y)
			min_y = m_positions[i][1];

		if(m_positions[i][2] > max_z)
			max_z = m_positions[i][2];
		if(m_positions[i][2] < min_z)
			min_z = m_positions[i][2];
	}
}

// Retourne tous les points de la mesh dans les limites donnes
std::vector<unsigned int> Mesh::allInside(float max_x, float max_y, float max_z,
								float min_x, float min_y, float min_z)
{
	std::vector<unsigned int> inside;

	for(unsigned int n = 0; n < m_positions.size(); n++) {
				
		if(	m_positions[n][0] < max_x && m_positions[n][0] > min_x &&
			m_positions[n][1] < max_y && m_positions[n][1] > min_y &&
			m_positions[n][2] < max_z && m_positions[n][2] > min_z)
		{
			inside.push_back(n);
		}
	}

	return inside;
	
}

//-----------------------------------------------------------------------------
//-----------------------------simplifySubOct----------------------------------
//-----------------------------------------------------------------------------

void Mesh::simplifySubOct(float max_x, float max_y, float max_z,
							float min_x, float min_y, float min_z,
							unsigned int n, std::vector<Vec3f>& newPositions, 
							std::vector<Vec3f>& newNormals, std::vector<unsigned int>& oldToNew)
{
	// Si octree ne peut pas etre subdivisee
	std::vector<unsigned int> inside = allInside(max_x, max_y, max_z, min_x, min_y, min_z);
	if(inside.size() < n && inside.size() != 0)
	{
		Vec3f posMoyen, normMoyen;
		posMoyen = normMoyen = Vec3f(0.0,0.0,0.0);
	
		// faire moyenne des points dans la bbox
		for(std::vector<unsigned int>::iterator it = inside.begin();
				it != inside.end(); it++) {
			
			posMoyen += m_positions[*it];
			normMoyen += m_normals[*it];
			
			// Change point dans le vecteur des indices
			oldToNew[*it] = newPositions.size();
			
		}
		
		posMoyen /= (float) inside.size();
		normalize(normMoyen);
		
		// Ajouter comme point de la mesh
		newPositions.push_back(posMoyen);
		newNormals.push_back(normMoyen);
	}
	// Si octree doit etre subdivisee
	else if(inside.size() != 0)
	{
		// lowLefFron
		simplifySubOct(MED_X, MED_Y, MED_Z, min_x, min_y, min_z, n, newPositions, newNormals, oldToNew);
		// uppLefFron
		simplifySubOct(MED_X, max_y, MED_Z,	min_x, MED_Y, min_z, n, newPositions, newNormals, oldToNew);
		// uppRigFron
		simplifySubOct(max_x, max_y, MED_Z, MED_X, MED_Y, min_z, n, newPositions, newNormals, oldToNew);
		// lowRigFron
		simplifySubOct(max_x, MED_Y, MED_Z, MED_X, min_y, min_z, n, newPositions, newNormals, oldToNew);
		
		// lowLefBack
		simplifySubOct(MED_X, MED_Y, max_z, min_x, min_y, MED_Z, n, newPositions, newNormals, oldToNew);
		// uppLefBack
		simplifySubOct(MED_X, max_y, max_z,	min_x, MED_Y, MED_Z, n, newPositions, newNormals, oldToNew);
		// uppRigBack
		simplifySubOct(max_x, max_y, max_z, MED_X, MED_Y, MED_Z, n, newPositions, newNormals, oldToNew);
		// lowRigBack
		simplifySubOct(max_x, MED_Y, max_z, MED_X, min_y, MED_Z, n, newPositions, newNormals, oldToNew);
	}
}

// Simplifier la mesh utilisant un octree
void Mesh::simplifyAdaptive(unsigned int n)
{
	if(m_positions.size() < n)
		return;

	std::vector<unsigned int> oldToNew;
	oldToNew.resize(m_positions.size());
	float max_x, max_y, max_z, min_x, min_y, min_z;
	getMaxMin(max_x, max_y, max_z, min_x, min_y, min_z);
	
	max_x *= 1.01;
	min_x *= 1.01;
	
	max_y *= 1.01;
	min_y *= 1.01;
	
	max_z *= 1.01;
	min_z *= 1.01;
	
	std::vector<Vec3f> newPositions, newNormals;
	
	// Remplit newNormals, newPositions et oldToNew
	simplifySubOct(max_x, max_y, max_z, min_x, min_y, min_z, n, newPositions, newNormals, oldToNew);
	
	// Apres calcul des nouveaux points et indices au dessus, echange indices anciens par les nouveaux
	for(unsigned int n = 0; n < m_triangles.size(); n++) {
		for(unsigned int corner = 0; corner < 3; corner++) {
			m_triangles[n][corner] = oldToNew[ m_triangles[n][corner] ];
		}
	}
	
	// efface triangles degeneres
	for(unsigned int n = 0; n < m_triangles.size(); n++) {		
		if(	m_triangles[n][0] == m_triangles[n][1] ||
			m_triangles[n][0] == m_triangles[n][2] ||
			m_triangles[n][1] == m_triangles[n][2])

			m_triangles.erase(m_triangles.begin() + n);
	}
	
	m_positions.clear();
	m_normals.clear();
	
	m_positions = newPositions;
	m_normals = newNormals;
	
}

//-----------------------------------------------------------------------------
//-------------------------------simplify--------------------------------------
//-----------------------------------------------------------------------------

void Mesh::simplify(unsigned int resolution)
{

	float max_x, max_y, max_z, min_x, min_y, min_z;
	std::vector<Vec3f> newPositions, newNormals;
	std::vector<unsigned int> oldToNew;
	newPositions.resize(resolution*resolution*resolution);
	newNormals.resize(resolution*resolution*resolution);
	oldToNew.resize(m_positions.size());
	getMaxMin(max_x, max_y, max_z, min_x, min_y, min_z);
	
	max_x *= 1.01;
	min_x *= 1.01;
	
	max_y *= 1.01;
	min_y *= 1.01;
	
	max_z *= 1.01;
	min_z *= 1.01;
	
	// cree cubes
	std::vector< std::vector< std::vector<Cube> > > cubes(resolution, std::vector< std::vector<Cube> > (resolution, std::vector<Cube>(resolution)));
	
	for(unsigned int i = 0; i < resolution; i++)
	{
		for(unsigned int j = 0; j < resolution; j++)
		{
			for(unsigned int k = 0; k < resolution; k++)
			{
				cubes[i][j][k].element = 0;
				
				cubes[i][j][k].max_x = max_x - ((float)i/(float)resolution) * (max_x - min_x);
				cubes[i][j][k].min_x = cubes[i][j][k].max_x - ((float)1/(float)resolution) * (max_x - min_x);
				
				cubes[i][j][k].max_y = max_y - ((float)j/(float)resolution) * (max_y - min_y);
				cubes[i][j][k].min_y = cubes[i][j][k].max_y - ((float)1/(float)resolution) * (max_y - min_y);
				
				cubes[i][j][k].max_z = max_z - ((float)k/(float)resolution) * (max_z - min_z);
				cubes[i][j][k].min_z = cubes[i][j][k].max_z - ((float)1/(float)resolution) * (max_z - min_z);
			}
		}
	}
	

	// calcule position, normale et elements de chaque cube
	for(unsigned int n = 0; n < m_positions.size(); n++) {
	
		for(unsigned int i = 0; i < resolution; i++)
		{
			for(unsigned int j = 0; j < resolution; j++)
			{
				for(unsigned int k = 0; k < resolution; k++)
				{
				
					if(	m_positions[n][0] < cubes[i][j][k].max_x &&
						m_positions[n][0] > cubes[i][j][k].min_x &&
						
						m_positions[n][1] < cubes[i][j][k].max_y &&
						m_positions[n][1] > cubes[i][j][k].min_y &&
						
						m_positions[n][2] < cubes[i][j][k].max_z &&
						m_positions[n][2] > cubes[i][j][k].min_z)
					{
						newPositions[i*resolution*resolution + j*resolution + k] += m_positions[n];
						newNormals[i*resolution*resolution + j*resolution + k] += m_normals[n];
						cubes[i][j][k].element++;
						
						// vecteur avec nouveaux indices en fonction des anciens
						oldToNew[n] = i*resolution*resolution + j*resolution + k;
						break;
					}
				}
			}
		}
	}
	
	// divise par le nombre d'elements et normalise
	for(unsigned int i = 0; i < resolution; i++)
	{
		for(unsigned int j = 0; j < resolution; j++)
		{
			for(unsigned int k = 0; k < resolution; k++)
			{
				if(cubes[i][j][k].element != 0)
					newPositions[i*resolution*resolution + j*resolution + k] /= (float) cubes[i][j][k].element;
					
				normalize(newNormals[i*resolution*resolution + j*resolution + k]);
			}
		}
	}

	// echange anciens indices par les nouveaux
	for(unsigned int n = 0; n < m_triangles.size(); n++) {
		for(unsigned int corner = 0; corner < 3; corner++) {
			m_triangles[n][corner] = oldToNew[ m_triangles[n][corner] ];
		}
	}
	
	// efface triangles degeneres
	for(unsigned int n = 0; n < m_triangles.size(); n++) {
		if(	m_triangles[n][0] == m_triangles[n][1] ||
			m_triangles[n][0] == m_triangles[n][2] ||
			m_triangles[n][1] == m_triangles[n][2])
				m_triangles.erase(m_triangles.begin() + n);
	}	
	
	m_positions = newPositions;
	m_normals = newNormals;
	
}

//-----------------------------------------------------------------------------
//-------------------------------subdivide-------------------------------------
//-----------------------------------------------------------------------------

void Mesh::subdivide()
{
	unsigned int posSize = m_positions.size();
	
	// remplit map avec le milieu de chaque arrete
	fillMiddle();
	// fait le moyennage de chaque ancien point
	triAveraging(posSize);
	
	// distribue les nouveaux points et cree les nouveaux triangles
	subdivideTri();
	
	// fait nouveau calcul pour les plus proches voisins
	recomputeNN();
	// reclacule les normales
	recomputeNormals();
}

void Mesh::recomputeNN()
{
	m_nneighbours.clear();
	m_nneighbours.resize(m_positions.size());

	// garde tous les triangles qui contiennent un point	
	for(unsigned int tri = 0; tri < m_triangles.size(); tri++) {
		for(unsigned int corner = 0; corner < 3 ; corner++)
		{
			m_nneighbours[m_triangles[tri][corner]].push_back(m_triangles[tri]);	
		}
	}
}

void Mesh::subdivideTri()
{

	std::vector<Triangle> newInd;
	
	// pour chaque triangle, garde tous les points au milieu de ses arretes
	// et garde les triangles necessaires
	for(unsigned int tri = 0; tri < m_triangles.size(); tri++) {
		unsigned int mid1 = middle_points[Edge(m_triangles[tri][0], m_triangles[tri][1])];
		unsigned int mid2 = middle_points[Edge(m_triangles[tri][1], m_triangles[tri][2])];
		unsigned int mid3 = middle_points[Edge(m_triangles[tri][0], m_triangles[tri][2])];
		
		newInd.push_back( Triangle(m_triangles[tri][0] , mid1 , mid3 ) );
		newInd.push_back( Triangle(m_triangles[tri][1] , mid2 , mid1 ) );
		newInd.push_back( Triangle(m_triangles[tri][2] , mid3 , mid2 ) );
		newInd.push_back( Triangle( mid1, mid2, mid3 ) );
		
	}
	m_triangles = newInd;
	
}

void Mesh::triAveraging(unsigned int limit)
{
	for(unsigned int p = 0; p < limit; p++) {
		Vec3f newPos = Vec3f(0.0,0.0,0.0);
		float alpha = ( 40 - pow( 3 + 2 * cos( 2 * M_PI / (float) m_positions.size() ) , 2 ) ) / 64.0;
		
		// prends points voisins au point en question
		std::set<unsigned int> voisins = getVoisins(m_nneighbours[p], p);
		
		// fait la moyenne des positions avec le poids alpha
		std::set<unsigned int>::iterator it;
		for (it = voisins.begin(); it != voisins.end(); it++)
		{
				newPos += m_positions[*it] * alpha / (float) voisins.size();
		}
	
		// recentre le point
		newPos += m_positions[p] * (1 - alpha);
	
		// Correction pour arretes extraoridnaires
		newPos = m_positions[p] + (newPos - m_positions[p]) * (5.0 / 3.0 - 8.0 / 3.0 * pow( 3.0/8.0 + 1.0/4.0 * cos(2 * 3.1416 / voisins.size()) , 2));
		//newPos = m_positions[p] + (newPos - m_positions[p]) * (8.0 / (voisins.size() + 2.0));
	
		m_positions[p] = newPos;
	}
}

void Mesh::fillMiddle()
{
	unsigned int limit = m_positions.size();
	for (unsigned int i = 0; i < limit; i++)
	{
		std::set<unsigned int> voisins = getVoisins(m_nneighbours[i], i);
	
		// pour chaque voisin du point, garde milieu des deux si pas fait avant
		std::set<unsigned int>::iterator it;
		for (it = voisins.begin(); it != voisins.end(); it++)
		{
			if(middle_points.find(Edge(i,*it)) == middle_points.end())
			{
				Vec3f middle = Vec3f(0.0,0.0,0.0);
		
				// prends triangles qui appartiennent aux deux ensembles
				// donc les triangles qui contiennent les deux points a la fois
				std::vector<Triangle> inter = intersectTriVect(m_nneighbours[*it], m_nneighbours[i]);
				
				// calcule contribution selon poids donnes
				std::vector<Triangle>::iterator itInter;
				for (itInter = inter.begin(); itInter != inter.end(); itInter++)
				{
					for(unsigned int j = 0; j < 3; j++)
					{
						if((*itInter)[j] == *it || (*itInter)[j] == i)
							middle += m_positions[(*itInter)[j]] * 3 / (inter.size()*8.0);
						else
							middle += m_positions[(*itInter)[j]] * 1 / (inter.size()*4.0);
					}
				}
		
				m_positions.push_back( middle );			
				middle_points[Edge(i,*it)] = m_positions.size() - 1;
			}
		}
	}
}

// donne intersection de deux ensembles de triangles
std::vector<Triangle> Mesh::intersectTriVect(std::vector<Triangle> v1, std::vector<Triangle> v2)
{
	std::vector<Triangle> inter;
	
	std::vector<Triangle>::iterator it1;
	for (it1 = v1.begin(); it1 != v1.end(); it1++)
	{
		std::vector<Triangle>::iterator it2;
		for (it2 = v2.begin(); it2 != v2.end(); it2++)
		{
			if(((*it1)[0] == (*it2)[0]) && ((*it1)[1] == (*it2)[1]) && ((*it1)[2] == (*it2)[2]))
				inter.push_back(*it1);
		}
    }
    
    return inter;
}

// donne points voisins en sachant les triangles qui contiennent le point en question
std::set<unsigned int> Mesh::getVoisins(std::vector<Triangle> tri, unsigned int point) {
	std::set<unsigned int> voisins;
	
	std::vector<Triangle>::iterator it;
	for (it = tri.begin(); it != tri.end(); it++)

	{
		Triangle t = (*it);

		if(t[0] == point)
		{
			voisins.insert(t[1]);
			voisins.insert(t[2]);
		}
		else if(t[1] == point)
		{
			voisins.insert(t[0]);
			voisins.insert(t[2]);
		}
		else
		{
			voisins.insert(t[1]);
			voisins.insert(t[0]);
		}
	}

	return voisins;
}

// calcule barycentre de tous les points
Vec3f Mesh::calculerBarycentre(std::set<unsigned int> points) {
	float sumX = 0;
	float sumY = 0;
	float sumZ = 0;
	
	std::set<unsigned int>::iterator it;
	for (it = points.begin(); it != points.end(); it++)
	{
		sumX += m_positions[*it][0];
		sumY += m_positions[*it][1];
		sumZ += m_positions[*it][2];
	}

	return Vec3f(sumX / points.size(), sumY / points.size(), sumZ / points.size());
}

// calcule barycentre d'un triangle
Vec3f Mesh::calculerBarycentreTri(Triangle tri) {
	float sumX = 0;
	float sumY = 0;
	float sumZ = 0;

	for (unsigned int i = 0; i <= 3; i++)
	{
		sumX += m_positions[tri[i]][0];
		sumY += m_positions[tri[i]][1];
		sumZ += m_positions[tri[i]][2];
	}

	return Vec3f(sumX / 3.0, sumY / 3.0, sumZ / 3.0);
}

// calcule barycentre geometrique pour laplacien geometrique
Vec3f Mesh::calculerBarycentreGeom(unsigned int point, std::set<unsigned int> points, std::vector<Triangle> triangles) {
	float sumX = 0.0f;
	float sumY = 0.0f;
	float sumZ = 0.0f;
	float sumCot = 0.0f;
	
	std::set<unsigned int>::iterator it;
	for (it = points.begin(); it != points.end(); it++)
	{
		std::vector<Triangle> triContain = containPoint(*it, triangles);
		float cot = 0.0f;
		
		// calcul du poids avec les cotangentes
		std::vector<Triangle>::iterator triIt;
		for (triIt = triContain.begin(); triIt != triContain.end(); ++triIt)
			if(!std::isnan(cotan(point, *it, *triIt)) && cotan(point, *it, *triIt) < pow(10,3))
				cot += cotan(point, *it, *triIt);
			
		sumX += cot*m_positions[*it][0];
		sumY += cot*m_positions[*it][1];
		sumZ += cot*m_positions[*it][2];
		
		sumCot += cot;
	}

	return Vec3f(sumX , sumY, sumZ) / sumCot;
}

// retourne triangles qui contiennent le point
std::vector<Triangle> Mesh::containPoint(unsigned int point, std::vector<Triangle> triangles)
{

	std::vector<Triangle> result;
	
	std::vector<Triangle>::iterator it;
	for (it = triangles.begin(); it != triangles.end(); it++)
		if((*it)[0] == point || (*it)[1] == point || (*it)[2] == point)
			result.push_back(*it);
			
	return result;
}

// calcule cotangente
float Mesh::cotan(unsigned int point1, unsigned int point2, Triangle triangle)
{
	Vec3f e1, e2;

	if(point1 != triangle[0] && point2 != triangle[0])
	{
		e1 = m_positions[triangle[0]] - m_positions[triangle[1]];
		e2 = m_positions[triangle[0]] - m_positions[triangle[2]];
	}
	else if(point1 != triangle[1] && point2 != triangle[1])
	{
		e1 = m_positions[triangle[1]] - m_positions[triangle[0]];
		e2 = m_positions[triangle[1]] - m_positions[triangle[2]];
	}
	else
	{
		e1 = m_positions[triangle[2]] - m_positions[triangle[1]];
		e2 = m_positions[triangle[2]] - m_positions[triangle[0]];	
	}
	
	return dot(-e1, e2) / sqrt(e1.squaredLength() * e2.squaredLength() - pow(dot(e1, e2) , 2));
	
}

//-----------------------------------------------------------------------------
//---------------------------------addFloor------------------------------------
//-----------------------------------------------------------------------------

// Ajoute un des deux types de sol
void Mesh::addFloor (unsigned int type)	{

	float minY, temp;
	getMaxMin(temp, temp, temp, temp, minY, temp);
	unsigned int offset = m_positions.size();

	if(type == 1) {
		for(unsigned int i = 0 ; i < 20 ; i++)	{
			for(unsigned int j = 0 ; j < 20 ; j++)	{
					m_positions.push_back( Vec3f( (j * 0.2) - 2 , minY , 2 - (i * 0.2) ) );
			}		
		}
		for(unsigned int i = 0 ; i < 20 - 1 ; i++)	{
			for(unsigned int j = 0 ; j < 20 - 1 ; j++)	{
					m_triangles.push_back( Triangle(offset + i * 20 + j, offset + i * 20 + j + 1, offset + (i + 1) * 20 + j) );
					m_triangles.push_back( Triangle(offset + i * 20 + j + 1, offset + (i + 1) * 20 + j + 1, offset + (i + 1) * 20 + j) );
			}		
		}
		
	}
	else if(type == 2) {
		for(unsigned int i = 0 ; i < 200 ; i++)	{
			for(unsigned int j = 0 ; j < 200 ; j++)	{
					m_positions.push_back( Vec3f( (j * 0.02) - 2 , minY , 2 - (i * 0.02) ) );
			}		
		}
		for(unsigned int i = 0 ; i < 200 - 1 ; i++)	{
			for(unsigned int j = 0 ; j < 200 - 1 ; j++)	{
					m_triangles.push_back( Triangle(offset + i * 200 + j, offset + i * 200 + j + 1, offset + (i + 1) * 200 + j) );
					m_triangles.push_back( Triangle(offset + i * 200 + j + 1, offset + (i + 1) * 200 + j + 1, offset + (i + 1) * 200 + j) );
			}		
		}
	}	
	
	recomputeNormals ();
}

// Enleve le sol
void Mesh::removeFloor (unsigned int type)	{
	if(type == 1) {
		m_positions.resize(m_positions.size() - 20 * 20);
		m_triangles.resize(m_triangles.size() - 2 * (20 - 1) * (20 - 1));
	}
	else if(type == 2) {
		m_positions.resize(m_positions.size() - 200 * 200);
		m_triangles.resize(m_triangles.size() - 2 * (200 - 1) * (200 - 1));
	}
}



void Mesh::centerAndScaleToUnit () {
    Vec3f c;
    for  (unsigned int i = 0; i < m_positions.size (); i++)
        c += m_positions[i];
    c /= m_positions.size ();
    float maxD = dist (m_positions[0], c);
    for (unsigned int i = 0; i < m_positions.size (); i++){
        float m = dist (m_positions[i], c);
        if (m > maxD)
            maxD = m;
    }
    for  (unsigned int i = 0; i < m_positions.size (); i++)
        m_positions[i] = (m_positions[i] - c) / maxD;
}
