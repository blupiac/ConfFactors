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
#include <fstream>
#include <cstdlib>
#include <string>
#include <algorithm>
#include <cmath>
#include <assert.h>

//#define DEBUG
#define TIMER

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
}

// Utilise pour indexes qui representent des coordonees de position
enum Dim { x, y, z };

//-----------------------------------------------------------------------------
//--------------------------------loadOFF--------------------------------------
//-----------------------------------------------------------------------------

void Mesh::loadOFF (const std::string & filename) {

	#ifdef TIMER
	tInit = std::clock();
	std::cout << "Computation begins." << std::endl;
	#endif

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
    m_areas.resize (sizeT);
    m_sortedAreasIdx.resize (sizeT);
    m_confFacts.resize (sizeV);
    m_gaussDiff.resize (sizeV);
    #ifdef DEBUG
    m_gausscurv.resize (sizeV);
    #endif
    m_nneighbours.resize(sizeV);
    totalArea = 0; totalCurv = 0; totalConf = 0;
    for (unsigned int i = 0; i < sizeV; i++)
    {
        in >> m_positions[i];
    }
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
            m_nneighbours[m_triangles[i][j]].push_back(i);
        }
        
        float area = getArea(m_triangles[i]);

        m_areas[i] = area;
		totalArea += area;
    }
    in.close ();

    // sort area idx vector
    std::size_t n(0);
    std::generate(std::begin(m_sortedAreasIdx), std::end(m_sortedAreasIdx), [&]{ return n++; });

    std::sort(  std::begin(m_sortedAreasIdx), 
                std::end(m_sortedAreasIdx),
                [&](int i1, int i2) { return m_areas[i1] < m_areas[i2]; } );

	// fill confFactor vector
	#ifdef TIMER
	std::cout << "Initialization ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif
    calculateConfFact();
    #ifdef TIMER
	std::cout << "ConfFact calculation ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif
    calculateSignature();
    #ifdef TIMER
	std::cout << "Signature generation ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

    centerAndScaleToUnit ();
    recomputeNormals ();
}

//-----------------------------------------------------------------------------
//--------------------------calculateConfFact----------------------------------
//-----------------------------------------------------------------------------

void Mesh::calculateConfFact () 
{
	minConf = 1e10;
	maxConf = -1e10;

	#ifdef DEBUG
	minGauss = 1e10;
	maxGauss = -1e10;
	#endif

	// calculate gaussian curv
	for (unsigned int i = 0; i < m_positions.size (); i++)
    {
    	// gaussDiff = target - gaussian
    	// cannot calculate target without total gauss curv
    	// so we start filling gaussDiff now to avoid recalculating
    	// the gaussDiff below
    	m_gaussDiff[i] = -getGaussCurv(i);
        totalCurv -= m_gaussDiff[i];
    }

    #ifdef TIMER
	std::cout << "Gauss curv ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

    // calculate target curv, and join the diff in a vector
	for (unsigned int i = 0; i < m_positions.size (); i++) 
	{
		float targCurv = getTargetCurv(i);

		// see explanation just above
		m_gaussDiff[i] += targCurv;
	}

	#ifdef TIMER
	std::cout << "Target curv ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	// laplacian matrix
	Eigen::SparseMatrix<float> lapMatrix;
	// rough estimation of non-zero elements: each vertex will have 6 neighbours plus itself
	lapMatrix.reserve(m_positions.size() * 7);
	lapMatrix = getLapMatrix();

	#ifdef TIMER
	std::cout << "Lap matrix fill ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	m_confFacts = solveConfFactor(m_gaussDiff, lapMatrix);

	#ifdef TIMER
	std::cout << "ConfFactor equation solving ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	auto result = std::minmax_element(m_confFacts.begin(), m_confFacts.end());

	minConf = m_confFacts[result.first - m_confFacts.begin()];
	maxConf = m_confFacts[result.second - m_confFacts.begin()];
}

float Mesh::normalizeConf(unsigned int confIdx)
{
	//std::cout << "this conf: " << m_confFacts[confIdx] << " min conf: " << minConf << " max conf: " << maxConf << std::endl;
	//std::cout << "normalized: " << (m_confFacts[confIdx] - minConf) / (maxConf - minConf) << std::endl;

	return (m_confFacts[confIdx] - minConf) / (maxConf - minConf);
}

#ifdef DEBUG

float Mesh::normalizeGausscurv(unsigned int gaussIdx)
{
	return (m_gausscurv[gaussIdx] - minGauss) / (maxGauss - minGauss);
}

#endif

float Mesh::getGaussCurv(unsigned int pointIdx)
{
	std::vector<unsigned int> neighbours = m_nneighbours[pointIdx];
	float totalAngle = 0;
	float curv;

	for(unsigned int i = 0; i < neighbours.size(); i++)
	{
		totalAngle += getAngle(m_triangles[neighbours[i]], pointIdx);
	}

	// interior
	if(!isBorder(pointIdx))
	{
		curv = 2 * M_PI - totalAngle;
	}
	// edge
	else
	{
		curv = M_PI - totalAngle;
	}

	#ifdef DEBUG
	m_gausscurv[pointIdx] = curv;

	if(minGauss > curv)
	{
		minGauss = curv; 
	}
	if(maxGauss < curv)
	{
		maxGauss = curv;
	}
	#endif

	return curv;
}

bool Mesh::isBorder(unsigned int ptIdx)
{
	std::set<unsigned int> neighPoint = getVoisins(m_nneighbours[ptIdx], ptIdx);

	std::set<unsigned int>::iterator ptIt;
	for (ptIt = neighPoint.begin(); ptIt != neighPoint.end(); ++ptIt)
	{
		unsigned int occur = 0;
		std::vector<unsigned int>::iterator triIt;
		for (triIt = m_nneighbours[ptIdx].begin(); triIt != m_nneighbours[ptIdx].end(); ++triIt)
		{
			Triangle t = m_triangles[*triIt];
			if(*ptIt == t[0] || *ptIt == t[1] || *ptIt == t[2])
			{
				occur++;
			}
		}

		if(occur < 2)
		{
			return true; 
		}
	}

	return false;
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

	return acos(res);
}

float Mesh::getTargetCurv(unsigned int i)
{
	std::vector<unsigned int> tris = m_nneighbours[i];
	float area = 0;

	for(unsigned int i = 0; i < tris.size(); i++)
	{
		area += getArea(m_triangles[tris[i]]) / 3;
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

// TODO: change all float to double

// based, but far from copied from:
// https://github.com/Sunwinds/3D-Geometric-Features/blob/master/ConformalFactor/CalConformalFactor.h
Eigen::SparseMatrix<float> Mesh::getLapMatrix()
{
	// triples to be inserted later on sparse matrix
	std::vector<Eigen::Triplet<float>> triplets;
	// another rough estimation: each vertex has 6 neighbours
	triplets.reserve(m_positions.size() * 6 * 4);

	for (unsigned int i = 0; i < m_positions.size (); i++) 
	{
		// get 1-neighborhood points
		std::set<unsigned int> neighPoint = getVoisins(m_nneighbours[i], i);

		std::set<unsigned int>::iterator ptIt;
		for (ptIt = neighPoint.begin(); ptIt != neighPoint.end(); ++ptIt)
		{
			float cots = 0;

			// get 2 triangles that contain main point plus the one of this iteration
			std::vector<unsigned int> triContain = containPoint(*ptIt, m_nneighbours[i]);
			
			if(triContain.size() == 2)
			{
				// calculate cot on both angles wanted
				float cotA = cotan(i, *ptIt, triContain[0]);
				float cotB = cotan(i, *ptIt, triContain[1]);
				float cotSum = cotA + cotB;
				if(!std::isnan(cotSum) && (cotSum) < pow(10,5))
				{
					cots += cotSum;
				}
				else
				{
					std::cerr << "Point with NaN or impossible cotan: " << i << std::endl;
				}
			}

			if(abs(cots) == 0)
			{
				continue;
			}

			// divide by 2 because each cotan gets added twice 
			// TODO: maybe there's a way to only do this once?
			cots /= 2.0;

			// add cots to both matrix elements
			triplets.push_back(Eigen::Triplet<float>(i, *ptIt, -cots));
			triplets.push_back(Eigen::Triplet<float>(*ptIt, i, -cots));
			// diag elems
			triplets.push_back(Eigen::Triplet<float>(i, i, cots));
			triplets.push_back(Eigen::Triplet<float>(*ptIt, *ptIt, cots));
		}

		if(neighPoint.size() == 0)
		{
			std::cerr << "Point " << i << " has no neighbours." << std::endl;
		}
	}

	// initializing sparse matrix
	Eigen::SparseMatrix<float> lap(m_positions.size(), m_positions.size());
	// inserts all elements in vector in matrix, if duplicated, they are added
	lap.setFromTriplets(triplets.begin(), triplets.end());

	return lap;
}

// https://scicomp.stackexchange.com/questions/21343/solving-linear-equations-using-eigen
std::vector<float> Mesh::solveConfFactor(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian)
{
	Eigen::SparseLU<Eigen::SparseMatrix<float> > solverLap;
	Eigen::VectorXf gaussDiffVect = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(gaussDiff.data(), gaussDiff.size());
	Eigen::VectorXf confFact;

	laplacian.makeCompressed();
	solverLap.analyzePattern(laplacian);
	solverLap.factorize(laplacian);

	if(solverLap.info() != Eigen::Success)
	{
		std::cerr << "Failure to solve confFact linear equations: " << solverLap.lastErrorMessage() << std::endl;
		exit(1);
	}

	confFact = solverLap.solve(gaussDiffVect);
	std::vector<float> res(&confFact[0], confFact.data()+confFact.cols()*confFact.rows());
	return res;
}

// TODO: multi threading
// recommended solver for 3D poisson eq. : https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
std::vector<float> Mesh::solveConfFactorCG(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian)
{
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> cgLap;
	Eigen::VectorXf gaussDiffVect = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(gaussDiff.data(), gaussDiff.size());
	Eigen::VectorXf confFact;

	cgLap.compute(laplacian);

	//std::cout << "#iterations:     " << cgLap.iterations() << std::endl;
	//std::cout << "estimated error: " << cgLap.error()      << std::endl;
	
	confFact = cgLap.solve(gaussDiffVect);
	std::vector<float> res(&confFact[0], confFact.data()+confFact.cols()*confFact.rows());
	return res;
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
		float randConfFactor = interpConfFactor(randPoint, triangleIdx);

		incrSignature(randConfFactor, minConf, maxConf, -99, 100);
	}

	#ifdef TIMER
	std::cout << "Point generation ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	for (unsigned int i = 0; i < m_positions.size (); i++) 
	{
		incrSignature(m_confFacts[i], minConf, maxConf, -99, 100);
	}

	#ifdef TIMER
	std::cout << "Bin filling ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	printSignature(signature, m_positions.size());
}

void Mesh::initializeSignature(int min, int max)
{
	signature.resize (max - min + 1);

	for(unsigned int i = 0; i < signature.size(); i ++)
	{
		Bin b;
		b.occur = 0.0f;
		b.idx = min + i;
		signature[i] = b;
	}
}

void Mesh::printSignature(std::vector<Bin> sig, unsigned int totalItems)
{
	ofstream outfile;
    outfile.open ("output/confFact.csv");

    outfile << "binNumber, occurence\n";

	for(unsigned int i = 0; i < signature.size(); i ++)
	{
		float occurPercent = sig[i].occur / (float) totalItems;
		outfile << sig[i].idx << ", " << occurPercent << "\n";
		//std::cout << "Bin " << sig[i].idx << ": " << occurPercent << " occurences: " << sig[i].occur << std::endl;
	}

	outfile.close();
}

int Mesh::getRandTri()
{
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float randArea = r * totalArea;
	unsigned int i = 0;

	while((randArea - m_areas[m_sortedAreasIdx[i]])>= 0)
	{
		randArea -= m_areas[m_sortedAreasIdx[i]];
		i++;
	}
	return i;
}

// http://www.cs.princeton.edu/~funk/tog02.pdf
// p8, equation (1)
Vec3f Mesh::getRandPoint(Triangle tri)
{
	float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	return Vec3f((1 - sqrt(r1)) * m_positions[tri[0]] + 
				sqrt(r1) * (1 - r2) * m_positions[tri[1]] + 
				sqrt(r1) * r2 * m_positions[tri[2]]);
}

// https://stackoverflow.com/questions/18755251/linear-interpolation-of-three-3d-points-in-3d-space
// section 3.4, use barycentric coordinates to find weight associated to each vertex
float Mesh::interpConfFactor(Vec3f point, unsigned int triIdx)
{
	unsigned int a = m_triangles[triIdx][0];
	unsigned int b = m_triangles[triIdx][1];
	unsigned int c = m_triangles[triIdx][2];

	Vec3f v0 = m_positions[b] - m_positions[a];
	Vec3f v1 = m_positions[c] - m_positions[a];
	Vec3f v2 = point - m_positions[a];

	float d00 = dot(v0, v0);
	float d01 = dot(v0, v1);
	float d11 = dot(v1, v1);
	float d20 = dot(v2, v0);
	float d21 = dot(v2, v1);

	float denom = d00 * d11 - d01 * d01;

	float v = (d11 * d20 - d01 * d21) / denom;
	float w = (d00 * d21 - d01 * d20) / denom;
	float u = 1.0 - v - w;

	return v * m_confFacts[a] + w * m_confFacts[b] + u * m_confFacts[c];
}

void Mesh::incrSignature(float confFact, float min, float max, int binMin, int binMax)
{
	float normIdx = (confFact - min) / (max - min);
	int equivBin = normIdx * (binMax - binMin);

	//std::cout << "confFact: " << confFact << " min: " << min << " max: " << max << " normIdx: " << normIdx << " equivBin: " << equivBin << std::endl;

	signature[equivBin].occur += 1;
}

//-----------------------------------------------------------------------------
//-------------------------------helpers---------------------------------------
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


// donne points voisins en sachant les triangles qui contiennent le point en question
std::set<unsigned int> Mesh::getVoisins(std::vector<unsigned int> tri, unsigned int point) {
	std::set<unsigned int> voisins;
	
	std::vector<unsigned int>::iterator it;
	for (it = tri.begin(); it != tri.end(); it++)

	{
		Triangle t = m_triangles[*it];

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

// retourne triangles qui contiennent le point
std::vector<unsigned int> Mesh::containPoint(unsigned int point, std::vector<unsigned int> triangles)
{

	std::vector<unsigned int> result;
	
	std::vector<unsigned int>::iterator it;
	for (it = triangles.begin(); it != triangles.end(); it++)
	{
		Triangle t = m_triangles[*it];
		if(t[0] == point || t[1] == point || t[2] == point)
			result.push_back(*it);
	}	
	return result;
}

// calcule cotangente
float Mesh::cotan(unsigned int point1, unsigned int point2, unsigned int triangle)
{
	Vec3f e1, e2;

	Triangle t = m_triangles[triangle];

	if(point1 != t[0] && point2 != t[0])
	{
		e1 = m_positions[t[0]] - m_positions[t[1]];
		e2 = m_positions[t[0]] - m_positions[t[2]];
	}
	else if(point1 != t[1] && point2 != t[1])
	{
		e1 = m_positions[t[1]] - m_positions[t[0]];
		e2 = m_positions[t[1]] - m_positions[t[2]];
	}
	else
	{
		e1 = m_positions[t[2]] - m_positions[t[1]];
		e2 = m_positions[t[2]] - m_positions[t[0]];	
	}
	
	return dot(-e1, e2) / sqrt(e1.squaredLength() * e2.squaredLength() - pow(dot(e1, e2) , 2));
	
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
