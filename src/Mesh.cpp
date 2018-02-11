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
#include <assert.h>

//#define DEBUG
#define TIMER

#define MED_X min_x + (max_x - min_x)/2
#define MED_Y min_y + (max_y - min_y)/2
#define MED_Z min_z + (max_z - min_z)/2

#define PRUNEFACT 0.01
#define NOISEFACT 0.0

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
    m_triangles.clear ();
    m_nneighbours.clear();
}

// Utilise pour indexes qui representent des coordonees de position
enum Dim { x, y, z };

//-----------------------------------------------------------------------------
//--------------------------------loadOFF--------------------------------------
//-----------------------------------------------------------------------------


/**
    Parses and does the pre-processing of an off file

    @param filename filepath of an off file to be loaded
*/
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
    totalArea = 0; totalCurv = 0;
    for (unsigned int i = 0; i < sizeV; i++)
    {
        in >> m_positions[i];
    }
    if(NOISEFACT != 0) {
    	makeNoise(NOISEFACT);
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

        // cumulative area array to make binary search possible later
        if(i == 0)
        {
        	m_areas[i] = area;
        }
        else
        {
        	m_areas[i] = area + m_areas[i - 1];
        }
		totalArea += area;
    }
    in.close ();

    // sort area idx vector
    std::size_t n(0);
    std::generate(std::begin(m_sortedAreasIdx), std::end(m_sortedAreasIdx), [&]{ return n++; });

    std::sort(  std::begin(m_sortedAreasIdx), 
                std::end(m_sortedAreasIdx),
                [&](int i1, int i2) { return m_areas[i1] < m_areas[i2]; } );

    std::sort(  std::begin(m_areas), 
                std::end(m_areas),
                [&](float f1, float f2) { return f1 < f2; } );

	// fill confFactor vector
	#ifdef TIMER
	std::cout << "Initialization ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif
    calculateConfFact();
    #ifdef TIMER
	std::cout << "ConfFact calculation ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif
	if(PRUNEFACT == 0) {
		calculateSignature();
	}
	else {
		calculatePruneSignature(PRUNEFACT);
	}
    #ifdef TIMER
	std::cout << "Signature generation ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

    centerAndScaleToUnit ();
}

/**
    Introduces noise into the mesh

    @param noiseCoef noise coefficient, values between
    0.01-0.1 give acceptable results, read readme for more
    information
*/
void Mesh::makeNoise(float noiseCoef) {
	std::srand(std::time(nullptr));
	for (unsigned int i = 0; i < m_positions.size (); i++) {
		float randCoef = ((double) std::rand() / (RAND_MAX)) + 1;
		m_positions[i] += m_positions[i] * randCoef * NOISEFACT;
    }
}

//-----------------------------------------------------------------------------
//--------------------------calculateConfFact----------------------------------
//-----------------------------------------------------------------------------

/**
    Fills the m_confFact vector with the confFactors of the mesh
*/
void Mesh::calculateConfFact () {
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

/**
    Returns the given conformal factor after being normalized 

    @param confIdx index of conformal factor to be returned
    @return normalized conformal factor
*/
float Mesh::normalizeConf(unsigned int confIdx) {
	return (m_confFacts[confIdx] - minConf) / (maxConf - minConf);
}

#ifdef DEBUG
/**
    Returns the given gaussian curvature after being normalized 

    @param gaussIdx index of conformal factor to be returned
    @return normalized gausian curvature
*/
float Mesh::normalizeGausscurv(unsigned int gaussIdx) {
	return (m_gausscurv[gaussIdx] - minGauss) / (maxGauss - minGauss);
}

#endif

/**
    Returns gaussian curvature of a vertex

    @param i vertex index
    @return gaussian curvature of vertex
*/
float Mesh::getGaussCurv(unsigned int pointIdx) {
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

/**
    Returns whether the vertex is inside a mesh or in its border

    @param ptIdx vertex index
    @return boolean saying whether vertex is in border or not
*/
bool Mesh::isBorder(unsigned int ptIdx) {
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

/**
    Returns angle of a vertex in a triangle

    @param tri triangle inside of which the vertex is in
    @param pointIdx vertex index
    @return angle of vertex
*/
float Mesh::getAngle(Triangle tri, int pointIdx) {
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

/**
    Returns target curvature of a vertex

    @param i vertex index
    @return target curvature of vertex
*/
float Mesh::getTargetCurv(unsigned int i) {
	std::vector<unsigned int> tris = m_nneighbours[i];
	float area = 0;

	for(unsigned int i = 0; i < tris.size(); i++)
	{
		area += getArea(m_triangles[tris[i]]) / 3;
	}

	return totalCurv * area / totalArea;
}

/**
    Returns area of a triangle

    @param i triangle index
    @return area of triangle
*/
float Mesh::getArea(unsigned int i) {
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

/**
    Returns area of a triangle

    @param tri triangle
    @return area of triangle
*/
float Mesh::getArea(Triangle tri) {
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

// based, but far from copied from:
// https://github.com/Sunwinds/3D-Geometric-Features/blob/master/ConformalFactor/CalConformalFactor.h
/**
    Returns laplacian matrix of the mesh

    @return laplacian matrix
*/
Eigen::SparseMatrix<float> Mesh::getLapMatrix() {
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

/**
    Solves the conformal factor equation and fills the 
    m_confFactor vector using the SparseLU solver

    @param gaussDiff vector with the difference between
    target curvature and gaussian curvature
    @param laplacian laplacian matrix of the mesh
    @return vector with conformal factors of each vertex
*/
// https://scicomp.stackexchange.com/questions/21343/solving-linear-equations-using-eigen
std::vector<float> Mesh::solveConfFactor(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian) {
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


/**
    Solves the conformal factor equation and fills the 
    m_confFactor vector using the ConjugateGradient solver

    @param gaussDiff vector with the difference between
    target curvature and gaussian curvature
    @param laplacian laplacian matrix of the mesh
    @return vector with conformal factors of each vertex
*/
// TODO? multi threading
// recommended solver for 3D poisson eq. : https://eigen.tuxfamily.org/dox/group__TopicSparseSystems.html
std::vector<float> Mesh::solveConfFactorCG(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian) {
	Eigen::ConjugateGradient<Eigen::SparseMatrix<float>, Eigen::Lower|Eigen::Upper> cgLap;
	Eigen::VectorXf gaussDiffVect = Eigen::Map<Eigen::VectorXf, Eigen::Unaligned>(gaussDiff.data(), gaussDiff.size());
	Eigen::VectorXf confFact;

	cgLap.setMaxIterations(10000);
	cgLap.compute(laplacian);

	std::cout << "#iterations:     " << cgLap.iterations() << std::endl;
	std::cout << "estimated error: " << cgLap.error()      << std::endl;
	
	confFact = cgLap.solve(gaussDiffVect);
	std::vector<float> res(&confFact[0], confFact.data()+confFact.cols()*confFact.rows());
	return res;
}

//-----------------------------------------------------------------------------
//-------------------------calculateSignature----------------------------------
//-----------------------------------------------------------------------------

/**
    Outputs the signature as confFact.csv on the output folder
*/
void Mesh::calculateSignature () {
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

	printSignature(signature, m_positions.size() * 6); // *6 to account for added vertexes
}

/**
    Outputs the signature as confFact.csv on the output folder
    Prunes m_confFact after point insertion

    @param pruneFact total percentage of valued pruned on the
    extremities of the m_confFact vector, ex. 0.01 -> 1% 
    -> 0.5% of highest and 0.5% of lowest values pruned
*/
void Mesh::calculatePruneSignature (float pruneFact) {
	initializeSignature(-99, 100);

	std::deque<float> allConfFact;

	for (unsigned int i = 0; i < 5 * m_positions.size (); i++) 
	{
		int triangleIdx = getRandTri();
		Vec3f randPoint = getRandPoint(m_triangles[triangleIdx]);
		float randConfFactor = interpConfFactor(randPoint, triangleIdx);

		allConfFact.push_back(randConfFactor);
	}

	#ifdef TIMER
	std::cout << "Point generation ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	for (unsigned int i = 0; i < m_positions.size (); i++) 
	{
		allConfFact.push_back(m_confFacts[i]);
	}
	pruneDeq(allConfFact, pruneFact);

	auto result = std::minmax_element(allConfFact.begin(), allConfFact.end());

	minConf = allConfFact[result.first - allConfFact.begin()];
	maxConf = allConfFact[result.second - allConfFact.begin()];

	for (unsigned int i = 0; i < allConfFact.size (); i++) 
	{
		incrSignature(allConfFact[i], minConf, maxConf, -99, 100);
	}

	#ifdef TIMER
	std::cout << "Bin filling ended: " << (std::clock() - tInit) / (double)CLOCKS_PER_SEC << " seconds" << std::endl;
	#endif

	printSignature(signature, allConfFact.size()); // *6 to account for added vertexes
}

bool Mesh::greater (float i,float j) { return (i>j); }
bool Mesh::smaller (float i,float j) { return (i<j); }

/**
    Prunes conformal factor vector

    @param deq deque to be pruned
    @param pruneFactor total percentage of valued pruned on the
    extremities of the vector, ex. 0.01 -> 1% 
    -> 0.5% of highest and 0.5% of lowest values pruned
*/
void Mesh::pruneDeq(std::deque<float>& deq, float pruneFactor) {
	unsigned int pruneNum = deq.size () * pruneFactor;
	std::nth_element (deq.begin(), deq.begin() + pruneNum, deq.end(), greater);
	deq.erase(deq.begin(), deq.begin() + pruneNum);

	std::nth_element (deq.begin(), deq.begin() + pruneNum, deq.end(), smaller);
	deq.erase(deq.begin(), deq.begin() + pruneNum);
}

/**
    Initializes signature with correct number of bins and indexes

    @param min lowest bin index
    @param max highest bin index
*/
void Mesh::initializeSignature(int min, int max) {
	signature.resize (max - min + 1);

	for(unsigned int i = 0; i < signature.size(); i ++)
	{
		Bin b;
		b.occur = 0.0f;
		b.idx = min + i;
		signature[i] = b;
	}
}

/**
    Prints signature to .csv

    @param sig bin vector containing signature
    @param totalItens total number of occurences
*/
void Mesh::printSignature(std::vector<Bin> sig, unsigned int totalItems) {
	ofstream outfile;
    outfile.open ("output/confFact.csv");

    outfile << "binNumber, occurence\n";

	for(unsigned int i = 0; i < sig.size(); i ++)
	{
		float occurPercent = sig[i].occur / (float) totalItems;
		outfile << sig[i].idx << ", " << occurPercent << "\n";
	}

	outfile.close();
}

/**
    Returns a random triangle in the mesh, with
    prbability of being chosen being proportional
    to area

    @return triangle index
*/
// https://stackoverflow.com/questions/1761626/weighted-random-numbers
// http://forums.codeguru.com/showthread.php?320298-Searching-a-list-for-quot-closest-quot-floating-point-value-using-STL
// TODO? https://en.wikipedia.org/wiki/Alias_method
int Mesh::getRandTri() {
	float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float randArea = r * totalArea;


	std::vector<float>::iterator upIt = std::upper_bound(m_areas.begin(), m_areas.end(), randArea, FloatLessThan());;
	return m_sortedAreasIdx[upIt - m_areas.begin()];
}

/**
    Returns random vertex inside of triangle

    @param tri triangle
    @return random point inside given triangle
*/
// http://www.cs.princeton.edu/~funk/tog02.pdf
// p8, equation (1)
Vec3f Mesh::getRandPoint(Triangle tri) {
	float r1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
	float r2 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

	return Vec3f((1 - sqrt(r1)) * m_positions[tri[0]] + 
				sqrt(r1) * (1 - r2) * m_positions[tri[1]] + 
				sqrt(r1) * r2 * m_positions[tri[2]]);
}

/**
    Gives the interpolated value of the conformal factor
    using the distance to other points as weight

    @param sig bin vector containing signature
    @param totalItens total number of occurences
    @return interpolated conformal factor
*/
// https://stackoverflow.com/questions/18755251/linear-interpolation-of-three-3d-points-in-3d-space
// section 3.4, use barycentric coordinates to find weight associated to each vertex
float Mesh::interpConfFactor(Vec3f point, unsigned int triIdx) {
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

/**
    Increments signature according to conformal
    factor of vertex

    @param confFact conformal factor of vertex
    @param min lowest conformal factor
    @param max highest conformal factor
    @param binMin lowest bin index
    @param binMax highest bin index
*/
void Mesh::incrSignature(float confFact, float min, float max, int binMin, int binMax) {
	float normIdx = (confFact - min) / (max - min);
	int equivBin = normIdx * (binMax - binMin);
	signature[equivBin].occur += 1;
}

//-----------------------------------------------------------------------------
//-------------------------------helpers---------------------------------------
//-----------------------------------------------------------------------------

/**
    Returns verices in a 1-neighborhood to given vertex

    @param tri all indexes of triangles a vertex is contianed in
    @param point vertex index
    @return set with neighbouring points
*/
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

/**
    Return all indexes of triangles the vertex is contianed in
    Useful to find all triangles two points can have in common

    @param point vertex index
    @param triangles 
    @return set with neighbouring points
*/
std::vector<unsigned int> Mesh::containPoint(unsigned int point, std::vector<unsigned int> triangles) {

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

/**
    Return all indexes of triangles the vertex is contianed in
    Useful to find all triangles two points can have in common

    @param point1 first point in the triangle
    @param second point in the tirangle
    @param triangle triangl eindex 
    @return cotangent of the vertex of the triangle that's not
    point1 or point2
*/
float Mesh::cotan(unsigned int point1, unsigned int point2, unsigned int triangle) {
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

/**
    scale to the unit cube and center at original
*/
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
