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

#pragma once
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/Eigenvalues>
#include <eigen3/Eigen/Sparse>
#include "Vec3.h"
#include "Triangle.h"

#define DEBUG

struct Edge {
 unsigned int i , j;
 Edge(unsigned int a , unsigned int b) : i(std::min(a,b)) , j(std::max(a,b)) {}
 bool operator == (Edge const & other) const { return ( ( (i == other.i) && (j == other.j) ) ||
 														( (i == other.j) && (j == other.i) ) ); }
 bool operator < (Edge const & other) const { 	if(i != other.i)	return (i < other.i);
 												else				return (j < other.j); }
};

struct Bin 
{
	int idx;
    float occur;
};

/// A Mesh class, storing a list of vertices and a list of triangles indexed over it.
class Mesh {
public:
    inline Mesh () {}
    inline virtual ~Mesh () {}

    inline std::vector<Vec3f> & positions () { return m_positions; }
    inline const std::vector<Vec3f> & positions () const { return m_positions; }
    inline  std::vector<Vec3f> & normals () { return m_normals; }
    inline const std::vector<Vec3f> & normals () const { return m_normals; }
    inline std::vector<Triangle> triangles () { return m_triangles; }
    inline const std::vector<Triangle> & triangles () const { return m_triangles; }
    inline  std::vector<float> & confFacts () { return m_confFacts; }

    /// Empty the positions, normals and triangles arrays.
    void clear ();

	/// Loads the mesh from a <file>.off
	void loadOFF (const std::string & filename);

	void calculateConfFact ();
	void calculateSignature ();
    float normalizeConf(unsigned int confIdx);
    #ifdef DEBUG
    float normalizeLaplacian(unsigned int lapIdx);
    float normalizeGausscurv(unsigned int gaussIdx);
    #endif

    /// Compute smooth per-vertex normals
    void recomputeNormals ();
    
    void laplacianFilter (float alpha);
    void laplacianFilterGeom(float alpha);

	void simplify(unsigned int resolution);
	void simplifyAdaptive(unsigned int n);

	void subdivide();

    /// scale to the unit cube and center at original
    void centerAndScaleToUnit ();

	void addFloor (unsigned int type);
    void removeFloor (unsigned int type);

private:
    std::vector<Vec3f> m_positions;
    std::vector<Vec3f> m_normals;
    std::vector<float> m_confFacts;
    std::vector<float> m_areas;
    std::vector<unsigned int> m_sortedAreasIdx;
    #ifdef DEBUG
    std::vector<float> m_gausscurv;
    std::vector<float> m_laplacian;
    #endif
    std::vector<Bin> signature;
    std::vector<Triangle> m_triangles;
    // TODO: send indexes instead of triangles, this way getTargetCurv can be optimized
    std::vector<std::vector<Triangle> > m_nneighbours;
    std::map <Edge, unsigned int> middle_points;
    std::map <std::pair<Triangle, unsigned int>, float> cotans;
    float totalArea, totalCurv, totalConf, minConf, maxConf;
    #ifdef DEBUG
    float minLap, maxLap, minGauss, maxGauss;
    #endif
    
    float getGaussCurv(unsigned int i);
    float getAngle(Triangle tri, int pointIdx);
	float getTargetCurv(unsigned int i);
    float getArea(unsigned int i);
    float getArea(Triangle tri);
	float getLaplacian(unsigned int i);
    float getAMixed(unsigned int i);
    float voronoiRegion(unsigned int ptIdx);
    float voronoiArea(unsigned int ptIdx, Triangle t);
    bool isObtuse(Triangle t);
    bool isBorder(unsigned int ptIdx);
    void purgeConf(float coef);
    Eigen::SparseMatrix<float> getLapMatrix();
    std::vector<float> solveConfFactor(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian);

	void initializeSignature(int min, int max);
    void printSignature(std::vector<Bin> sig, unsigned int totalItems);
	int getRandTri();
	Vec3f getRandPoint(Triangle tri);
	float interpConfFactor(Vec3f point, unsigned int triIdx);
	void incrSignature(float confFact, float min, float max, int binMin, int binMax);

    void getMaxMin(float& max_x, float& max_y, float& max_z, float& min_x,
				float& min_y, float& min_z);
				
	std::vector<unsigned int> allInside(float max_x, float max_y, float max_z,
										float min_x, float min_y, float min_z);
										
	void simplifySubOct(float max_x, float max_y, float max_z,
							float min_x, float min_y, float min_z,
							unsigned int n, std::vector<Vec3f>& newPositions, 
							std::vector<Vec3f>& newNormals, std::vector<unsigned int>& oldToNew);
	
	std::set<unsigned int> getVoisins(std::vector<Triangle> tri, unsigned int point);
    Vec3f calculerBarycentre(std::set<unsigned int> points);
    Vec3f calculerBarycentreTri(Triangle tri);
    Vec3f calculerBarycentreGeom(unsigned int point, std::set<unsigned int> points, std::vector<Triangle> triangles);
    
    std::vector<Triangle> containPoint(unsigned int point, std::vector<Triangle> triangles);
    float cotan(unsigned int point1, unsigned int point2, Triangle triangle);
    
    void fillMiddle();
    void triAveraging(unsigned int limit);
    void subdivideTri();
    void recomputeNN();

    std::vector<Triangle> intersectTriVect(std::vector<Triangle> v1, std::vector<Triangle> v2);
    
};
