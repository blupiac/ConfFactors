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
#include <ctime>
#include <deque>

#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/IterativeLinearSolvers>

#include "Vec3.h"
#include "Triangle.h"

//#define DEBUG
#define TIMER

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

struct FloatLessThan
{
    bool operator() (const float & left, const float & right)
    {
        return left < right;
    }
};


/// A Mesh class, storing a list of vertices and a list of triangles indexed over it.
class Mesh {
public:
    inline Mesh () {}
    inline virtual ~Mesh () {}

    inline std::vector<Vec3f> & positions () { return m_positions; }
    inline const std::vector<Vec3f> & positions () const { return m_positions; }
    inline std::vector<Triangle> triangles () { return m_triangles; }
    inline const std::vector<Triangle> & triangles () const { return m_triangles; }
    inline  std::vector<float> & confFacts () { return m_confFacts; }

    /// Empty the positions, triangles arrays.
    void clear ();

	/**
        Parses and does the pre-processing of an off file

        @param filename filepath of an off file to be loaded
    */
	void loadOFF (const std::string & filename);

    /**
        Fills the m_confFact vector with the confFactors of the mesh
    */
	void calculateConfFact ();
    /**
        Outputs the signature as confFact.csv on the output folder
    */
	void calculateSignature ();
    /**
        Outputs the signature as confFact.csv on the output folder
        Prunes m_confFact after point insertion

        @param pruneFact total percentage of valued pruned on the
        extremities of the m_confFact vector, ex. 0.01 -> 1% 
        -> 0.5% of highest and 0.5% of lowest values pruned
    */
    void calculatePruneSignature (float pruneFact);
    /**
        Returns the given conformal factor after being normalized 

        @param confIdx index of conformal factor to be returned
        @return normalized conformal factor
    */
    float normalizeConf(unsigned int confIdx);
    #ifdef DEBUG
    /**
        Returns the given gaussian curvature after being normalized 

        @param gaussIdx index of conformal factor to be returned
        @return normalized gausian curvature
    */
    float normalizeGausscurv(unsigned int gaussIdx);
    #endif

    /**
        scale to the unit cube and center at original
    */
    void centerAndScaleToUnit ();


private:
    std::vector<Vec3f> m_positions;
    std::vector<float> m_confFacts;
    std::vector<float> m_gaussDiff;
    std::vector<float> m_areas;
    std::vector<unsigned int> m_sortedAreasIdx;
    #ifdef DEBUG
    std::vector<float> m_gausscurv;
    #endif
    std::vector<Bin> signature;
    std::vector<Triangle> m_triangles;
    std::vector<std::vector<unsigned int> > m_nneighbours;
    float totalArea, totalCurv, minConf, maxConf;
    #ifdef DEBUG
    float minGauss, maxGauss;
    #endif
    #ifdef TIMER
    std::clock_t tInit;
    #endif
    
    /**
        Returns gaussian curvature of a vertex

        @param i vertex index
        @return gaussian curvature of vertex
    */
    float getGaussCurv(unsigned int i);
    /**
        Returns angle of a vertex in a triangle

        @param tri triangle inside of which the vertex is in
        @param pointIdx vertex index
        @return angle of vertex
    */
    float getAngle(Triangle tri, int pointIdx);
    /**
        Returns target curvature of a vertex

        @param i vertex index
        @return target curvature of vertex
    */
	float getTargetCurv(unsigned int i);
    /**
        Returns area of a triangle

        @param i triangle index
        @return area of triangle
    */
    float getArea(unsigned int i);
    /**
        Returns area of a triangle

        @param tri triangle
        @return area of triangle
    */
    float getArea(Triangle tri);
    /**
        Returns whether the vertex is inside a mesh or in its border

        @param ptIdx vertex index
        @return boolean saying whether vertex is in border or not
    */
    bool isBorder(unsigned int ptIdx);
    /**
        Returns laplacian matrix of the mesh

        @return laplacian matrix
    */
    Eigen::SparseMatrix<float> getLapMatrix();
    /**
        Solves the conformal factor equation and fills the 
        m_confFactor vector using the SparseLU solver

        @param gaussDiff vector with the difference between
        target curvature and gaussian curvature
        @param laplacian laplacian matrix of the mesh
        @return vector with conformal factors of each vertex
    */
    std::vector<float> solveConfFactor(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian);
    /**
        Solves the conformal factor equation and fills the 
        m_confFactor vector using the ConjugateGradient solver

        @param gaussDiff vector with the difference between
        target curvature and gaussian curvature
        @param laplacian laplacian matrix of the mesh
        @return vector with conformal factors of each vertex
    */
    std::vector<float> solveConfFactorCG(std::vector<float> gaussDiff, Eigen::SparseMatrix<float> laplacian);

    /**
        Initializes signature with correct number of bins and indexes

        @param min lowest bin index
        @param max highest bin index
    */
	void initializeSignature(int min, int max);
    /**
        Prints signature to .csv

        @param sig bin vector containing signature
        @param totalItens total number of occurences
    */
    void printSignature(std::vector<Bin> sig, unsigned int totalItems);
    /**
        Returns a random triangle in the mesh, with
        prbability of being chosen being proportional
        to area

        @return triangle index
    */
	int getRandTri();
    /**
        Returns random vertex inside of triangle

        @param tri triangle
        @return random point inside given triangle
    */
	Vec3f getRandPoint(Triangle tri);
    /**
        Gives the interpolated value of the conformal factor
        using the distance to other points as weight

        @param sig bin vector containing signature
        @param totalItens total number of occurences
        @return interpolated conformal factor
    */
	float interpConfFactor(Vec3f point, unsigned int triIdx);
    /**
        Increments signature according to conformal
        factor of vertex

        @param confFact conformal factor of vertex
        @param min lowest conformal factor
        @param max highest conformal factor
        @param binMin lowest bin index
        @param binMax highest bin index
    */
	void incrSignature(float confFact, float min, float max, int binMin, int binMax);

    static bool greater (float i,float j);
    static bool smaller (float i,float j);
    /**
        Prunes conformal factor vector

        @param deq deque to be pruned
        @param pruneFactor total percentage of valued pruned on the
        extremities of the vector, ex. 0.01 -> 1% 
        -> 0.5% of highest and 0.5% of lowest values pruned
    */
    void pruneDeq(std::deque<float>& deq, float pruneFactor);
    /**
        Introduces noise into the mesh

        @param noiseCoef noise coefficient, values between
        0.01-0.1 give acceptable results, read readme for more
        information
    */
    void makeNoise(float noiseCoef);

    // helpers
    /**
        Returns verices in a 1-neighborhood to given vertex

        @param tri all indexes of triangles a vertex is contianed in
        @param point vertex index
        @return set with neighbouring points
    */
	std::set<unsigned int> getVoisins(std::vector<unsigned int> tri, unsigned int point);
    /**
        Return all indexes of triangles the vertex is contianed in
        Useful to find all triangles two points can have in common

        @param point vertex index
        @param triangles 
        @return set with neighbouring points
    */
    std::vector<unsigned int> containPoint(unsigned int point, std::vector<unsigned int> triangles);
    /**
        Return all indexes of triangles the vertex is contianed in
        Useful to find all triangles two points can have in common

        @param point1 first point in the triangle
        @param second point in the tirangle
        @param triangle triangl eindex 
        @return cotangent of the vertex of the triangle that's not
        point1 or point2
    */
    float cotan(unsigned int point1, unsigned int point2, unsigned int triangle);
};
