       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *      File:   spline.h
 *
 *
 *
 *	A general spline class.
 *
 *	Given a set of control points it allows the evaluation of the spline.
 *
 *	The CATMULLROM spline is guaranteed to pass through the control
 *	points.
 */


#ifndef	SPLINE_H
#define SPLINE_H

#ifdef	USEBOOSTPYTHON
#include <Python.h>
#endif

#include <vector>
#include <map>
#include <string>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>



//#define	SPLINE_BSPLINE		1
#define	SPLINE_CATMULLROM	2
//#define	SPLINE_FULLBSPLINE	3

#include <cando/chem/chemPackage.h>


namespace chem {



class	Dumb_Spline;
    typedef	gctools::smart_ptr<Dumb_Spline>	RPSpline;
class	Dumb_Spline{

public:
    friend	RPSpline	Spline( int iSplineType, const VectorVector3s& points );

    Vector3			evaluateAtIndex( double index );
    double			distanceFromPoint( double index, Vector3& pnt );



#if 0
    Vector3			evaluateAtDistance( double x );
#endif
	void	describe(core::Lisp_sp lisp);
	Dumb_Spline(int iSplineType, const VectorVector3s& points );
double	domainMax()	{return this->gmMatrices.size();};
double	domainMin()	{return 0.0;};


private:
	int		iType;
	vector<Matrix>	gmMatrices;
	Matrix		mbAMatrix;
	Matrix		mbBMatrix;
	Matrix		mbsMatrix;
	Matrix		mbYMatrix;
	Matrix		mbZMatrix;

	Matrix	gMatrixForPoints( const VectorVector3s& points, int i4 );
void	defineMatrices( int iSplineType );
void	definePoints( const VectorVector3s& points );

};


inline	RPSpline	Spline(int iSplineType,const VectorVector3s& points) {
    	return RPSpline(new Dumb_Spline(iSplineType,points));
}

};
#endif  /* ifdef SPLINE_H*/
