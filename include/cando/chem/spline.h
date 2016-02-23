/*
    File: spline.h
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
       
       
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
	void	describe();
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
