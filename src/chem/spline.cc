/*
    File: spline.cc
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


#if 1
#define	TURN_DEBUG_OFF
#define	MATHEMATICA_DEBUG_ON	1
#define	PARMDEBUG	1
#endif


#include <math.h>
#include <stdio.h>
#include <iomanip>

#include <clasp/core/common.h>
#include <clasp/core/lispStream.h>
#include <cando/geom/matrix.h>
#include <cando/chem/spline.h>



#if 0
#define		DEBUGINNER	1
#define		DEBUGTORSION	1
#define		DEBUGANGLE	1
#define		DEBUGBOND	1
#define		DEBUGTORSION	1
#endif
#define		DEBUGTORSION	1


using namespace chem;

void	Dumb_Spline::defineMatrices(int iSplineType )
{
    if ( iSplineType == SPLINE_CATMULLROM ) {
	double mSpline[] = { 	-1,  2, -1,  0,
				 3, -5,  0,  2,
				-3,  4,  1,  0,
				 1, -1,  0,  0 };
	Matrix m;
	m.setFromDoubleArray(16,mSpline);
	m = m*0.5;
	this->mbAMatrix = m;
	this->mbBMatrix = m;
	this->mbsMatrix = m;
	this->mbYMatrix = m;
	this->mbZMatrix = m;
    }
}

/*
 *	Constructor
 *
 *
 */
Dumb_Spline::Dumb_Spline(int iSplineType, const VectorVector3s& points )
{
    this->iType = iSplineType;
    this->defineMatrices(iSplineType);
    this->definePoints(points);
}



Matrix	Dumb_Spline::gMatrixForPoints( const VectorVector3s& points, int i4 )
{
Matrix	m;
Vector3	p1, p2, p3, p4;

   p1 = points[i4-3];
   p2 = points[i4-2];
   p3 = points[i4-1];
   p4 = points[i4];
   m.at(0, 0) = p1.getX();
   m.at(1, 0) = p1.getY();
   m.at(2, 0) = p1.getZ();
   m.at(3, 0) = 1.0;

   m.at(0, 1) = p2.getX();
   m.at(1, 1) = p2.getY();
   m.at(2, 1) = p2.getZ();
   m.at(3, 1) = 1.0;

   m.at(0, 2) = p3.getX();
   m.at(1, 2) = p3.getY();
   m.at(2, 2) = p3.getZ();
   m.at(3, 2) = 1.0;

   m.at(0, 3) = p4.getX();
   m.at(1, 3) = p4.getY();
   m.at(2, 3) = p4.getZ();
   m.at(3, 3) = 1.0;
   return m;
}




void	Dumb_Spline::definePoints( const VectorVector3s& points )
{
int	gmSize;
Matrix	m, g, gm;
int	i;
    gmSize = points.size() - 3;
    this->gmMatrices.resize(gmSize,m);
    for ( i=0; i<gmSize; i++ ) {
	m = this->mbsMatrix;
	if ( i == 1 ) {
	    m = this->mbBMatrix;
	} else if (i == 0 ) {
	    m = this->mbAMatrix;
	} else if (i == gmSize-2) {
	    m = this->mbYMatrix;
	} else if (i == gmSize-1) {
	    m = this->mbZMatrix;
	}
	g = this->gMatrixForPoints( points, i+3 );
	this->gmMatrices[i] = g * m;
    }
}

void	Dumb_Spline::describe()
{
    Matrix	m, gm;
    uint	i;
    core::clasp_write_string(fmt::format( "spline type = {}\n" , this->iType ));
    core::clasp_write_string("gmMatrices ======>\n");
    for ( i=0; i<this->gmMatrices.size(); i++ ) {
	this->gmMatrices[i].dump();
        core::clasp_write_string("--------------\n");
    }
    core::clasp_write_string("mbAMatrix ======>\n");
    this->mbAMatrix.dump();
    core::clasp_write_string("mbBMatrix ======>\n");
    this->mbBMatrix.dump();
    core::clasp_write_string("mbsMatrix ======>\n");
    this->mbsMatrix.dump();
    core::clasp_write_string("mbYMatrix ======>\n");
    this->mbYMatrix.dump();
    core::clasp_write_string("mbZMatrix ======>\n");
    this->mbZMatrix.dump();
}


//
//	evaluateAtIndex
//
//	Evaluate the spline at the point ti.
//
//	0<=ti<1 evaluates over the first spline

Vector3	Dumb_Spline::evaluateAtIndex( double ti )
{
uint	i;
double	t;
double	t2,t3;
Vector3	tv;
Matrix	m;
    i = floor(ti);
    if ( i>=this->gmMatrices.size() ) {
	i = this->gmMatrices.size()-1;
	t = 1.0;
    } else {
	t = ti - i;
    }
    m = this->gmMatrices[i];
    t2 = t*t;
    t3 = t2*t;
    tv.set(t3,t2,t);
    return(m*tv);
}



//
//	distanceFromPoint
//
//	Return the vector distance between a point on the spline and
//	a given point.
double	Dumb_Spline::distanceFromPoint(	double	index, Vector3& pnt )
{
    Vector3	sp, diff;

    sp = this->evaluateAtIndex(index);
    diff = pnt - sp;
    return diff.length();
}
