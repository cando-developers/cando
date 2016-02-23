/*
    File: vector4.cc
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




#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <iostream>
#include <cando/geom/vector4.h>
//#include "core/omath.h"


//
// Constructors
//
Vector4::Vector4()
{

    this->coords[0] = 0.0;
    this->coords[1] = 0.0;
    this->coords[2] = 0.0;
    this->coords[3] = 0.0;
}

//
// Destructor
//
Vector4::~Vector4()
{
}


Vector4::Vector4( double w, double x, double y, double z )
{

    this->coords[0] = w;
    this->coords[1] = x;
    this->coords[2] = y;
    this->coords[3] = z;
}


void Vector4::dump()
{
    _lisp->print(BF("<%lf,%lf,%lf,%lf>")
		% this->coords[0]
		% this->coords[1]
		% this->coords[2]
		% this->coords[3] );
}


string	Vector4::asString()
{
stringstream ss;
    ss << "(" << this->coords[0] << ",";
    ss << this->coords[1] << ",";
    ss << this->coords[2] << ",";
    ss << this->coords[3] << ")";
    return ss.str();
}

void Vector4::writeToStream( std::ostream& out )
{
	out << "(";
	out << this->coords[0] << "@";
	out << this->coords[1] << "@";
	out << this->coords[2] << "@";
	out << this->coords[3] << ")";
}





Vector4 Vector4::normalized()
{
double	l;
Vector4	v;

    l = this->length();
    if ( l > 0.0 ) {
	v.coords[0] = this->coords[0]/l;
	v.coords[1] = this->coords[1]/l;
	v.coords[2] = this->coords[2]/l;
	v.coords[3] = this->coords[3]/l;
	return v;
    }
    SIMPLE_ERROR(BF("Attempted to normalize a zero vector"));
}


double Vector4::dotProduct( Vector4 v )
{
	return(  this->getW()*v.getW()
		+this->getX()*v.getX()
		+this->getY()*v.getY()
		+this->getZ()*v.getZ());
}


double Vector4::length()
{
	return(sqrt(
		 this->getW()*this->getW()
		+this->getX()*this->getX()
		+this->getY()*this->getY()
		+this->getZ()*this->getZ()));
}

int	Vector4::indexOfLargestElement()
{
int	iMax,i;
double	dMax;
    iMax = 0;
    dMax = -9999999.9;
    for (i=0;i<4;i++) {
	if ( dMax < this->coords[i] ) {
	    dMax = this->coords[i];
	    iMax = i;
	}
    }
    return iMax;
}


#if 0
#ifndef WIN32
ostream& operator<<( std::ostream& out, Vector4 v )
{
	v.dump(out);
	return out;
}
#endif
#endif
