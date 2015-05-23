       
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


void Vector4::dump(const core::Lisp_sp& lisp )
{_G();
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





Vector4 Vector4::normalized(const core::Lisp_sp& lisp)
{_G();
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


