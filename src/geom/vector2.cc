/*
    File: vector2.cc
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
       
#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//


#include <iostream>
#include <sstream>
#include <clasp/core/foundation.h>
#include <cando/geom/vector2.h>
#include <math.h>




namespace geom
{



double	calculate2Distance( const Vector2& va,
			   const Vector2& vb)
{
    Vector2 vc = va - vb;
    return vc.length();
}

double	calculate2DistanceSquared( const Vector2& va,
				  const Vector2& vb)
{
    Vector2 vc = va - vb;
    return vc.dotProduct(vc);
}


/*! Return the angle in radians
 */
double calculate2Angle( const Vector2& va,
			const Vector2& vb,
		       const Vector2& vc)
{
    Vector2	vab = (va-vb).normalized();
    Vector2	vcb = (vc-vb).normalized();
    double ang = acos(vab.dotProduct(vcb));
    return ang;
}



//
// Constructors
//
Vector2::Vector2()
{

    this->coords[0] = 0.0;
    this->coords[1] = 0.0;
}

//
// Destructor
//
Vector2::~Vector2()
{
}


Vector2::Vector2( double x, double y )
{

    this->coords[0] = x;
    this->coords[1] = y;
}

double	Vector2::get(int i) const
{
    HARD_ASSERT(i==0||i==1);
    return this->coords[i];
};



#ifdef XML_ARCHIVE
    void	Vector2::archive(core::ArchiveP node)
{
    node->attribute("x",this->coords[0]);
    node->attribute("y",this->coords[1]);
    LOG("Serialized vector@{} = {}, {}" , this , this->coords[0] , this->coords[1]  );
};
#endif


void Vector2::dump()
{
    std::cout << this->asString() << std::endl;
}


bool Vector2::isDefined()
{
    if ( this->coords[0] != 0.0 ) return true;
    if ( this->coords[1] != 0.0 ) return true;
    return false;
}



    Vector2 Vector2::normalized() const
{
double	l;
Vector2	v;

    l = this->length();
    if ( l > 0.0 ) {
	v.coords[0] = this->coords[0]/l;
	v.coords[1] = this->coords[1]/l;
	return v;
    }
    THROW_HARD_ERROR("Tried to normalize a zero vector") ;
}


double Vector2::dotProduct( Vector2 v )
{
    return(this->getX()*v.getX() +this->getY()*v.getY() );
}


double Vector2::length() const
{
    return(sqrt(this->getX()*this->getX()+this->getY()*this->getY() ));
}



string	Vector2::asString() const
{
stringstream	ss;
    ss.str("");
    ss << "(" << this->getX()
	<< "," <<this->getY()
	<< ")";
    return ss.str();
}


};
