       
#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//


#include <iostream>
#include <sstream>
#include <clasp/core/foundation.h>
#include <cando/candoBase/vector2.h>
#include <math.h>




namespace candoBase
{



double	calculateDistance( const Vector2& va,
			   const Vector2& vb,
			   const core::Lisp_sp& lisp)
{
    Vector2 vc = va - vb;
    return vc.length();
}

double	calculateDistanceSquared( const Vector2& va,
				  const Vector2& vb)
{
    Vector2 vc = va - vb;
    return vc.dotProduct(vc);
}


/*! Return the angle in radians
 */
double calculateAngle( const Vector2& va,
			const Vector2& vb,
		       const Vector2& vc,
		       const core::Lisp_sp& lisp)
{
    Vector2	vab = (va-vb).normalized(lisp);
    Vector2	vcb = (vc-vb).normalized(lisp);
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
{_G();
    node->attribute("x",this->coords[0]);
    node->attribute("y",this->coords[1]);
    LOG(BF("Serialized vector@%X = %lf, %lf") % this % this->coords[0] % this->coords[1]  );
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



    Vector2 Vector2::normalized(core::Lisp_sp lisp) const
{_G();
double	l;
Vector2	v;

    l = this->length();
    if ( l > 0.0 ) {
	v.coords[0] = this->coords[0]/l;
	v.coords[1] = this->coords[1]/l;
	return v;
    }
    THROW_HARD_ERROR(BF("Tried to normalize a zero vector")) ;
}


double Vector2::dotProduct( Vector2 v )
{
    return(this->getX()*v.getX() +this->getY()*v.getY() );
}


double Vector2::length() const
{
    return(sqrt(this->getX()*this->getX()+this->getY()*this->getY() ));
}



string	Vector2::asString()
{
stringstream	ss;
    ss.str("");
    ss << "(" << this->getX()
	<< "," <<this->getY()
	<< ")";
    return ss.str();
}


};
