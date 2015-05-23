#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <iostream>
#include <sstream>
#include <cando/geom/vector3.h>
#include <clasp/core/cons.h>
#include <clasp/core/serialize.h>
#include <math.h>
//#include "core/archiveNode.h"
#include <clasp/core/lispStream.h>
#include <cando/geom/matrix.h>
#include <clasp/core/wrappers.h>



//
// Constructors
//
Vector3::Vector3()
{
    this->coords[0] = 0.0;
    this->coords[1] = 0.0;
    this->coords[2] = 0.0;
}

//
// Destructor
//
Vector3::~Vector3()
{
    //nothing
}


Vector3::Vector3( double x, double y, double z )
{

    this->coords[0] = x;
    this->coords[1] = y;
    this->coords[2] = z;
}

Vector3 Vector3::inNanometers() const
{
    return Vector3(this->coords[0]*0.1,
		   this->coords[1]*0.1,
		   this->coords[2]*0.1);
}


void Vector3::dump(const core::Lisp_sp& lisp)
{
    lisp->print(BF("<%lf,%lf,%lf>") % this->coords[0] % this->coords[1] % this->coords[2] );
}

void Vector3::writeToStream( std::ostream& out )
{
	out << "(";
	out << this->coords[0] << ",";
	out << this->coords[1] << ",";
	out << this->coords[2] << ")";
}




void Vector3::write(core::T_sp stream)
{
    stringstream ss;
    this->writeToStream(ss);
    core::clasp_write_string(ss.str(),stream);
}




Vector3 Vector3::crossProduct( const Vector3& v ) const
{
double	x,y,z;

	x = this->getY()*v.getZ()-this->getZ()*v.getY();
	y = this->getZ()*v.getX()-this->getX()*v.getZ();
	z = this->getX()*v.getY()-this->getY()*v.getX();
	return Vector3(x,y,z);
}

Vector3 Vector3::normalized() const
{
double	l;
Vector3	v;

    l = this->length();
    if ( l > 0.0 ) {
	v.coords[0] = this->coords[0]/l;
	v.coords[1] = this->coords[1]/l;
	v.coords[2] = this->coords[2]/l;
	return v;
    }
    THROW_HARD_ERROR(BF("Attempted to normalize the vector (0,0,0)"));
}




Vector3 Vector3::normalized(const core::Lisp_sp& lisp) const
{_G();
double	l;
Vector3	v;

    l = this->length();
    if ( l > 0.0 ) {
	v.coords[0] = this->coords[0]/l;
	v.coords[1] = this->coords[1]/l;
	v.coords[2] = this->coords[2]/l;
	return v;
    }
    SIMPLE_ERROR(BF("Attempted to normalize the vector (0,0,0)"));
}

Vector3 Vector3::normalizedOrZero() const
{
double	l;
Vector3	v;

    l = this->length();
    if ( l > 0.0 ) {
	v.coords[0] = this->coords[0]/l;
	v.coords[1] = this->coords[1]/l;
	v.coords[2] = this->coords[2]/l;
	return v;
    }
    v.coords[0] = 0.0;
    v.coords[1] = 0.0;
    v.coords[2] = 0.0;
    return v;
}



double Vector3::dotProduct( const Vector3& v )
{
	return(this->getX()*v.getX()
		+this->getY()*v.getY()
		+this->getZ()*v.getZ());
}


double Vector3::length() const
{
	return(sqrt(this->getX()*this->getX()
		+this->getY()*this->getY()
		+this->getZ()*this->getZ()));
}



string	Vector3::asString() const
{
stringstream	ss;
    ss.str("");
    ss << std::setw(5) << std::setprecision(4);
    ss << "(" <<this->getX()
	<< "," <<this->getY()
	<< "," <<this->getZ()
	<< ") ";
    return ss.str();
}

void Vector3::parseFromString(string& str)
{
char	paren, comma;
stringstream	ss;
double		x,y,z;
    ss.str(str);
    ss >> std::skipws >> paren >> std::skipws;
    ss >> x;
    ss >> std::skipws >> comma >> std::skipws;
    ss >> y;
    ss >> std::skipws >> comma >> std::skipws;
    ss >> z;
    this->set(x,y,z);
    LOG(BF("Vector3::parseFromString got (%lf,%lf,%lf)") % x % y % z );
}


void Vector3::parseFromCharacterArray(char* str)
{
char	space, comma;
stringstream	ss;
double		x,y,z;
    ss.str(str);
    ss >> space;
    ss >> x;
    ss >> comma;
    ss >> y;
    ss >> comma;
    ss >> z;
    this->set(x,y,z);
    LOG(BF("Vector3::parseFromCharacterArray got (%lf,%lf,%lf)") % x % y % z );
}





#if 0 //[
QDomNode_sp	Vector3::asXml(string nm)
{
QDomNode_sp	nd;
    nd = QDomNode_O::create(_lisp,nm);
    nd->addAttributeDoubleScientific("x",this->getX());
    nd->addAttributeDoubleScientific("y",this->getY());
    nd->addAttributeDoubleScientific("z",this->getZ());
    return nd;
}


void	Vector3::parseFromXml(QDomNode_sp node)
{
    double x = node->getAttributeDouble("x");
    double y = node->getAttributeDouble("y");
    double z = node->getAttributeDouble("z");
    this->set(x,y,z);
}
#endif //]



#define	DIFF_TOL	0.0001
bool	Vector3::sameAs(const Vector3& v)
{
Vector3	diff;
    diff = this->sub(v);
    if ( fabs(diff.getX()) < DIFF_TOL
		&& fabs(diff.getY()) < DIFF_TOL
		&& fabs(diff.getZ()) < DIFF_TOL ) {
	return true;
    }
    return false;
}


double	Vector3::angleToVectorAboutNormal(const Vector3& toVector, const Vector3& aboutNormal, const core::Lisp_sp& lisp )
{_G();
	Vector3 men = this->normalizedOrZero();
	LOG(BF("me.normalized{%lf,%lf,%lf} ") % men.getX() % men.getY() % men.getZ() );
	Vector3 ton = toVector.normalizedOrZero();
	LOG(BF("to.normalized{%lf,%lf,%lf} ") % ton.getX() % ton.getY() % ton.getZ() );
	LOG(BF("about_normal{%lf,%lf,%lf} ") % aboutNormal.getX() % aboutNormal.getY() % aboutNormal.getZ() );
	double cosTheta = men.dotProduct(ton);
	if ( fabs(cosTheta-1.0)<0.00001 )  return 0.0;
	if ( fabs(cosTheta+1.0)<0.00001 )  return 3.14159;
	LOG(BF("cosTheta = %lf") % cosTheta  );
	Vector3 cr = men.crossProduct(ton);
	LOG(BF("status") );
	double dir = cr.dotProduct(aboutNormal);
	LOG(BF("dir = %lf ") % dir );
	double theta = (dir>=0.0)?acos(cosTheta):(6.2831853-acos(cosTheta));
	LOG(BF("theta = %lf deg") % theta/0.0174533 );
	return theta;
    }



void	Vector3::archive(core::ArchiveP node)
{_G();
    node->attribute("x",this->coords[0]);
    node->attribute("y",this->coords[1]);
    node->attribute("z",this->coords[2]);
    LOG(BF("Serialized vector@%X = %lf, %lf, %lf") % this % this->coords[0] % this->coords[1] % this->coords[2]  );
};

#if 0
void Vector3::serialize(core::SerializerNode node)
{
    node
	->orderedPOD(this->coords[0])
	->orderedPOD(this->coords[1])
	->orderedPOD(this->coords[2])
	;
}
#endif
void Vector3::fillFromCons(core::Cons_sp vals)
{
  double x = core::clasp_to_double(core::oFirst(vals).as<core::Number_O>());
  double y = core::clasp_to_double(core::oSecond(vals).as<core::Number_O>());
  double z = core::clasp_to_double(core::oThird(vals).as<core::Number_O>());
  this->coords[0] = x;
  this->coords[1] = y;
  this->coords[2] = z;
}
    


double	calculateDistance( const Vector3& va,
			   const Vector3& vb,
			   const core::Lisp_sp& lisp)
{
    Vector3 vc = va - vb;
    return vc.length();
}

double	calculateDistanceSquared( const Vector3& va,
				  const Vector3& vb)
{
    Vector3 vc = va - vb;
    return vc.dotProduct(vc);
}


/*! Return the angle in radians
 */
double calculateAngle( const Vector3& va,
			const Vector3& vb,
		       const Vector3& vc,
		       const core::Lisp_sp& lisp)
{
    Vector3	vab = (va-vb).normalized(lisp);
    Vector3	vcb = (vc-vb).normalized(lisp);
    double ang = acos(vab.dotProduct(vcb));
    return ang;
}


/*! Return the dihedral in radians
 */
double calculateDihedral( const Vector3& va,
			const Vector3& vb,
			const Vector3& vc,
			const Vector3& vd,
			  const core::Lisp_sp& lisp )
{
    Vector3 vab = (va - vb);
    Vector3 vcb = (vc - vb);
    Vector3 vdc = (vd - vc);
    Vector3 vacCross = (vab.crossProduct(vcb)).normalized(lisp);
    LOG(BF("vacCross = %lf,%lf,%lf") % (vacCross.getX()) % (vacCross.getY()) % (vacCross.getZ() ) );
    Vector3 vdcCross = (vdc.crossProduct(vcb)).normalized(lisp);
    LOG(BF("vdcCross = %lf,%lf,%lf") % (vdcCross.getX()) % (vdcCross.getY()) % (vdcCross.getZ() ) );
    Vector3 vCross = vacCross.crossProduct(vdcCross);
    LOG(BF("vCross = %lf,%lf,%lf") % (vCross.getX()) % (vCross.getY()) % (vCross.getZ() ) );
    double dih = acos(vacCross.dotProduct(vdcCross));
    LOG(BF("dih = %lf") % (dih ) );
    double sgn = (vCross.dotProduct(vcb))<0.0?-1.0:+1.0;
//    if ( enantiomer ) return -dih*sgn;
    return dih*sgn;
}




Vector3	buildOrigin()
{
    return Vector3(0.0,0.0,0.0);
}

Vector3	buildUsingBond( double distance, const Vector3& vb )
{
Vector3	vTarget;
    vTarget = Vector3(distance,0.0,0.0);
    return vTarget+vb;
}

//! Build a vector at distance from vb and angle from v
Vector3 buildUsingBondAngle( double distance, const Vector3& vb,
			     double angle, const Vector3& va, const core::Lisp_sp& lisp )
{
Vector3 vd, vdn, vr;
double	ca, sa;
    vd = va-vb;
    vdn = vd.normalized(lisp);
    ca = cos(angle);
    sa = sin(angle);
    vr = Vector3(vdn.getX()*cos(angle)+vdn.getY()*sin(angle),
		-vdn.getX()*sin(angle)+vdn.getY()*cos(angle), 0.0);
    vr = vr.multiplyByScalar(distance);
    vr = vr+vb;
    return vr;
}


Vector3 buildUsingBondAngleDihedral( double distance, const Vector3& vc,
					double angle, const Vector3& vb,
				     double dihedral, const Vector3& va, const core::Lisp_sp& lisp )
{
    Vector3 bcDir = vb-vc;
    if ( bcDir.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
    Vector3 bcDirNorm = bcDir.normalized(lisp);
    Vector3 dPosDist = bcDirNorm.multiplyByScalar(distance);
    //
    // Now find the axis around which to rotate the bond angle
    //	
    Vector3	abDir = va-vb;
    if ( abDir.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
    Vector3 abDirNorm = abDir.normalizedOrZero();
    Vector3 angleAxis = bcDirNorm.crossProduct(abDirNorm);
    if ( angleAxis.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
    Vector3 angleAxisNorm = angleAxis.normalizedOrZero();
    Matrix angleRotation;
    angleRotation.rotationAxis(angle,&angleAxisNorm);
    Vector3 dPosAngle = angleRotation.multiplyByVector3(dPosDist);

	    //
	    // Now rotate around the dihedral bond
	    //
    Matrix dihedralRotation;
    dihedralRotation.rotationAxis(-dihedral,&bcDirNorm);
    Vector3 dPosDihedral = dihedralRotation.multiplyByVector3(dPosAngle);
	    //
	    // Now translate it to atom C
	    //
    Vector3 dPos = dPosDihedral.add(vc);
    return dPos;
}


