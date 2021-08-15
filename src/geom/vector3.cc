/*
    File: vector3.cc
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


#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cando/geom/vector3.h>
#include <clasp/core/cons.h>
#include <clasp/core/serialize.h>
#include <math.h>
//#include "core/archiveNode.h"
#include <clasp/core/lispStream.h>
#include <cando/geom/matrix.h>
#include <clasp/core/wrappers.h>



Vector3 Vector3::inNanometers() const
{
    return Vector3(this->coords[0]*0.1,
		   this->coords[1]*0.1,
		   this->coords[2]*0.1);
}


void Vector3::dump()
{
    _lisp->print(BF("<%lf,%lf,%lf>") % this->coords[0] % this->coords[1] % this->coords[2] );
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
	<< " " <<this->getY()
	<< " " <<this->getZ()
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


double	Vector3::angleToVectorAboutNormal(const Vector3& toVector, const Vector3& aboutNormal)
{
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


void Vector3::fillFromCons(core::Cons_sp vals)
{
  double x = core::clasp_to_double(core::oFirst(vals).as<core::Number_O>());
  double y = core::clasp_to_double(core::oSecond(vals).as<core::Number_O>());
  double z = core::clasp_to_double(core::oThird(vals).as<core::Number_O>());
  this->coords[0] = x;
  this->coords[1] = y;
  this->coords[2] = z;
}
    



namespace geom {
CL_DEFUN double	calculateDistance( const Vector3& va,
                                   const Vector3& vb)
{
  Vector3 vc = va - vb;
  return vc.length();
}

CL_DEFUN double	calculateDistanceSquared( const Vector3& va,
                                          const Vector3& vb)
{
  Vector3 vc = va - vb;
  return vc.dotProduct(vc);
}


/*! Return the angle in radians
 */
CL_DEFUN double calculateAngle( const Vector3& va,
                                const Vector3& vb,
                                const Vector3& vc )
{
  Vector3	vab = (va-vb).normalized();
  Vector3	vcb = (vc-vb).normalized();
  double ang = acos(vab.dotProduct(vcb));
  return ang;
}


/*! Return the dihedral in radians
 */
CL_DEFUN double calculateDihedral( const Vector3& va,
                                   const Vector3& vb,
                                   const Vector3& vc,
                                   const Vector3& vd)
{
  Vector3 vab = (va - vb);
  Vector3 vcb = (vc - vb);
  Vector3 vdc = (vd - vc);
  Vector3 vacCross = (vab.crossProduct(vcb)).normalized();
  LOG(BF("vacCross = %lf,%lf,%lf") % (vacCross.getX()) % (vacCross.getY()) % (vacCross.getZ() ) );
  Vector3 vdcCross = (vdc.crossProduct(vcb)).normalized();
  LOG(BF("vdcCross = %lf,%lf,%lf") % (vdcCross.getX()) % (vdcCross.getY()) % (vdcCross.getZ() ) );
  Vector3 vCross = vacCross.crossProduct(vdcCross);
  LOG(BF("vCross = %lf,%lf,%lf") % (vCross.getX()) % (vCross.getY()) % (vCross.getZ() ) );
  double dih = acos(vacCross.dotProduct(vdcCross));
  LOG(BF("dih = %lf") % (dih ) );
  double sgn = (vCross.dotProduct(vcb))<0.0?-1.0:+1.0;
//    if ( enantiomer ) return -dih*sgn;
  return dih*sgn;
}



CL_DOCSTRING(R"doc(Return a vector (0 0 0))doc")
CL_DEFUN Vector3 geom__build_origin()
{
  return Vector3(0.0,0.0,0.0);
}


CL_DOCSTRING(R"doc(Return a vector along the x axis distance away from vb)doc")
CL_DEFUN Vector3 geom__build_using_bond( double distance, const Vector3& vb )
{
  Vector3	vTarget;
  vTarget = Vector3(distance,0.0,0.0);
  return vTarget+vb;
}

//! Build a vector at distance from vb and angle from va
CL_DEFUN Vector3 geom__build_using_bond_angle( double distance, const Vector3& vb,
                                               double angle, const Vector3& va)
{
  Vector3 vd, vdn, vr;
  double	ca, sa;
  vd = va-vb;
  vdn = vd.normalized();
  ca = cos(angle);
  sa = sin(angle);
  vr = Vector3(vdn.getX()*cos(angle)+vdn.getY()*sin(angle),
               -vdn.getX()*sin(angle)+vdn.getY()*cos(angle), 0.0);
  vr = vr.multiplyByScalar(distance);
  vr = vr+vb;
  return vr;
}


CL_DEFUN Vector3 geom__build_using_bond_angle_dihedral( double distance, const Vector3& vc,
                                                        double angle, const Vector3& vb,
                                                        double dihedral, const Vector3& va)
{
  Vector3 bcDir = vb-vc;
  if ( bcDir.length() == 0.0 ) return Vector3(0.0,0.0,0.0);
  Vector3 bcDirNorm = bcDir.normalized();
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

CL_DEFUN double geom__planeVectorAngle(double dx, double dy)
{
  double dlen = std::sqrt(dx*dx+dy*dy);
  if (dlen<SMALL_NUMBER) return 0.0;
  dx = dx/dlen;
  dy = dy/dlen;
  if (std::fabs(dx)<SMALL_NUMBER) {
    return (MY_PI/2.0) * (dy<0.0? -1.0 : 1.0);
  } else if (dy<0.0) {
    return - std::acos(dx);
  } else {
    return std::acos(dx);
  }
}


CL_DOCSTRING(R"doc(Extract a geom:vec from a nvector at the particular index.)doc")
CL_LISPIFY_NAME("vec-extract");
CL_DEFUN void geom__vec_extract(Vector3& vec, chem::NVector_sp coordinates, size_t index0)
{
  if ((index0+2)<coordinates->length()) {
  vec.set((*coordinates)[index0],
         (*coordinates)[index0+1],
         (*coordinates)[index0+2]);
  return;
  }
  SIMPLE_ERROR(BF("Out of bounds extraction of geom:vec from nvector. Trying to extract starting at %lu and the nvector length is %lu") % index0 % coordinates->length());
}

CL_DOCSTRING(R"doc(Put a geom:vec into a nvector at the particular index.)doc")
CL_LISPIFY_NAME("vec-put");
CL_DEFUN void geom__vec_put(chem::NVector_sp coordinates, const Vector3& pos, size_t index0)
{
  if ((index0+2)<coordinates->length()) {
    (*coordinates)[index0]   = pos.getX();
    (*coordinates)[index0+1] = pos.getY();
    (*coordinates)[index0+2] = pos.getZ();
  return;
  }
  SIMPLE_ERROR(BF("Out of bounds extraction of geom:vec from nvector. Trying to put starting at %lu and the nvector length is %lu") % index0 % coordinates->length());
}

CL_DOCSTRING(R"doc(Extract a geom:vec from a nvector at the particular index.)doc")
CL_LISPIFY_NAME("vec-extract-transformed");
CL_DEFUN void geom__vec_extract_transformed(Vector3& vec, chem::NVector_sp coordinates, size_t index0, const Matrix& transform)
{
  if ((index0+2)<coordinates->length()) {
    double xp, yp, zp;
    transform.transform_nvector_point(vec.getX(),vec.getY(),vec.getZ(),coordinates,index0);
    return;
  }
  SIMPLE_ERROR(BF("Out of bounds extraction of geom:vec from nvector. Trying to extract starting at %lu and the nvector length is %lu") % index0 % coordinates->length());
}


};
