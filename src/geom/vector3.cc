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
#define	DEBUG_LEVEL_FULL


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
  core::writeln_bf_stream(fmt::sprintf("<%lf,%lf,%lf>" , this->coords[0] , this->coords[1] , this->coords[2] ));
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

void Vector3::crossProductSet( const Vector3& v, Vector3& result ) const
{
double	x,y,z;

	x = this->getY()*v.getZ()-this->getZ()*v.getY();
	y = this->getZ()*v.getX()-this->getX()*v.getZ();
	z = this->getX()*v.getY()-this->getY()*v.getX();
        result.set(x,y,z);
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
    SIMPLE_ERROR(("Attempted to normalize the vector (0,0,0)"));
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
    LOG("Vector3::parseFromString got (%lf,%lf,%lf)" , x , y , z );
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
    LOG("Vector3::parseFromCharacterArray got (%lf,%lf,%lf)" , x , y , z );
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
	LOG("me.normalized{%lf,%lf,%lf} " , men.getX() , men.getY() , men.getZ() );
	Vector3 ton = toVector.normalizedOrZero();
	LOG("to.normalized{%lf,%lf,%lf} " , ton.getX() , ton.getY() , ton.getZ() );
	LOG("about_normal{%lf,%lf,%lf} " , aboutNormal.getX() , aboutNormal.getY() , aboutNormal.getZ() );
	double cosTheta = men.dotProduct(ton);
	if ( fabs(cosTheta-1.0)<0.00001 )  return 0.0;
	if ( fabs(cosTheta+1.0)<0.00001 )  return 3.14159;
	LOG("cosTheta = %lf" , cosTheta  );
	Vector3 cr = men.crossProduct(ton);
	LOG("status" );
	double dir = cr.dotProduct(aboutNormal);
	LOG("dir = %lf " , dir );
	double theta = (dir>=0.0)?acos(cosTheta):(6.2831853-acos(cosTheta));
	LOG("theta = %lf deg" , theta/0.0174533 );
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
DOCGROUP(cando);
CL_DEFUN double	calculateDistance( const Vector3& va,
                                   const Vector3& vb)
{
  Vector3 vc = va - vb;
  return vc.length();
}

DOCGROUP(cando);
CL_DEFUN double	calculateDistanceSquared( const Vector3& va,
                                          const Vector3& vb)
{
  Vector3 vc = va - vb;
  return vc.dotProduct(vc);
}


/*! Return the angle in radians
 */
DOCGROUP(cando);
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
DOCGROUP(cando);
CL_DEFUN double calculateDihedral( const Vector3& va,
                                   const Vector3& vb,
                                   const Vector3& vc,
                                   const Vector3& vd)
{
  Vector3 vab = (va - vb);
  Vector3 vcb = (vc - vb);
  Vector3 vdc = (vd - vc);
  Vector3 vacCross = (vab.crossProduct(vcb)).normalized();
  LOG("vacCross = %lf,%lf,%lf" , (vacCross.getX()) , (vacCross.getY()) , (vacCross.getZ() ) );
  Vector3 vdcCross = (vdc.crossProduct(vcb)).normalized();
  LOG("vdcCross = %lf,%lf,%lf" , (vdcCross.getX()) , (vdcCross.getY()) , (vdcCross.getZ() ) );
  Vector3 vCross = vacCross.crossProduct(vdcCross);
  LOG("vCross = %lf,%lf,%lf" , (vCross.getX()) , (vCross.getY()) , (vCross.getZ() ) );
  double dih = acos(vacCross.dotProduct(vdcCross));
  LOG("dih = %lf" , (dih ) );
  double sgn = (vCross.dotProduct(vcb))<0.0?-1.0:+1.0;
//    if ( enantiomer ) return -dih*sgn;
  return dih*sgn;
}



CL_DOCSTRING(R"dx(Return a vector (0 0 0))dx");
DOCGROUP(cando);
CL_DEFUN Vector3 geom__build_origin()
{
  return Vector3(0.0,0.0,0.0);
}


CL_DOCSTRING(R"dx(Return a vector along the x axis distance away from vb)dx");
DOCGROUP(cando);
CL_DEFUN Vector3 geom__build_using_bond( double distance, const Vector3& vb )
{
  Vector3	vTarget;
  vTarget = Vector3(distance,0.0,0.0);
  return vTarget+vb;
}

//! Build a vector at distance from vb and angle from va
DOCGROUP(cando);
CL_DEFUN Vector3 geom__build_using_bond_angle( double distance, const Vector3& vb,
                                               double angle, const Vector3& va)
{
  Vector3 vd, vdn, vr;
  vd = va-vb;
  vdn = vd.normalized();
  vr = Vector3(vdn.getX()*cos(angle)+vdn.getY()*sin(angle),
               -vdn.getX()*sin(angle)+vdn.getY()*cos(angle), 0.0);
  vr = vr.multiplyByScalar(distance);
  vr = vr+vb;
  return vr;
}


DOCGROUP(cando);
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


/*
 *      myAcos
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the acos of a number, but don't crap
 *      out on domain errors.
 */
double
myAcos( double d )
{
  if ( d >= 1.0 )
    return(0.0);
  if ( d <= -1.0 )
    return(core::numerics::pi);
  return(acos(d));
}


/*
 *      dVectorAbsAngle
 *
 *      Author: Christian Schafmeister (1991)
 *
 *      Return the absolute angle between two vectors in Radians
 *      the sign of the angle is determined by the direction of
 *      the cross product.
 */
double
dVectorAbsAngle( const Vector3& vX, const Vector3& vY, const Vector3& vRef )
{
  Vector3  vT1, vT2, vT;
  double  dAngle;

  vT1 = vX.normalized();
  vT2 = vY.normalized();

  vT = vT1.crossProduct(vT2);
  dAngle = myAcos(vT1.dotProduct(vT2));
  if ( vT.dotProduct(vRef) < 0.0 ) dAngle = -dAngle;
  return(dAngle);
}



/*
 *      zvZMatrixCalculatePositionFromAngles
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Use NEWTON-RAPHSON method for finding the coordinate for the
 *      vector(vC) which is dAngleA from the vector vA (on the X axis)
 *	and dAngleB from a vector (vB) which lies in the XY plane
 *	dAngleC from the X axis.
 *
 *      The point is dBond from the origin.
 *
 */
#define VERYSMALL       0.000000000001
#define MAXNEWTONSTEPS  20
const Vector3 vXAxis = { 1.0, 0.0, 0.0 };
const Vector3 vYAxis = { 0.0, 1.0, 0.0 };
const Vector3 vZAxis = { 0.0, 0.0, 1.0 };

Vector3 zvZMatrixCalculatePositionFromAngles( double dAngleA, double dAngleB,
                                              double dAngleC, double dBond )
{
  int             iCount;
  double		dCosA, dSinA;
  double		dCosB;
  double		dCosC, dSinC;
  double		dCosX, dSinX;
  double		dX, dXNew;
  double		dF1, dF2;
  Vector3		vNew;

  dCosA = cos(dAngleA);
  dSinA = sin(dAngleA);
  dCosB = cos(dAngleB);

  dCosC = cos(dAngleC);
  dSinC = sin(dAngleC);

		/* The idea is to minimize the function: */
		/* E = ( DOT(vC,vB) - cos(dAngleB) )^2 */
		/* using NEWTONS method */
		/* The vector vC is constrained to make the angle */
		/* dAngleA with vA */
		/* The vector vC makes the angle (dX) with the XY plane */
		/* and the parameter that is optimized is dX */


		/* A reasonable starting point */

  dX = dAngleB;
  iCount = 0;
  while ( iCount <MAXNEWTONSTEPS ) {

    dCosX = cos(dX);
    dSinX = sin(dX);

    dF1 = -2*dSinA*dSinC*(-dCosB + dCosA*dCosC +
                          dCosX*dSinA*dSinC)*dSinX;

    dF2 = -2.0*dCosX*dSinA*dSinC*
              (-dCosB + dCosA*dCosC + dCosX*dSinA*dSinC) +
              2.0*dSinA*dSinA*dSinC*dSinC*dSinX*dSinX;

    LOG( "Iteration %d dF1=%lf  dF2=%lf  dB=%lf\n",
         iCount, dF1, dF2, dX );
    if ( fabs(dF1) < VERYSMALL*10.0 ) break;
    if ( fabs(dF2) < VERYSMALL ) {
      SIMPLE_ERROR( "Could not optimize! dF1 = %lf, dF2 = %lfdX = %lf steps=%d", dF1, dF2, dX, iCount );
    }
    dXNew = dX - dF1/dF2;
    if ( fabs(dXNew - dX) < VERYSMALL ) break;
    dX = dXNew;
    iCount++;
  }

#ifdef  DEBUG
  if ( iCount > MAXNEWTONSTEPS )
    DDEBUG( ("Exceeded maximum number of Newton Raphson steps: %d\n",
             MAXNEWTONSTEPS) );
#endif

                /* Generate new coordinate */

  vNew.set( dBond*cos(dAngleA),
            dBond*sin(dAngleA)*cos(dX),
            dBond*sin(dAngleA)*sin(dX) );
  return(vNew);
}



/*
 *      ZMatrixBondTwoAnglesOrientation
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Build the external coordinate for the atom when
 *      the orientation, a bond length and two angles are supplied.
 *      The orientation is a positive or negative number which specifies
 *      the orientation of the new position.  It is calculated by:
 *              a=crossProduct( vPAtomA-vPCenter, vPAtomB-vPCenter );
 *              orientation = dotProduct( vPPos-vPCenter, a );
 *
 *      vPAtomC points to the position of the central atom.
 *
 */

CL_DEFUN Vector3 geom__build_using_bond_two_angles_orientation(
    const Vector3& vCenter,
    const Vector3& vAtomA,
    const Vector3& vAtomB,
    double dBond,
    double dAngleA,
    double dAngleB,
    double dOrient )
{
  Matrix          mT, mT1, mT2, mTX, mTY, mTZ, mTT;
  double          dAngleX, dAngleY, dAngleZ;
  double          dAngle;
  Vector3         vTrans, vTempAC, vTempBC, vTempXZ, vNew, vLab;

                /* The procedure for finding the the coordinate is: */
                /* Translate vCenter to the origin -> A'-C' */
                /* Find angle between PROJ((A'-C'),YZ plane) & Y axis */
                /* Rotate into XZ plane */
                /* Find angle between (A''-C'') and X axis */
                /* Rotate onto X axis */
                /* Find angle between PROJ((B'''-C'''),YZ plane) and Y axis */
                /* Rotate onto XY plane */
                /* Calculate coordinates in 3Space */
                /* Apply the reverse transformation to the new point */
                /* Actually, all that is done is the elements for the */
                /* forward transformations are calculated then used */
                /* to generate an inverse transform matrix */

  vTrans = vCenter;
  vTempAC = vAtomA-vCenter;
  vTempBC = vAtomB-vCenter; // vVectorSub( vPAtomB, vPAtomC );
  LOG( "AC= %s\n", vTempAC.asString() );
  LOG( "BC= %s\n", vTempBC.asString() );
  vTempXZ = vTempAC;
  vTempXZ.getY() = 0.0;
  LOG( "XZ= %s\n", vTempXZ.asString() );
  if ( vTempXZ.length() != 0.0 ) {
    dAngleY = dVectorAbsAngle(vTempXZ, vXAxis, vYAxis );
  } else dAngleY = 0.0;
  LOG( "dAngleY = %lf\n", dAngleY );
  mT.rotationY(-dAngleY );
  LOG( "Rotated around Y\n" );
  vTempAC = mT.multiplyByVector3(vTempAC);
  vTempBC = mT.multiplyByVector3(vTempBC);
  LOG( "New AC= %s\n", vTempAC.asString() );
  LOG( "New BC= %s\n", vTempBC.asString() );

  dAngleZ = dVectorAbsAngle( vTempAC, vXAxis, vZAxis );
  mT.rotationZ( -dAngleZ );
  vTempBC = mT.multiplyByVector3(vTempBC);
#ifdef DEBUG
  vTempAC = mT.multiplyByVector3(vTempAC);
#endif
  LOG( "Rotated around Z\n" );
  LOG( "New AC= %s\n", vTempAC.asString() );
  LOG( "New BC= %s\n", vTempBC.asString() );

  vTempBC.getX() = 0.0;

  dAngleX = dVectorAbsAngle( vTempBC, vYAxis, vXAxis );

                /* Build the transformation matrix to convert from */
                /* lab coordinates to molecule coordinates in mT*/

  mTX.rotationX( dAngleX );
  mTZ.rotationZ( dAngleZ ); // MatrixZRotate( mTZ, dAngleZ );
  mTY.rotationY( dAngleY ); // MatrixYRotate( mTY, dAngleY );
  mTT.translate( vTrans ); // MatrixTranslate( mTT, dVX(&vTrans), dVY(&vTrans), dVZ(&vTrans) );
  mT1 = mTZ.multiplyByMatrix(mTX); // MatrixMultiply mT1, mTZ, mTX );
  mT2 = mTY.multiplyByMatrix(mT1); // MatrixMultiply( mT2, mTY, mT1 );
  mT = mTT.multiplyByMatrix(mT2);  // MatrixMultiply( mT, mTT, mT2 );

                /* Calculate coordinates of new atom */
  dAngle = calculateAngle( vAtomA, vCenter, vAtomB ); // dVectorAtomAngle( vPAtomA, vPAtomC, vPAtomB );
  vLab = zvZMatrixCalculatePositionFromAngles( dAngleA, dAngleB, dAngle, dBond );

  if ( dOrient != 0.0 ) {
    vLab.getZ() = dOrient * vLab.getZ(); // VectorSetZ( &vLab, dOrient*dVZ(&vLab) );
  }

                /* If there is no chirality defined yet then just */
                /* leave it the way it is */

  vNew = mT.multiplyByVector3(vLab);
  LOG( "ZMatrix2Angle:  %s\n", vNew.asString() );
  return vNew;
}

DOCGROUP(cando);
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


CL_DOCSTRING(R"dx(Extract a geom:vec from a nvector at the particular index.)dx");
CL_LISPIFY_NAME("vec-extract");
DOCGROUP(cando);
CL_DEFUN void geom__vec_extract(Vector3& vec, chem::NVector_sp coordinates, size_t index0)
{
  if ((index0+2)<coordinates->length()) {
  vec.set((*coordinates)[index0],
         (*coordinates)[index0+1],
         (*coordinates)[index0+2]);
  return;
  }
  SIMPLE_ERROR(("Out of bounds extraction of geom:vec from nvector. Trying to extract starting at %lu and the nvector length is %lu") , index0 , coordinates->length());
}

CL_DOCSTRING(R"dx(Put a geom:vec into a nvector at the particular index.)dx");
CL_LISPIFY_NAME("vec-put");
DOCGROUP(cando);
CL_DEFUN void geom__vec_put(chem::NVector_sp coordinates, const Vector3& pos, size_t index0)
{
  if ((index0+2)<coordinates->length()) {
    (*coordinates)[index0]   = pos.getX();
    (*coordinates)[index0+1] = pos.getY();
    (*coordinates)[index0+2] = pos.getZ();
  return;
  }
  SIMPLE_ERROR(("Out of bounds extraction of geom:vec from nvector. Trying to put starting at %lu and the nvector length is %lu") , index0 , coordinates->length());
}

CL_DOCSTRING(R"dx(Extract a geom:vec from a nvector at the particular index.)dx");
CL_LISPIFY_NAME("vec-extract-transformed");
DOCGROUP(cando);
CL_DEFUN void geom__vec_extract_transformed(Vector3& vec, chem::NVector_sp coordinates, size_t index0, const Matrix& transform)
{
  if ((index0+2)<coordinates->length()) {
    transform.transform_nvector_point(vec.getX(),vec.getY(),vec.getZ(),coordinates,index0);
    return;
  }
  SIMPLE_ERROR(("Out of bounds extraction of geom:vec from nvector. Trying to extract starting at %lu and the nvector length is %lu") , index0 , coordinates->length());
}


};
