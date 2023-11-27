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
#include <clasp/core/evaluator.h>
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
#include <clasp/core/array.h>


Vector3::Vector3(core::SimpleVector_float_sp svf, size_t index) {
  if (index+2<svf->length()) {
    this->coords[0] = (*svf)[index];
    this->coords[1] = (*svf)[index+1];
    this->coords[2] = (*svf)[index+2];
    return;
  }
  SIMPLE_ERROR("The index {} is out of range for a 3D vector within the array of length {}",
               index, svf->length());
}

Vector3 Vector3::inNanometers() const
{
    return Vector3(this->coords[0]*0.1,
		   this->coords[1]*0.1,
		   this->coords[2]*0.1);
}


void Vector3::dump()
{
  core::clasp_write_string(fmt::format("<{},{},{}>\n" , this->coords[0] , this->coords[1] , this->coords[2] ));
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
    SIMPLE_ERROR("Attempted to normalize the vector (0,0,0)");
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
    LOG("Vector3::parseFromString got ({},{},{})" , x , y , z );
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
    LOG("Vector3::parseFromCharacterArray got ({},{},{})" , x , y , z );
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
	LOG("me.normalized{{},{},{}} " , men.getX() , men.getY() , men.getZ() );
	Vector3 ton = toVector.normalizedOrZero();
	LOG("to.normalized{{},{},{}} " , ton.getX() , ton.getY() , ton.getZ() );
	LOG("about_normal{{},{},{}} " , aboutNormal.getX() , aboutNormal.getY() , aboutNormal.getZ() );
	double cosTheta = men.dotProduct(ton);
	if ( fabs(cosTheta-1.0)<0.00001 )  return 0.0;
	if ( fabs(cosTheta+1.0)<0.00001 )  return 3.14159;
	LOG("cosTheta = {}" , cosTheta  );
	Vector3 cr = men.crossProduct(ton);
	LOG("status" );
	double dir = cr.dotProduct(aboutNormal);
	LOG("dir = {} " , dir );
	double theta = (dir>=0.0)?safe_acos(cosTheta):(6.2831853-safe_acos(cosTheta));
	LOG("theta = {} deg" , theta/0.0174533 );
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


/*
 *      safe_acos
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the acos of a number, but don't crap
 *      out on domain errors.
 */
double safe_acos( double d )
{
  double ac = std::acos(d);
  unlikely_if (isnan(ac)) {
    unlikely_if ( d >= 1.0 ) return(0.0);
    unlikely_if ( d <= -1.0 ) return(core::numerics::pi);
    SIMPLE_ERROR("Could not calculate acos of {}", d);
  }
  return ac;
}

float safe_acos( float d )
{
  float ac = std::acos(d);
  unlikely_if (isnan(ac)) {
    unlikely_if ( d >= 1.0 ) return(0.0);
    unlikely_if ( d <= -1.0 ) return(core::numerics::pi);
    SIMPLE_ERROR("Could not calculate acos of {}", d);
  }
  return ac;
}

/*
 *      safe_asin
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return the asin of a number, but don't crap
 *      out on domain errors.
 */
double safe_asin( double d )
{
  double as = std::asin(d);
  unlikely_if(isnan(as)) {
    unlikely_if ( d >= 1.0 ) return(core::numerics::pi/2.0);
    unlikely_if ( d <= -1.0 ) return(-core::numerics::pi/2.0);
    SIMPLE_ERROR("Could not calculate asin of {}", d);
  }
  return as;
}

float safe_asin( float d )
{
  float as = std::asin(d);
  unlikely_if(isnan(as)) {
    unlikely_if ( d >= 1.0 ) return(core::numerics::pi/2.0);
    unlikely_if ( d <= -1.0 ) return(-core::numerics::pi/2.0);
    SIMPLE_ERROR("Could not calculate asin of {}", d);
  }
  return as;
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
  double   dAngle;

  vT1 = vX.normalized();
  vT2 = vY.normalized();

  vT = vT1.crossProduct(vT2);
  dAngle = safe_acos(vT1.dotProduct(vT2));
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
  double		dCosB, dSinB;
  double		dCosC, dSinC;
  double		dCosX, dSinX;
  double		dX, dXNew;
  double		dF1, dF2;
  Vector3		vNew;

  dCosA = cos(dAngleA);
  dSinA = sin(dAngleA);
  dCosB = cos(dAngleB);
  dSinB = sin(dAngleB);
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

    LOG( "Iteration {} dF1={}  dF2={}  dB={}\n",
         iCount, dF1, dF2, dX );
    if ( fabs(dF1) < VERYSMALL*10.0 ) break;
    if ( fabs(dF2) < VERYSMALL ) {
      SIMPLE_ERROR( "Could not optimize! dF1 = {}, dF2 = {}dX = {} steps={}", dF1, dF2, dX, iCount );
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
  LOG( "AC= {}\n", vTempAC.asString() );
  LOG( "BC= {}\n", vTempBC.asString() );
  vTempXZ = vTempAC;
  vTempXZ.getY() = 0.0;
  LOG( "XZ= {}\n", vTempXZ.asString() );
  if ( vTempXZ.length() != 0.0 ) {
    dAngleY = dVectorAbsAngle(vTempXZ, vXAxis, vYAxis );
  } else dAngleY = 0.0;
  LOG( "dAngleY = {}\n", dAngleY );
  mT.rotationY(-dAngleY );
  LOG( "Rotated around Y\n" );
  vTempAC = mT.multiplyByVector3(vTempAC);
  vTempBC = mT.multiplyByVector3(vTempBC);
  LOG( "New AC= {}\n", vTempAC.asString() );
  LOG( "New BC= {}\n", vTempBC.asString() );

  dAngleZ = dVectorAbsAngle( vTempAC, vXAxis, vZAxis );
  mT.rotationZ( -dAngleZ );
  vTempBC = mT.multiplyByVector3(vTempBC);
#ifdef DEBUG
  vTempAC = mT.multiplyByVector3(vTempAC);
#endif
  LOG( "Rotated around Z\n" );
  LOG( "New AC= {}\n", vTempAC.asString() );
  LOG( "New BC= {}\n", vTempBC.asString() );

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
  LOG( "ZMatrix2Angle:  {}\n", vNew.asString() );
  return vNew;
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
  SIMPLE_ERROR("Out of bounds extraction of geom:vec from nvector. Trying to extract starting at {} and the nvector length is {}" , index0 , coordinates->length());
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
  SIMPLE_ERROR("Out of bounds extraction of geom:vec from nvector. Trying to put starting at {} and the nvector length is {}" , index0 , coordinates->length());
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
  SIMPLE_ERROR("Out of bounds extraction of geom:vec from nvector. Trying to extract starting at {} and the nvector length is {}" , index0 , coordinates->length());
}

CL_DEFUN core::T_sp geom__vecreal_type() {
#if VECREAL==VECREAL_DOUBLE
  return cl::_sym_double_float;
#else
  return cl::_sym_single_float;
  #endif
}

CL_DEFUN core::T_sp geom__vecreal(core::T_sp val) {
#if VECREAL==VECREAL_DOUBLE
  if (val.single_floatp()) {
    return mk_double_float(val.unsafe_single_float());
  } else if (gc::IsA<core::DoubleFloat_sp>(val)) {
    return val;
  } else {
    TYPE_ERROR(val,core::Cons_O::createList(cl::_sym_double_float, cl::_sym_single_float));
  }
#else
  if (val.single_floatp()) {
    return val;
  } else if (gc::IsA<core::DoubleFloat_sp>(val)) {
    return core::make_single_float(gc::As<core::DoubleFloat_sp>(val)->get());
  } else {
    TYPE_ERROR(val,core::Cons_O::createList(cl::_sym_double_float, cl::_sym_single_float));
  }
#endif
}
Vector3 transform_rotor3( float r_scalar, float r_xy, float r_yz, float r_zx, const Vector3& v)
{
    const float S_x = r_scalar*v.getX() + r_xy*v.getY() - r_zx*v.getZ();
    const float S_y = r_scalar*v.getY() - r_xy*v.getX() + r_yz*v.getZ();
    const float S_z = r_scalar*v.getZ() - r_yz*v.getY() + r_zx*v.getX();
    const float S_xyz = r_xy*v.getZ() + r_yz*v.getX() + r_zx*v.getY();

    float result_x = S_x*r_scalar +   S_y*r_xy + S_xyz*r_yz -   S_z*r_zx;
    float result_y = S_y*r_scalar -   S_x*r_xy +   S_z*r_yz + S_xyz*r_zx;
    float result_z = S_z*r_scalar + S_xyz*r_xy -   S_y*r_yz +   S_x*r_zx;
    return Vector3( result_x, result_y, result_z );
}

};
