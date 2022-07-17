/*
    File: coordSys.cc
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

//
// coordSys.cc
//
#include <clasp/core/common.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/coordSys.h>
#include <clasp/core/numerics.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/omatrix.h>
//#include "core/archiveNode.h"
#include <cando/geom/color.h>
#include <cando/chem/atom.h>
#include <clasp/core/wrappers.h>

#define	PI		3.1415926535897932
#define	DEGTORAD	0.0174533
#define	PARALLEL_CUTOFF	0.00001


namespace chem {

void	CoordinateSystem_O::initialize()
{
    this->Base::initialize();
    this->canonical();
}


CoordinateSystem_O::CoordinateSystem_O(const CoordinateSystem_O& orig) : core::CxxObject_O(orig)
{
    this->origin = orig.origin;
    this->x = orig.x;
    this->y = orig.y;
    this->z = orig.z;
}


//
//	matrixFromCanonical
//
//	Generate the matrix for transformation
//	from the canonical coordinate system to this one
//
//	This means that a point at 0,0,0 will be transformed to
//	this coordinate systems origin.
//	A point at 1,0,0 will be transformed to a point
//	one unit along this ones X-axis
//
CL_LISPIFY_NAME("matrixFromCanonical");
CL_DEFMETHOD Matrix	CoordinateSystem_O::matrixFromCanonical()
{
Matrix	mt, m;

    mt.at(0, 0) = this->x.getX();
    mt.at(1, 0) = this->x.getY();
    mt.at(2, 0) = this->x.getZ();

    mt.at(0, 1) = this->x.getX();
    mt.at(1, 1) = this->x.getY();
    mt.at(2, 1) = this->x.getZ();

    mt.at(0, 2) = this->x.getX();
    mt.at(1, 2) = this->x.getY();
    mt.at(2, 2) = this->x.getZ();

    m.translate(this->origin);
    m = m*mt;
    return m;
}


//
//	matrixToCanonical
//
//	Generate the matrix for transformation
//	to the canonical coordinate system from this one
//
//	This means that a point on this coordinate systems origin will be
//	transformed to 0,0,0.  A point one unit along this coordinate systems
//	X-axis will be transformed to 1,0,0
//
CL_LISPIFY_NAME("matrixToCanonical");
CL_DEFMETHOD Matrix	CoordinateSystem_O::matrixToCanonical()
{
Matrix  mt, m, mn;
Vector3 vt;
    m.at(0, 0) = this->x.getX();
    m.at(0, 1) = this->x.getY();
    m.at(0, 2) = this->x.getZ();

    m.at(1, 0) = this->y.getX();
    m.at(1, 1) = this->y.getY();
    m.at(1, 2) = this->y.getZ();

    m.at(2, 0) = this->z.getX();
    m.at(2, 1) = this->z.getY();
    m.at(2, 2) = this->z.getZ();

    vt = this->origin.multiplyByScalar(-1.0);
    mt.translate(vt);
    LOG("mt = %s" , (mt.asString().c_str() ) );
    mn = m*mt;
    LOG("mn = %s" , (mn.asString().c_str() ) );
    return mn;
}


//
//	canonical
//
//	Set (this) to the canonical coordinate set.
//	origin is 0,0,0
//	x is 1,0,0
//	y is 0,1,0
//	z is 0,0,1
CL_LISPIFY_NAME("canonical");
CL_DEFMETHOD void CoordinateSystem_O::canonical()
{
    this->origin = Vector3(0.0,0.0,0.0);
    this->x = Vector3(1.0,0.0,0.0);
    this->y = Vector3(0.0,1.0,0.0);
    this->z = Vector3(0.0,0.0,1.0);
}


//
//	defineForAtoms
//
//	Build a coordinate set using the three atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the z axis is orthogonal to (axy) and x axis
//	the y axis is orthogonal to the new x and z axes
CL_LISPIFY_NAME("defineForAtoms");
CL_DOCSTRING(R"dx(Define a coordinate system using the origin, x-axis and another vector in the xy plane)dx")
CL_DEFMETHOD void	CoordinateSystem_O::defineForAtoms( Atom_sp aorigin, Atom_sp ax, Atom_sp axy )
{
    Vector3	vo, vx, vxy, vz, vy;
    vo = aorigin->getPosition();
    LOG("vo=%lf,%lf,%lf" , vo.getX() , vo.getY() , vo.getZ() );
    vx = (ax->getPosition() - vo).normalized();
    LOG("vx=%lf,%lf,%lf" , vx.getX() , vx.getY() , vx.getZ() );
    vxy = (axy->getPosition() - vo).normalized();
    LOG("vxy=%lf,%lf,%lf" , vxy.getX() , vxy.getY() , vxy.getZ() );
    vz = (vx.crossProduct(vxy)).normalized();
    LOG("vz=%lf,%lf,%lf" , vz.getX() , vz.getY() , vz.getZ() );
    vy = (vz.crossProduct(vx)).normalized();
    LOG("vy=%lf,%lf,%lf" , vy.getX() , vy.getY() , vy.getZ() );
    this->origin = vo;
    this->x = vx;
    this->y = vy;
    this->z = vz;
}


//
//	defineForAtomVectors
//
//	Build a coordinate set using the three atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the z axis is orthogonal to (axy) and x axis
//	the y axis is orthogonal to the new x and z axes
CL_LISPIFY_NAME("defineForAtomVectors");
CL_DOCSTRING(R"dx(Define a coordinate system using the three vectors. vo is the origin
vax is along the x axis and vaxy is in the xy-plane.)dx")
CL_DEFMETHOD void	CoordinateSystem_O::defineForAtomVectors( const Vector3& vo,
						  const Vector3& vax,
						  const Vector3& vaxy )
{
    Vector3	vx, vxy, vz, vy;
vx = (vax.sub(vo)).normalized();
vxy = (vaxy.sub(vo)).normalized();
    LOG("vo=%lf,%lf,%lf" , vo.getX() , vo.getY() , vo.getZ() );
    LOG("vx=%lf,%lf,%lf" , vx.getX() , vx.getY() , vx.getZ() );
    LOG("vxy=%lf,%lf,%lf" , vxy.getX() , vxy.getY() , vxy.getZ() );

    vz = (vx.crossProduct(vxy)).normalized();
    LOG("vz=%lf,%lf,%lf" , vz.getX() , vz.getY() , vz.getZ() );
    vy = (vz.crossProduct(vx)).normalized();
    LOG("vy=%lf,%lf,%lf" , vy.getX() , vy.getY() , vy.getZ() );
    this->origin = vo;
    this->x = vx;
    this->y = vy;
    this->z = vz;
}

//
//	defineForTwoAtomVectors
//
//	Build a coordinate set using the two atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the y axis is whereever it ends up
//	the z axis is orthogonal to y and x axis
void	CoordinateSystem_O::defineForTwoAtomVectors( const Vector3& vo,
						     const Vector3& vax )
{
    Vector3	vx, vz, vy;
vx = (vax.sub(vo)).normalized();
    LOG("vo=%lf,%lf,%lf" , vo.getX() , vo.getY() , vo.getZ() );
    LOG("vx=%lf,%lf,%lf" , vx.getX() , vx.getY() , vx.getZ() );
    Vector3 vOffset;
    vOffset.set(0.0,0.5,0.0);
    vy = (vx.add(vOffset)).normalized();
    	// Check if vy is a little offset from vx
    vz = ( vx.crossProduct(vy) );
    if ( vz.dotProduct(vz) < 0.0001 )
    {
	vOffset.set(0.5,0.0,0.0);
	vy = (vx.add(vOffset)).normalized();
	vz = ( vx.crossProduct(vy) ).normalized();
    } else
    {
	vz = vz.normalized();
    }
    vy = (vz.crossProduct(vx)).normalized();
    LOG("vy=%lf,%lf,%lf" , vy.getX() , vy.getY() , vy.getZ() );
    LOG("vz=%lf,%lf,%lf" , vz.getX() , vz.getY() , vz.getZ() );
    this->origin = vo;
    this->x = vx;
    this->y = vy;
    this->z = vz;
}

//
//	defineForOneAtomVectors
//
//	Build a coordinate set using the two atoms
//	the origin is on atom (aorigin)
//	the x axis is along (ax-aorigin)
//	the y axis is whereever it ends up
//	the z axis is orthogonal to y and x axis
void	CoordinateSystem_O::defineForOneAtomVector( const Vector3& vo )
{
    Vector3	vx, vz, vy;
    this->origin = vo;
    this->x.set(1.0,0.0,0.0);
    this->y.set(0.0,1.0,0.0);
    this->z.set(0.0,0.0,1.0);
}



//
//	defineForVectorsOriginXDirYDir
//
//	Build a coordinate set using the origin and the direction
//	of the X-axis and Y-axis.  The Z-axis is calculated from
//	the cross-product of xDir and yDir.
//
CL_LISPIFY_NAME("defineForVectorsOriginXDirYDir");
CL_DEFMETHOD void	CoordinateSystem_O::defineForVectorsOriginXDirYDir(
    const Vector3& orig,
    const Vector3& xDir,
    const Vector3& yDir )
{
    Vector3	vzt;
    this->origin = orig;
    this->x = xDir.normalized();
    this->y = yDir.normalized();
    vzt = (this->x).crossProduct(this->y);
    this->z = vzt.normalized();
}


//
//	defineForVectorsOriginXDirXYPlane
//
//	Build a coordinate set using the origin and the direction
//	of the X-axis and a vector in the XY plane.
//	The Z-axis is calculated from the cross-product of xDir and xyPlane.
//
CL_LISPIFY_NAME("defineForVectorsOriginXDirXYPlane");
CL_DEFMETHOD void	CoordinateSystem_O::defineForVectorsOriginXDirXYPlane(
    const Vector3& orig,
    const Vector3& xDir,
    const Vector3& xyPlane )
{
    Vector3	vzt,vxy,vyt;
    this->origin = orig;
    this->x = xDir.normalized();
    vxy = xyPlane.normalized();
    vzt = (this->x).crossProduct(vxy);
    this->z = vzt.normalized();
    vyt = (this->z).crossProduct(this->x);
    this->y = vyt.normalized();
}

//
//	defineForVectorsOriginZDirXZPlane
//
//	Build a coordinate set using the origin and the direction
//	of the Z-axis and a vector in the XZ plane.
//	The Y-axis is calculated from the cross-product of xDir and xyPlane.
//
CL_LISPIFY_NAME("defineForVectorsOriginZDirXZPlane");
CL_DEFMETHOD void	CoordinateSystem_O::defineForVectorsOriginZDirXZPlane(
    const Vector3& orig,
    const Vector3& zDir,
    const Vector3& xzPlane )
{
    Vector3		vyt,vxz,vxt;
    this->origin = orig;
    this->z = zDir.normalized();
    vxz = xzPlane.normalized();
    vyt = (this->z).crossProduct(vxz);
    this->y = vyt.normalized();
    vxt = (this->y).crossProduct(this->z);
    this->x = vxt.normalized();
}




//
//	defineForVectorsOriginXDirZDir
//
//	Build a coordinate set using the origin and the direction
//	of the X-axis and Z-axis.  The Y-axis is calculated from
//	the cross-product of xDir and zDir.
//
CL_LISPIFY_NAME("defineForVectorsOriginXDirZDir");
CL_DEFMETHOD void	CoordinateSystem_O::defineForVectorsOriginXDirZDir(
    const Vector3& orig,
    const Vector3& xDir,
    const Vector3& zDir )
{
    Vector3	vyt;
    this->origin = orig;
    this->x = xDir.normalized();
    this->z = zDir.normalized();
    vyt = (this->z).crossProduct(this->x);
    this->y = vyt.normalized();
}


CL_LISPIFY_NAME("defineRandom");
CL_DEFMETHOD void	CoordinateSystem_O::defineRandom( )
{_OF();
double x,y,z;
#define	RANGE	20.0
x = core::randomNumber01()*RANGE - (RANGE/2.0);
y = core::randomNumber01()*RANGE - (RANGE/2.0);
z = core::randomNumber01()*RANGE - (RANGE/2.0);
    this->origin.set(x,y,z);
    while (1) {
	x = core::randomNumber01()*2.0-1.0;
	y = core::randomNumber01()*2.0-1.0;
	z = core::randomNumber01()*2.0-1.0;
	this->x.set(x,y,z);
	if ( this->x.length() > 0.00001 ) break;
    }
    this->x = this->x.normalized();
    while (1) {
	x = core::randomNumber01()*2.0-1.0;
	y = core::randomNumber01()*2.0-1.0;
	z = core::randomNumber01()*2.0-1.0;
	this->y.set(x,y,z);
	if ( this->y.length() > 0.00001 ) break;
    }
    this->z = this->x.crossProduct(this->y);
    this->z = this->z.normalized();
    this->y = this->z.crossProduct(this->x);
}


CL_LISPIFY_NAME("defineFromCoordinateSystem");
CL_DEFMETHOD void	CoordinateSystem_O::defineFromCoordinateSystem( CoordinateSystem_sp c )
{
    this->origin = c->getOrigin();
    this->x = c->getX();
    this->y = c->getY();
    this->z = c->getZ();
}


CL_LISPIFY_NAME("defineCanonical");
CL_DEFMETHOD void	CoordinateSystem_O::defineCanonical( )
{
    this->origin.set(0.0,0.0,0.0);
    this->x.set(1.0,0.0,0.0);
    this->y.set(0.0,1.0,0.0);
    this->z.set(0.0,0.0,1.0);
}








//
//	matrixForTransformTo
//
//	Return the matrix to transform from (this) to (cs)
Matrix	CoordinateSystem_O::matrixForTransformTo(CoordinateSystem_sp cs )
{
Matrix	m1, m2;

    m1 = this->matrixToCanonical();
    m2 = cs->matrixFromCanonical();
    return(m1*m2);
}




//
//	transformWithMatrix
//
//	Transform this coordinate system with the matrix
//
CL_LISPIFY_NAME("transformWithMatrix");
CL_DEFMETHOD void	CoordinateSystem_O::transformWithMatrix(const Matrix& m )
{
    Vector3	vox, voy, voz;

    vox = m*(this->x + this->origin);
    voy = m*(this->y + this->origin);
    voz = m*(this->z + this->origin);
    this->origin = m*this->origin;
    this->x = vox - this->origin;
    this->y = voy - this->origin;
    this->z = voz - this->origin;
}

//
//	transformWithMatrixPointer
//
//	Transform this coordinate system with the matrix
//
void	CoordinateSystem_O::transformWithMatrixPointer( Matrix* m )
{
    Vector3	vox, voy, voz;

    vox = (*m)*(this->x + this->origin);
    voy = (*m)*(this->y + this->origin);
    voz = (*m)*(this->z + this->origin);
    this->origin = (*m)*this->origin;
    this->x = vox - this->origin;
    this->y = voy - this->origin;
    this->z = voz - this->origin;
}

#if 0
//
//	renderXml
//
//	Render the coordinate system in XML format
//
adapt::QDomNode_sp	CoordinateSystem_O::renderXml( double width,
							double axisLength)
{
IMPLEMENT_ME();
adapt::QDomNode_sp	graalphacs, line;
#if 0 //[

    graalphacs = xmlBlock();
    line = xmlLine( XML_RED, (int)width, this->origin,
				this->origin+this->x*axisLength );
    graalphacs->addChild(line);
    line = xmlLine( XML_GREEN, (int)width, this->origin,
				this->origin+this->y*axisLength );
    graalphacs->addChild(line);
    line = xmlLine( XML_BLUE, (int)width, this->origin,
				this->origin+this->z*axisLength );
    graalphacs->addChild(line);
#endif //]
    return graalphacs;
}
#endif



//
//	dump
//
//	Dump the coordinate system
//
CL_LISPIFY_NAME("dump");
CL_DEFMETHOD void	CoordinateSystem_O::dump()
{
  core::writeln_bf_stream(fmt::sprintf("origin: ( %lf, %lf, %lf )" , this->origin.getX() , this->origin.getY() , this->origin.getZ() ));
  core::writeln_bf_stream(fmt::sprintf("x-axis: ( %lf, %lf, %lf )" , this->x.getX() , this->x.getY() , this->x.getZ() ));
  core::writeln_bf_stream(fmt::sprintf("y-axis: ( %lf, %lf, %lf )" , this->y.getX() , this->y.getY() , this->y.getZ() ));
  core::writeln_bf_stream(fmt::sprintf("z-axis: ( %lf, %lf, %lf )" , this->z.getX() , this->z.getY() , this->z.getZ() ));
}

#ifdef XML_ARCHIVE
void	CoordinateSystem_O::archive(core::ArchiveP node)
{
    node->archivePlainObject<Vector3>( "origin", "Vector", this->origin );
    node->archivePlainObject<Vector3>( "x", "Vector", this->x );
    node->archivePlainObject<Vector3>( "y", "Vector", this->y );
    node->archivePlainObject<Vector3>( "z", "Vector", this->z );
}
#endif


CL_DEFMETHOD CoordinateSystem_sp	CoordinateSystem_O::copyCoordinateSystem()
{
  auto  n  = gctools::GC<CoordinateSystem_O>::copy( *this); // = RP_Copy<CoordinateSystem_O>(this);
    return n;
}


CL_LISPIFY_NAME("sameAs");
CL_DEFMETHOD bool	CoordinateSystem_O::sameAs(CoordinateSystem_sp c)
{
    if ( !this->origin.sameAs(c->getOrigin()) ) return false;
    if ( !this->x.sameAs(c->getX()) ) return false;
    if ( !this->y.sameAs(c->getY()) ) return false;
    if ( !this->z.sameAs(c->getZ()) ) return false;
    return true;
}

#ifdef RENDER
geom::Render_sp CoordinateSystem_O::rendered(core::List_sp koptions)
{
    geom::GrLines_sp lines = geom::GrLines_O::create();
    uint idxOrigin = lines->appendVertex(this->origin,geom::Color_O::white(_lisp));
    uint idxX = lines->appendVertex(this->origin+this->x,geom::Color_O::red(_lisp));
    uint idxY = lines->appendVertex(this->origin+this->y,geom::Color_O::green(_lisp));
    uint idxZ = lines->appendVertex(this->origin+this->z,geom::Color_O::blue(_lisp));
    lines->appendLine(idxOrigin,idxX);
    lines->appendLine(idxOrigin,idxY);
    lines->appendLine(idxOrigin,idxZ);
    geom::DisplayList_sp prims = geom::DisplayList_O::create();
    prims->append(lines);
    return prims;
}
#endif

string	CoordinateSystem_O::asString()
{
stringstream	ss;
ss << "Origin: " << this->origin.asString() << std::endl;
ss << "xdir: " << this->x.asString() << std::endl;
ss << "ydir: " << this->y.asString() << std::endl;
ss << "zdir: " << this->z.asString() << std::endl;
    return ss.str();
}

Vector3	CoordinateSystem_O::vectorRelativeToYou(const Vector3& pos)
{
    Vector3	newPos;
Matrix	transformToCanonical;
    transformToCanonical = this->matrixToCanonical();
    newPos = transformToCanonical.multiplyByVector3(pos);
    return newPos;
}







}; // namespace chem
