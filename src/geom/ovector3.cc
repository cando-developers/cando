/*
    File: ovector3.cc
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
#include <clasp/core/common.h>
#include <clasp/core/symbolToEnumConverter.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/multipleValues.h>
#include <cando/geom/vector3.h>
//#include "core/serializerNode.h"
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>

namespace geom {




    OVector3_sp OVector3_O::createFromVector3(const Vector3& pnt)
{
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll3(pnt.getX(),pnt.getY(),pnt.getZ());
    return ov;
}

    OVector3_sp OVector3_O::create(Vector3 const& pnt)
{
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll3(pnt.getX(),pnt.getY(),pnt.getZ());
    return ov;
}

    OVector3_sp OVector3_O::create(core::Cons_sp pnt)
    {
    ASSERTF(pnt->length()==3,BF("Poorly formed Vector"));
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll3(core::clasp_to_double(core::oCar(pnt).as<core::Number_O>()),
	       core::clasp_to_double(core::oCadr(pnt).as<core::Number_O>()),
	       core::clasp_to_double(core::oCaddr(pnt).as<core::Number_O>()));
    return ov;
}



CL_NAME("VEC");
CL_DEFUN OVector3_sp OVector3_O::make(double x, double y, double z)
    {
	GC_ALLOCATE(OVector3_O,ov);
	ov->_Value[0] = x;
	ov->_Value[1] = y;
	ov->_Value[2] = z;
	return ov;
    }



#ifdef XML_ARCHIVE
    void	OVector3_O::archiveBase(core::ArchiveP node)
{
    node->archivePlainObjectIfDefined<Vector3>( "pos","V3", 
    					this->_Value.isDefined(), this->_Value);
}
#endif

string OVector3_O::__repr__() const {
  stringstream ss;
  ss << "'( ";
  ss << this->_Value.getX() << " ";
  ss << this->_Value.getY() << " ";
  ss << this->_Value.getZ();
  ss << ")";
  return ss.str();
}

CL_NAME("V-");
CL_DEFMETHOD     Vector3 OVector3_O::sub(const Vector3& other)
    {
	Vector3 s = this->_Value.sub(other);
	return s;
    }

CL_NAME("V+");
CL_DEFMETHOD Vector3 OVector3_O::_PLUS_(const Vector3& other)
{
  Vector3 s = this->_Value.add(other);
  return s;
}

CL_NAME("VLENGTH");
CL_DEFMETHOD     double OVector3_O::magnitude()
    {
	double l = this->_Value.length();
	return l;
    }

CL_LISPIFY_NAME("VCROSS");
CL_DEFMETHOD     Vector3 OVector3_O::crossProduct(const Vector3& other)
    {
	Vector3 s = this->_Value.crossProduct(other);
	return s;
    }

CL_NAME("VDOT");
CL_DEFMETHOD     double OVector3_O::dotProduct(const Vector3& other)
    {
	double d = this->_Value.dotProduct(other);
	return d;
    }


OVector3_sp OVector3_O::sum(core::List_sp args)
{
  Vector3 sum(0.0,0.0,0.0);
  while ( args.notnilp() )
  {
    OVector3_sp one = gc::As<OVector3_sp>(args);
    sum = sum.add(one->get());
  }
  return OVector3_O::createFromVector3(sum);
}


    core::T_sp OVector3_O::deepCopy() const
    {
	return OVector3_O::create(this->_Value);
    }



    
    
CL_NAME("VSUM");
CL_DEFMETHOD core::T_sp OVector3_O::add(core::List_sp points)
{
  Vector3 result = Vector3(this->getX(),this->getY(),this->getZ());
  for ( auto cur : points ) {
    OVector3_sp o = oCar(cur).as<OVector3_O>();
    result = result.add(o->get());
    LOG(BF("Adding %s") % _rep_(o) );
    LOG(BF("Intermediate result %s") % result.asString() );
  }
  return OVector3_O::createFromVector3(result);
}



CL_NAME("V*");
CL_DEFMETHOD Vector3 OVector3_O::timesScalar(double d)
{
    Vector3 p = this->_Value.multiplyByScalar(d);
    return p;
}


CL_NAME("VNORMALIZED");
CL_DEFMETHOD Vector3 OVector3_O::normalized()
{
  return this->_Value.normalized();
}

CL_NAME("VDIHEDRAL");
CL_DEFMETHOD double	OVector3_O::dihedral( const Vector3& vb, const Vector3& vc, const Vector3& vd )
{
  return geom::calculateDihedral(this->_Value,vb,vc,vd);
}

CL_NAME("VANGLE");
CL_DEFMETHOD double	OVector3_O::angle( const Vector3& vb, const Vector3& vc)
{
  return geom::calculateAngle(this->_Value,vb,vc);
}

CL_NAME("VDISTANCE");
CL_DEFMETHOD double	OVector3_O::distance( const Vector3& vb )
{
  return geom::calculateDistance(this->_Value,vb);
}

CL_LISPIFY_NAME("v_setUsingBond");
CL_DEFMETHOD void OVector3_O::setUsingBond(double distance, OVector3_sp v)
{
  this->_Value = buildUsingBond(distance,v->_Value);
}

CL_LISPIFY_NAME("v_setUsingBondAngle");
CL_DEFMETHOD void OVector3_O::setUsingBondAngle(double bond, OVector3_sp bondPos,
                                     double angle, OVector3_sp anglePos)
{
  this->_Value = buildUsingBondAngle(bond,bondPos->_Value,angle,anglePos->_Value);
}

CL_LISPIFY_NAME("v_setUsingBondAngleDihedral");
CL_DEFMETHOD void OVector3_O::setUsingBondAngleDihedral(double bond, OVector3_sp bondPos,
                                             double angle, OVector3_sp anglePos,
                                             double dihedral, OVector3_sp dihedralPos)
{
  this->_Value = buildUsingBondAngleDihedral(bond,bondPos->_Value,
                                                      angle,anglePos->_Value,
                                                      dihedral,dihedralPos->_Value
                                                      );
}

void OVector3_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,x),this->_Value[0]);
  node->field(INTERN_(kw,y),this->_Value[1]);
  node->field(INTERN_(kw,z),this->_Value[2]);
}



};
