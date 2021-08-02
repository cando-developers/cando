/*
    File: ovector2.cc
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
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/multipleValues.h>
#include <cando/geom/ovector2.h>
#include <cando/geom/vector2.h>
//#include "core/archiveNode.h"
#include <clasp/core/environment.h>
#include <clasp/core/wrappers.h>

namespace geom {

#if 0
    T_sp OVector2_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::LispPtr lisp)
{
    double x = translate::from_object<double>::convert(bargs->lookup(CorePkg,"tx")->object());
    double y = translate::from_object<double>::convert(bargs->lookup(CorePkg,"ty")->object());
    LOG(BF("Got values for vector x,y,z = %lf, %lf") % x % y );
    this->_Value[0] = x;
    this->_Value[1] = y;
    return _Nil<core::T_O>();
}
#endif

#ifdef XML_ARCHIVE
    void	OVector2_O::archiveBase(core::ArchiveP node)
{
    node->archivePlainObjectIfDefined<Vector2>( "pos","Vector2", 
    					this->_Value.isDefined(), this->_Value);
}
#endif

string	OVector2_O::__repr__() const
{
    stringstream ss;
    ss << this->_Value.getX() << " ";
    ss << this->_Value.getY();
    return ss.str();
}

CL_LISPIFY_NAME("sub");
CL_DEFMETHOD     Vector2 OVector2_O::sub(const Vector2& other)
    {
	Vector2 s = this->_Value.sub(other);
	return s;
    }

CL_LISPIFY_NAME("vector2-magnitude");
CL_DEFMETHOD     double OVector2_O::magnitude()
    {
	double l = this->_Value.length();
	return l;
    }


CL_LISPIFY_NAME("dotProduct");
CL_DEFMETHOD     double OVector2_O::dotProduct(const Vector2& other)
    {
	double d = this->_Value.dotProduct(other);
	return d;
    }


    OVector2_sp OVector2_O::sum(core::List_sp args)
    {
      Vector2 sum(0.0,0.0);
      while ( args.notnilp() )
      {
        OVector2_sp one = gc::As<OVector2_sp>(oCar(args));
        sum = sum.add(one->get());
      }
      return OVector2_O::createFromVector2(sum);
    }



#define ARGS_OVector2_O_add "((self ovector2) &rest points)"
#define DECL_OVector2_O_add ""
#define DOCS_OVector2_O_add "Add a collection of points together"
CL_LISPIFY_NAME("add");
CL_DEFMETHOD     core::T_sp OVector2_O::add(core::List_sp points)
    {
	Vector2 result = Vector2(this->getX(),this->getY());
	LOG(BF("First add = %s") % result.asString() );
	for ( auto cur : points ) {
	    OVector2_sp o = cur->car<OVector2_O>();
	    result = result.add(o->get());
	    LOG(BF("Adding %s") % _rep_(o) );
	    LOG(BF("Intermediate result %s") % result.asString() );
	}
	return OVector2_O::createFromVector2(result);
    }



CL_LISPIFY_NAME("timesScalar");
CL_DEFMETHOD Vector2 OVector2_O::timesScalar(double d)
{
    Vector2 p = this->_Value.multiplyByScalar(d);
    return p;
}


CL_LISPIFY_NAME("normalized");
CL_DEFMETHOD     Vector2 OVector2_O::normalized()
    {
	return this->_Value.normalized();
    }


CL_LISPIFY_NAME("vector2-angle");
CL_DEFMETHOD double	OVector2_O::angle( const Vector2& vb, const Vector2& vc)
{
  return geom::calculate2Angle(this->_Value,vb,vc);
}

CL_LISPIFY_NAME("vector2-distance");
CL_DEFMETHOD double	OVector2_O::distance( const Vector2& vb )
{
    return calculate2Distance(this->_Value,vb);
}






};
