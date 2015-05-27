#define	DEBUG_LEVEL_FULL
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
    T_sp OVector2_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::Lisp_sp lisp)
{_G();
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
string	OVector2_O::__str__()
{
    return this->__repr__();
}

    Vector2 OVector2_O::sub(const Vector2& other)
    {
	Vector2 s = this->_Value.sub(other);
	return s;
    }

    double OVector2_O::magnitude()
    {
	double l = this->_Value.length();
	return l;
    }


    double OVector2_O::dotProduct(const Vector2& other)
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
    core::T_sp OVector2_O::add(core::List_sp points)
    {_G();
	Vector2 result = Vector2(this->getX(),this->getY());
	LOG(BF("First add = %s") % result.asString() );
	for ( auto cur : points ) {
	    OVector2_sp o = cur->car<OVector2_O>();
	    result = result.add(o->get());
	    LOG(BF("Adding %s") % o->__repr__() );
	    LOG(BF("Intermediate result %s") % result.asString() );
	}
	return OVector2_O::createFromVector2(result);
    }



Vector2 OVector2_O::timesScalar(double d)
{
    Vector2 p = this->_Value.multiplyByScalar(d);
    return p;
}


    Vector2 OVector2_O::normalized()
    {
	return this->_Value.normalized();
    }


double	OVector2_O::angle( const Vector2& vb, const Vector2& vc)
{
    return calculateAngle(this->_Value,vb,vc);
}

double	OVector2_O::distance( const Vector2& vb )
{
    return calculateDistance(this->_Value,vb);
}



    void OVector2_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<OVector2_O>()
		.def("setAll", &OVector2_O::setAll )
		.def("x", &OVector2_O::getX )
		.def("y", &OVector2_O::getY )
		.def("vector2-angle", &OVector2_O::angle)    // need translator
		.def("vector2-distance", &OVector2_O::distance)
	    .def("vector2-magnitude",&OVector2_O::magnitude)
	    .def("normalized",&OVector2_O::normalized)
	    .def("sub",&OVector2_O::sub)
	    .def("dotProduct",&OVector2_O::dotProduct)
	    .DEF(OVector2_O,add)
	    .def("timesScalar",&OVector2_O::timesScalar)
	;
    }
    void OVector2_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,OVector2,"","",_lisp)
	    ;
#endif
    }


    EXPOSE_CLASS(geom,OVector2_O);

};

