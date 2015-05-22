#define	DEBUG_LEVEL_FULL
#include "core/common.h"
#include "core/symbolToEnumConverter.h"
#include "ovector3.h"
#include "core/multipleValues.h"
#include "vector3.h"
//#include "core/serializerNode.h"
#include "core/wrappers.h"

namespace candoBase {




    OVector3_sp OVector3_O::createFromVector3(const Vector3& pnt)
{
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll(pnt.getX(),pnt.getY(),pnt.getZ());
    return ov;
}

    OVector3_sp OVector3_O::create(Vector3 const& pnt)
{
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll(pnt.getX(),pnt.getY(),pnt.getZ());
    return ov;
}

    OVector3_sp OVector3_O::create(core::Cons_sp pnt)
    {_G();
    ASSERTF(pnt->length()==3,BF("Poorly formed Vector"));
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll(core::oCar(pnt).as<core::Number_O>()->as_double(),
	       core::oCadr(pnt).as<core::Number_O>()->as_double(),
	       core::oCaddr(pnt).as<core::Number_O>()->as_double());
    return ov;
}



#define ARGS_OVector3_O_make "(x y z)"
#define DECL_OVector3_O_make ""
#define DOCS_OVector3_O_make "make OVector3 args: x y z"
    OVector3_sp OVector3_O::make(double x, double y, double z)
    {_G();
	GC_ALLOCATE(OVector3_O,ov);
	ov->_Value[0] = x;
	ov->_Value[1] = y;
	ov->_Value[2] = z;
	return ov;
    }

#if defined(OLD_SERIALIZE)
    void OVector3_O::serialize(serialize::SNode node)
    {
	node->orderedPOD(this->_Value[0]);
	CR_HINT(node,false);
	node->orderedPOD(this->_Value[1]);
	CR_HINT(node,false);
	node->orderedPOD(this->_Value[2]);
	CR_HINT(node,false);
    }
#endif


#ifdef XML_ARCHIVE
    void	OVector3_O::archiveBase(core::ArchiveP node)
{
    node->archivePlainObjectIfDefined<Vector3>( "pos","Vector3", 
    					this->_Value.isDefined(), this->_Value);
}
#endif


string	OVector3_O::__repr__() const
{
    stringstream ss;
    ss << this->_Value.getX() << " ";
    ss << this->_Value.getY() << " ";
    ss << this->_Value.getZ();
    return ss.str();
}
string	OVector3_O::__str__()
{
    return this->__repr__();
}

    Vector3 OVector3_O::sub(const Vector3& other)
    {
	Vector3 s = this->_Value.sub(other);
	return s;
    }

    double OVector3_O::magnitude()
    {
	double l = this->_Value.length();
	return l;
    }

    Vector3 OVector3_O::crossProduct(const Vector3& other)
    {
	Vector3 s = this->_Value.crossProduct(other);
	return s;
    }

    double OVector3_O::dotProduct(const Vector3& other)
    {
	double d = this->_Value.dotProduct(other);
	return d;
    }


    OVector3_sp OVector3_O::sum(core::Cons_sp args)
    {
	Vector3 sum(0.0,0.0,0.0);
	while ( args.notnilp() )
	{
	    OVector3_sp one = args.as<OVector3_O>();
	    sum = sum.add(one->get());
	}
	return OVector3_O::createFromVector3(sum);
    }


    core::T_sp OVector3_O::deepCopy() const
    {_G();
	return OVector3_O::create(this->_Value);
    }



    
    
#define ARGS_OVector3_O_add "((self ovector3) &rest points)"
#define DECL_OVector3_O_add ""
#define DOCS_OVector3_O_add "OVector3_O_add"
    core::T_sp OVector3_O::add(core::Cons_sp points)
    {_G();
	Vector3 result = Vector3(this->getX(),this->getY(),this->getZ());
	core::Cons_sp cur = points;
	while ( cur.notnilp() )
	{
	    OVector3_sp o = cur->car<OVector3_O>();
	    result = result.add(o->get());
	    LOG(BF("Adding %s") % o->__repr__() );
	    LOG(BF("Intermediate result %s") % result.asString() );
	    cur = cur->cdr();
	}
	return OVector3_O::createFromVector3(result);
    }



Vector3 OVector3_O::timesScalar(double d)
{
    Vector3 p = this->_Value.multiplyByScalar(d);
    return p;
}


    Vector3 OVector3_O::normalized()
    {
	return this->_Value.normalized(_lisp);
    }

double	OVector3_O::dihedral( const Vector3& vb, const Vector3& vc, const Vector3& vd )
{
    return calculateDihedral(this->_Value,vb,vc,vd,_lisp);
}

double	OVector3_O::angle( const Vector3& vb, const Vector3& vc)
{
    return calculateAngle(this->_Value,vb,vc,_lisp);
}

double	OVector3_O::distance( const Vector3& vb )
{
    return calculateDistance(this->_Value,vb,_lisp);
}



    void OVector3_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<OVector3_O>()
		.def("setAll", &OVector3_O::setAll )
		.def("x", &OVector3_O::getX )
		.def("y", &OVector3_O::getY )
		.def("z", &OVector3_O::getZ )
		.def("dihedral", &OVector3_O::dihedral)
		.def("vector3-angle", &OVector3_O::angle)
		.def("vector3-distance", &OVector3_O::distance)
	    .def("vector-magnitude",&OVector3_O::magnitude)
	    .def("vector-normalized",&OVector3_O::normalized)
	    .def("sub",&OVector3_O::sub)
	    .def("dotProduct",&OVector3_O::dotProduct)
	    .def("crossProduct",&OVector3_O::crossProduct)
	    .DEF(OVector3_O,add)
	    .def("timesScalar",&OVector3_O::timesScalar)
	;
	Defun_maker(CandoBasePkg,OVector3);
    }
    void OVector3_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,OVector3,"","",_lisp)
		.def("setAll", &OVector3_O::setAll )
		.def("x", &OVector3_O::getX )
		.def("y", &OVector3_O::getY )
		.def("z", &OVector3_O::getZ )
		.def("dihedral", &OVector3_O::dihedral)
		.def("angle", &OVector3_O::angle)
		.def("distance", &OVector3_O::distance)
	    .def("length",&OVector3_O::length)
	    .def("normalized",&OVector3_O::normalized)
	    .def("sub",&OVector3_O::sub)
	    .def("dotProduct",&OVector3_O::dotProduct)
	    .def("crossProduct",&OVector3_O::crossProduct)
	    .def("add",&OVector3_O::add)
	    .def("timesScalar",&OVector3_O::timesScalar)
	    ;
#endif
    }


    EXPOSE_CLASS(candoBase,OVector3_O);

};

