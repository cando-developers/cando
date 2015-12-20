#define	DEBUG_LEVEL_FULL
#include <clasp/core/common.h>
#include <clasp/core/symbolToEnumConverter.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/multipleValues.h>
#include <cando/geom/vector3.h>
//#include "core/serializerNode.h"
#include <cando/main/symbolTable.h>
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
    {_G();
    ASSERTF(pnt->length()==3,BF("Poorly formed Vector"));
    GC_ALLOCATE(OVector3_O,ov);
    ov->setAll3(core::clasp_to_double(core::oCar(pnt).as<core::Number_O>()),
	       core::clasp_to_double(core::oCadr(pnt).as<core::Number_O>()),
	       core::clasp_to_double(core::oCaddr(pnt).as<core::Number_O>()));
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

void OVector3_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,x),this->_Value[0]);
  node->field(INTERN_(kw,y),this->_Value[1]);
  node->field(INTERN_(kw,z),this->_Value[2]);
}


#ifdef XML_ARCHIVE
    void	OVector3_O::archiveBase(core::ArchiveP node)
{
    node->archivePlainObjectIfDefined<Vector3>( "pos","Vector3", 
    					this->_Value.isDefined(), this->_Value);
}
#endif

string OVector3_O::__repr__() const {
  stringstream ss;
  ss << "[ ";
  ss << this->_Value.getX() << " ";
  ss << this->_Value.getY() << " ";
  ss << this->_Value.getZ();
  ss << " ]";
  return ss.str();
}
string OVector3_O::__str__() {
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
    {_G();
	return OVector3_O::create(this->_Value);
    }



    
    
#define ARGS_OVector3_O_add "((self ovector3) &rest points)"
#define DECL_OVector3_O_add ""
#define DOCS_OVector3_O_add "OVector3_O_add"
core::T_sp OVector3_O::add(core::List_sp points)
{_G();
  Vector3 result = Vector3(this->getX(),this->getY(),this->getZ());
  for ( auto cur : points ) {
    OVector3_sp o = oCar(cur).as<OVector3_O>();
    result = result.add(o->get());
    LOG(BF("Adding %s") % o->__repr__() );
    LOG(BF("Intermediate result %s") % result.asString() );
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
  return this->_Value.normalized();
}

double	OVector3_O::dihedral( const Vector3& vb, const Vector3& vc, const Vector3& vd )
{
    return calculateDihedral(this->_Value,vb,vc,vd);
}

double	OVector3_O::angle( const Vector3& vb, const Vector3& vc)
{
    return calculateAngle(this->_Value,vb,vc);
}

double	OVector3_O::distance( const Vector3& vb )
{
    return calculateDistance(this->_Value,vb);
}

void OVector3_O::setUsingBond(double distance, OVector3_sp v)
{
  this->_Value = buildUsingBond(distance,v->_Value);
}

void OVector3_O::setUsingBondAngle(double bond, OVector3_sp bondPos,
                                     double angle, OVector3_sp anglePos)
{
  this->_Value = buildUsingBondAngle(bond,bondPos->_Value,angle,anglePos->_Value);
}

void OVector3_O::setUsingBondAngleDihedral(double bond, OVector3_sp bondPos,
                                             double angle, OVector3_sp anglePos,
                                             double dihedral, OVector3_sp dihedralPos)
{
  this->_Value = buildUsingBondAngleDihedral(bond,bondPos->_Value,
                                                      angle,anglePos->_Value,
                                                      dihedral,dihedralPos->_Value
                                                      );
}

core::List_sp OVector3_O::encode() {
  core::Vector_sp v = core::core__make_vector(cl::_sym_DoubleFloat_O,3);
  (*v)[0] = core::clasp_make_double_float(this->_Value[0]);
  (*v)[1] = core::clasp_make_double_float(this->_Value[1]);
  (*v)[2] = core::clasp_make_double_float(this->_Value[2]);
  printf("%s:%d encoded: %s\n", __FILE__, __LINE__, _rep_(v).c_str());
  return core::Cons_O::create(_Nil<T_O>(),v);
}

void OVector3_O::decode(core::List_sp r) {
  core::Vector_sp v = gc::As<core::Vector_sp>(oCdr(r));
  this->_Value[0] = core::clasp_to_double((*v)[0]);
  this->_Value[1] = core::clasp_to_double((*v)[1]);
  this->_Value[2] = core::clasp_to_double((*v)[2]);
}


void OVector3_O::exposeCando(core::Lisp_sp lisp)
{
  core::class_<OVector3_O>()
    .def("setAll3", &OVector3_O::setAll3 )
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
    .def("add",&OVector3_O::add,
         ARGS_OVector3_O_add,
         DECL_OVector3_O_add,
         DOCS_OVector3_O_add)
    .def("timesScalar",&OVector3_O::timesScalar)
    .def("setUsingBond",&OVector3_O::setUsingBond)
    .def("setUsingBondAngle",&OVector3_O::setUsingBondAngle)
    .def("setUsingBondAngleDihedral",&OVector3_O::setUsingBondAngleDihedral)
    ;
  Defun_maker(GeomPkg,OVector3);
}
    void OVector3_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,OVector3,"","",_lisp)
		.def("setAll3", &OVector3_O::setAll3 )
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


    EXPOSE_CLASS(geom,OVector3_O);

};

