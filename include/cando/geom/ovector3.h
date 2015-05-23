#ifndef	OVECTOR3_H
#define	OVECTOR3_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lispVector.h>
#include <clasp/core/numbers.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/vector3.h>

namespace geom {

    SMART(OVector3 );
    class OVector3_O : public core::T_O
    {
	friend class OMatrix_O;
	LISP_BASE1(core::T_O);
	LISP_CLASS(geom,GeomPkg,OVector3_O,"OVector3");
    public:
        core::Cons_sp encode() const;
        void decode(core::Cons_sp);
public:
	Vector3	_Value;
    public:
	static OVector3_sp make(double x, double y, double z);
public:
//	void	archiveBase(core::ArchiveP node);
//	void	serialize(serialize::SNode node);
public:
	static OVector3_sp createFromVector3(const Vector3& pos);
	static OVector3_sp create(Vector3 const& pos);
    /*! Create from 3 element Cons */
	static OVector3_sp create(core::Cons_sp pos);
    /*! Carry out a vector sum */
	static OVector3_sp sum(core::List_sp args);
public:
    /*! Return a deep copy of this value */
	core::T_sp deepCopy() const;
	core::T_sp shallowCopy() const { return this->deepCopy();};
    string __repr__() const;
	string __str__() ;
	void setAll3(double x, double y, double z) { this->_Value.set(x,y,z); };
	double getX() { return this->_Value.getX(); };
	double getY() { return this->_Value.getY(); };
	double getZ() { return this->_Value.getZ(); };
	Vector3& value() { return this->_Value;};
	Vector3	get() { return this->_Value; };
		/*! Calculate the dihedral angle from this-vb-vc-vd
		 */
	double dihedral( const Vector3& vb, const Vector3& vc, const Vector3& vd );
		/*! Calculate the angle from this-vb-vc
		 */
	double angle( const Vector3& vb, const Vector3& vc );
		/*! Calculate the distance from this-vb
		 */
	double distance( const Vector3& vb );
	double dotProduct( const Vector3& other);
	Vector3 crossProduct(const Vector3& other);

	Vector3 sub(const Vector3& other);

	double magnitude();
	Vector3 normalized();

	/*! add two vectors
	 */
	core::T_sp add(core::List_sp points);
	Vector3 timesScalar(double s);

	DEFAULT_CTOR_DTOR(OVector3_O);
    };



}; // namespace geom
TRANSLATE(geom::OVector3_O);




namespace translate
{
    template <>
    struct	from_object<Vector3>
    {
	typedef	Vector3		ExpectedType;
	typedef	Vector3		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
          geom::OVector3_sp v = gc::As<geom::OVector3_sp>(o);
	    _v.set(v->getX(),v->getY(),v->getZ());
	}
    };
    template <>
    struct	from_object<const Vector3&>
    {
	typedef	Vector3		ExpectedType;
	typedef	Vector3		DeclareType;
	DeclareType _v;
	from_object(core::T_sp o)
	{
	    if ( geom::OVector3_sp ov3 = o.asOrNull<geom::OVector3_O>() ) {
              _v.set(ov3->getX(),ov3->getY(),ov3->getZ());
	    } else if ( core::Vector_sp vec = o.asOrNull<core::Vector_O>() ) {
              if ( core::cl_length(vec) != 3 ) {
		    SIMPLE_ERROR(BF("ovector3 can only have three elements"));
		}
	        _v.set(clasp_to_double(gc::As<core::Number_sp>((*vec)[0])),
		       clasp_to_double(gc::As<core::Number_sp>((*vec)[1])),
		       clasp_to_double(gc::As<core::Number_sp>((*vec)[2])));
	    } else if ( core::List_sp list = o.asOrNull<core::List_V>() ) {
              if ( core::cl_length(list) != 3 ) {
		    SIMPLE_ERROR(BF("ovector3 can only have three elements"));
		}
		core::T_sp e1 = oCar(list);
		list = oCdr(list);
		core::T_sp e2 = oCar(list);
		list = oCdr(list);
		core::T_sp e3 = oCar(list);
	        _v.set(core::clasp_to_double(e1.as<core::Number_O>()),
		       core::clasp_to_double(e2.as<core::Number_O>()),
		       core::clasp_to_double(e3.as<core::Number_O>()));
	    } else {
		SIMPLE_ERROR(BF("Illegal type for ovector3"));
	    }
	}
    };


    template <>
    struct	to_object<Vector3>
    {
	typedef	geom::OVector3_sp ExpectedType;
	typedef	geom::OVector3_sp DeclareType;
	static core::T_sp convert(Vector3 pos)
	{_G();
	    geom::OVector3_sp ov = geom::OVector3_O::create();
	    ov->setAll3(pos.getX(),pos.getY(),pos.getZ());
	    return (ov);
	}
    };




    template <>
    struct	to_object<const Vector3& >
    {
	typedef	geom::OVector3_sp ExpectedType;
	typedef	geom::OVector3_sp DeclareType;
	static core::T_sp convert(Vector3 pos)
	{_G();
	    geom::OVector3_sp ov = geom::OVector3_O::create();
	    ov->setAll3(pos.getX(),pos.getY(),pos.getZ());
	    return (ov);
	}
    };


};




#endif


