#ifndef	OVector2_H
#define	OVector2_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/vector2.h>

namespace geom {

    SMART(OVector2 );
    class OVector2_O : public core::CxxObject_O
    {
	LISP_BASE1(core::CxxObject_O);
	LISP_CLASS(geom,GeomPkg,OVector2_O,"OVector2");
DECLARE_INIT()
public:
	Vector2	_Value;
public:
//	void	archiveBase(core::ArchiveP node);
public:
	static OVector2_sp createFromVector2( const Vector2& pos);
    /*! Carry out a vector sum */
	static OVector2_sp sum(core::List_sp args);
public:
	string __repr__() const;
	string __str__() ;
CL_NAME("setAll");
CL_DEFMETHOD 	void setAll(double x, double y) { this->_Value.set(x,y); };
CL_NAME("x");
CL_DEFMETHOD 	double getX() { return this->_Value.getX(); };
CL_NAME("y");
CL_DEFMETHOD 	double getY() { return this->_Value.getY(); };
	Vector2	get() { return this->_Value; };
		/*! Calculate the angle from this-vb-vc
		 */
	double angle( const Vector2& vb, const Vector2& vc );
		/*! Calculate the distance from this-vb
		 */
	double distance( const Vector2& vb );
	double dotProduct( const Vector2& other);

	Vector2 sub(const Vector2& other);

	double magnitude();
	Vector2 normalized();

	/*! add two vectors
	 */
	core::T_sp add(core::List_sp points);
	Vector2 timesScalar(double s);

	DEFAULT_CTOR_DTOR(OVector2_O);
    };

    inline OVector2_sp OVector2_O::createFromVector2(const Vector2& pnt)
{
    GC_ALLOCATE(OVector2_O, ov);
    ov->setAll(pnt.getX(),pnt.getY());
    return ov;
}





}; // namespace geom
TRANSLATE(geom::OVector2_O);


namespace translate
{

template <>
struct	from_object<geom::Vector2>
{
    typedef	geom::Vector2		ExpectedType;
    typedef	geom::Vector2		DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
	geom::OVector2_sp v = o.asOrNull<geom::OVector2_O>();
	this->_v.set(v->getX(),v->getY());
    }
};


template <>
struct	from_object<const geom::Vector2&>
{
    typedef	geom::Vector2		ExpectedType;
    typedef	geom::Vector2		DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
	geom::OVector2_sp v = o.asOrNull<geom::OVector2_O>();
	this->_v.set(v->getX(),v->getY());
    }
};


template <>
struct	to_object<geom::Vector2>
{
    typedef	geom::OVector2_sp ExpectedType;
    typedef	geom::OVector2_sp DeclareType;
    DeclareType _v;
    static core::T_sp convert(geom::Vector2 pos)
    {_G();
	geom::OVector2_sp ov = geom::OVector2_O::create();
	ov->setAll(pos.getX(),pos.getY());
	return (ov);
    }
};


template <>
struct	to_object<const geom::Vector2& >
{
    typedef	geom::OVector2_sp ExpectedType;
    typedef	geom::OVector2_sp DeclareType;
    static core::T_sp convert(geom::Vector2 pos)
    {_G();
	geom::OVector2_sp ov = geom::OVector2_O::create();
	ov->setAll(pos.getX(),pos.getY());
	return (ov);
    }
};

};
#endif


