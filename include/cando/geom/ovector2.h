/*
    File: ovector2.h
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
	LISP_CLASS(geom,GeomPkg,OVector2_O,"OVector2",core::CxxObject_O);
public:
	Vector2	_Value;
public:
public:
	static OVector2_sp createFromVector2( const Vector2& pos);
    /*! Carry out a vector sum */
	static OVector2_sp sum(core::List_sp args);
public:
	string __repr__() const;
CL_LISPIFY_NAME("setAll");
CL_DEFMETHOD 	void setAll(double x, double y) { this->_Value.set(x,y); };
CL_LISPIFY_NAME("x");
CL_DEFMETHOD 	double getX() { return this->_Value.getX(); };
CL_LISPIFY_NAME("y");
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
  auto  ov = gctools::GC<OVector2_O>::allocate();
    ov->setAll(pnt.getX(),pnt.getY());
    return ov;
}





}; // namespace geom


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
