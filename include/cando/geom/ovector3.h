/*
    File: ovector3.h
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
#ifndef	OVECTOR3_H
#define	OVECTOR3_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/numbers.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/vector3.h>

namespace geom {

    SMART(OVector3 );
    class OVector3_O : public core::CxxObject_O
    {
	friend class OMatrix_O;
	LISP_CLASS(geom,GeomPkg,OVector3_O,"V3",core::CxxObject_O);
    public:
        bool fieldsp() const { return true; };
        void	fields(core::Record_sp node);
    public:
#if 0
        core::List_sp encode() ;
        void decode(core::List_sp);
#endif
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
CL_LISPIFY_NAME("setAll3");
CL_DEFMETHOD 	void setAll3(double x, double y, double z) { this->_Value.set(x,y,z); };
CL_LISPIFY_NAME("vx");
CL_DEFMETHOD 	double getX() { return this->_Value.getX(); };
CL_LISPIFY_NAME("vy");
CL_DEFMETHOD 	double getY() { return this->_Value.getY(); };
CL_LISPIFY_NAME("vz");
CL_DEFMETHOD 	double getZ() { return this->_Value.getZ(); };
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
        
        Vector3 _PLUS_(const Vector3& other);
        Vector3 _MINUS_(const Vector3& other);
        
	double magnitude();
	Vector3 normalized();

	/*! add two vectors
	 */
	core::T_sp add(core::List_sp points);
	Vector3 timesScalar(double s);

        void setUsingBond(double distance, OVector3_sp v);
        void setUsingBondAngle(double bond, OVector3_sp bondPos,
                                 double angle, OVector3_sp anglePos);
        void setUsingBondAngleDihedral(double bond, OVector3_sp bondPos,
                                        double angle, OVector3_sp anglePos,
                                        double dihedral, OVector3_sp dihedralPos);
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
          if ( core::Vector_sp vo = o.asOrNull<core::Vector_O>() ) {
            if ( cl__length(vo) == 3 ) {
              _v.set(core::clasp_to_double(vo->rowMajorAref(0)),
                     core::clasp_to_double(vo->rowMajorAref(1)),
                     core::clasp_to_double(vo->rowMajorAref(2)));
              return;
            } else {
              SIMPLE_ERROR(BF("Vector must have 3 elements"));
            }
          } else if (geom::OVector3_sp ov = o.asOrNull<geom::OVector3_O>()) {
	    _v.set(ov->getX(),ov->getY(),ov->getZ());
            return;
          }
          SIMPLE_ERROR(BF("You must provide an OVector3 or a Vector of three elements"));
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
              if ( core::cl__length(vec) != 3 ) {
		    SIMPLE_ERROR(BF("ovector3 can only have three elements"));
		}
              _v.set(clasp_to_double(gc::As<core::Number_sp>(vec->rowMajorAref(0))),
                     clasp_to_double(gc::As<core::Number_sp>(vec->rowMajorAref(1))),
                     clasp_to_double(gc::As<core::Number_sp>(vec->rowMajorAref(2))));
	    } else if ( core::List_sp list = o.asOrNull<core::List_V>() ) {
              if ( core::cl__length(list) != 3 ) {
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
