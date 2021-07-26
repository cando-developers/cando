/*
    File: coordinateArray.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *      File:   coordinateArray.h
 *
 *      Description:
 *              Routines to maintain a vector of Vector3's
 */

#ifndef SimpleVectorCoordinate_H
#define SimpleVectorCoordinate_H

#include <vector>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/array.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/vector3.h>
#include <cando/geom/ovector3.h>
#include <cando/geom/matrix.h>
#include <clasp/core/lispStream.fwd.h>


namespace geom {
  extern core::Symbol_sp& _sym_OVector3_O;
};

// ----------------------------------------------------------------------
// Arrays specialized for Coordinate
//

namespace geom {
  FORWARD(SimpleVectorCoordinate);
  FORWARD(ComplexVectorCoordinate);
};
template <>
struct gctools::GCInfo<geom::SimpleVectorCoordinate_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = atomic;
};
namespace geom {
  class SimpleVectorCoordinate_O;
  typedef core::template_SimpleVector<SimpleVectorCoordinate_O,Vector3,core::AbstractSimpleVector_O> specialized_SimpleVectorCoordinate;
  class SimpleVectorCoordinate_O : public specialized_SimpleVectorCoordinate {
    LISP_CLASS(geom, GeomPkg, SimpleVectorCoordinate_O, "SimpleVectorCoordinate",core::AbstractSimpleVector_O);
    virtual ~SimpleVectorCoordinate_O() {};
  public:
    typedef specialized_SimpleVectorCoordinate TemplatedBase;
    typedef typename TemplatedBase::leaf_type leaf_type;
    typedef typename TemplatedBase::value_type value_type;
    typedef typename TemplatedBase::simple_element_type simple_element_type;
    typedef typename TemplatedBase::vector_type vector_type;
    typedef typename TemplatedBase::iterator iterator;
    typedef typename TemplatedBase::const_iterator const_iterator;
    typedef value_type container_value_type;
  public:
    static value_type default_initial_element(void) { return Vector3(0.0,0.0,0.0); }
    static value_type from_object(core::T_sp obj) { return translate::from_object<Vector3>(obj)._v; };
    static core::T_sp to_object(const value_type& v) { return OVector3_O::create(v);};
  public:
  SimpleVectorCoordinate_O(size_t length, value_type initialElement=value_type(),
                           bool initialElementSupplied=false,
                           size_t initialContentsSize=0,
                           const value_type* initialContents=NULL)
    : TemplatedBase(length,initialElement,initialElementSupplied,initialContentsSize,initialContents) {};
    static SimpleVectorCoordinate_sp make(size_t length,
                                          value_type initialElement=value_type(),
                                          bool initialElementSupplied=false,
                                          size_t initialContentsSize=0,
                                          const value_type* initialContents=NULL,
                                          bool static_vector_p = false) {
      auto bs = gctools::GC<SimpleVectorCoordinate_O>::allocate_container(static_vector_p,length,initialElement,initialElementSupplied,initialContentsSize,initialContents);
      return bs;
    }
    smart_ptr_type copy(size_t length, value_type initialElement, bool initialElementSupplied) {
      return make(length, initialElement, initialElementSupplied,
                  MIN(length, this->length()), this->begin());
    }
    smart_ptr_type copy() {
      return copy(this->length(), default_initial_element(), false);
    }
    static SimpleVectorCoordinate_sp create(core::List_sp elements) {
      size_t sz = core::cl__length(elements);
      size_t i(0);
      SimpleVectorCoordinate_sp v = make(sz,Vector3(),false,0,NULL);
      for ( auto cur : elements ) {
        (*v)[i++] = gc::As<OVector3_sp>(core::oCar(cur))->get();
      }
      return v;
    }
    
  public:
    // Specific to SimpleVectorCoordinate_O
//    virtual void __write__(T_sp stream) const final;
  public:
    virtual core::T_sp element_type() const override { return geom::_sym_OVector3_O;};
  public:
    Vector3& getElement(size_t i) { return (*this)[i];};
    const Vector3& getElement(size_t i) const { return (*this)[i];};
    core::T_sp getElementObject(size_t i) const { return this->rowMajorAref(i); };
    core::T_sp setElementObject(size_t i, core::T_sp obj) {this->rowMajorAset(i,obj); return obj;};
    void zeroAllCoordinates();
    void addPoints(SimpleVectorCoordinate_sp y);
    double distanceTo(SimpleVectorCoordinate_sp y);
  };
};

namespace geom {
  FORWARD(SimpleMDArrayCoordinate);
  class SimpleMDArrayCoordinate_O : public core::template_SimpleArray<SimpleMDArrayCoordinate_O,SimpleVectorCoordinate_O,core::SimpleMDArray_O> {
    LISP_CLASS(geom, GeomPkg, SimpleMDArrayCoordinate_O, "SimpleMDArrayCoordinate",core::SimpleMDArray_O);
    virtual ~SimpleMDArrayCoordinate_O() {};
  public:
    typedef core::template_SimpleArray<SimpleMDArrayCoordinate_O,SimpleVectorCoordinate_O,core::SimpleMDArray_O> TemplatedBase;
    typedef typename TemplatedBase::simple_element_type simple_element_type;
    typedef typename TemplatedBase::simple_type simple_type;
  public: // make array
  SimpleMDArrayCoordinate_O(size_t rank,
                  core::List_sp dimensions,
                  core::Array_sp data) : TemplatedBase(rank,dimensions,data) {};
    static SimpleMDArrayCoordinate_sp make_multi_dimensional(core::List_sp dim_desig, simple_element_type initialElement, core::T_sp data) {
      ASSERT(dim_desig.consp()||dim_desig.nilp());
      size_t rank;
      size_t arrayTotalSize = calculateArrayTotalSizeAndValidateDimensions(dim_desig,rank);
      LIKELY_if (data.nilp()) {
        data = SimpleVectorCoordinate_O::make(arrayTotalSize,initialElement,true);
      }
      SimpleMDArrayCoordinate_sp array = gctools::GC<SimpleMDArrayCoordinate_O>::allocate_container(false,rank,dim_desig,gc::As<core::Array_sp>(data));
      return array;
    }
  };
};



namespace geom {
  FORWARD(MDArrayCoordinate);
};

namespace geom {
  class MDArrayCoordinate_O : public core::template_Array<MDArrayCoordinate_O,SimpleMDArrayCoordinate_O,SimpleVectorCoordinate_O,core::MDArray_O> {
    LISP_CLASS(geom, GeomPkg, MDArrayCoordinate_O, "MDArrayCoordinate",core::MDArray_O);
    virtual ~MDArrayCoordinate_O() {};
  public:
    typedef core::template_Array<MDArrayCoordinate_O,SimpleMDArrayCoordinate_O,SimpleVectorCoordinate_O,core::MDArray_O> TemplatedBase;
    typedef typename TemplatedBase::simple_element_type simple_element_type;
    typedef typename TemplatedBase::simple_type simple_type;
  public: // make array
  MDArrayCoordinate_O(size_t rank,
                  core::List_sp dimensions,
                  core::Array_sp data,
                  bool displacedToP,
                  core::Fixnum_sp displacedIndexOffset) : TemplatedBase(rank,dimensions,data,displacedToP,displacedIndexOffset) {};
    static MDArrayCoordinate_sp make_multi_dimensional(core::List_sp dim_desig, simple_element_type initialElement, core::T_sp dataOrDisplacedTo, bool displacedToP, core::Fixnum_sp displacedIndexOffset) {
      ASSERT(dim_desig.consp()||dim_desig.nilp());
      size_t rank;
      size_t arrayTotalSize = calculateArrayTotalSizeAndValidateDimensions(dim_desig,rank);
      LIKELY_if (dataOrDisplacedTo.nilp()) {
        dataOrDisplacedTo = simple_type::make(arrayTotalSize,initialElement,true);
      }
      MDArrayCoordinate_sp array = gctools::GC<MDArrayCoordinate_O>::allocate_container(false,rank,dim_desig,gc::As<core::Array_sp>(dataOrDisplacedTo),displacedToP,displacedIndexOffset);
      return array;
    }
    size_t vectorPush_Vector3(const Vector3& newElement);
    size_t vectorPushExtend_Vector3(const Vector3& newElement, size_t extension = 0);
  public:
//    virtual bool equalp(T_sp o) const final;
  };
};

namespace geom {
  class ComplexVectorCoordinate_O : public core::template_Vector<ComplexVectorCoordinate_O,SimpleVectorCoordinate_O,core::ComplexVector_O> {
    LISP_CLASS(geom, GeomPkg, ComplexVectorCoordinate_O, "ComplexVectorCoordinate",core::ComplexVector_O);
    virtual ~ComplexVectorCoordinate_O() {};
  public:
    typedef core::template_Vector<ComplexVectorCoordinate_O,SimpleVectorCoordinate_O,core::ComplexVector_O> TemplatedBase;
    typedef typename TemplatedBase::simple_element_type simple_element_type;
    typedef typename TemplatedBase::simple_type simple_type;
  public: // make vector
    ComplexVectorCoordinate_O(size_t rank1,
                              size_t dimension,
                              core::T_sp fillPointer,
                              core::Array_sp data,
                              bool displacedToP,
                              core::Fixnum_sp displacedIndexOffset) : TemplatedBase(dimension,fillPointer,data,displacedToP,displacedIndexOffset) {};
    static ComplexVectorCoordinate_sp make_vector(size_t dimension, simple_element_type initialElement/*=simple_element_type()*/, core::T_sp fillPointer/*=_Nil<T_O>()*/, core::T_sp dataOrDisplacedTo/*=_Nil<T_O>()*/, bool displacedToP/*=false*/, core::Fixnum_sp displacedIndexOffset/*=clasp_make_fixnum(0)*/ ) {
      LIKELY_if (dataOrDisplacedTo.nilp()) {
        dataOrDisplacedTo = simple_type::make(dimension,initialElement,true);
      }
      ComplexVectorCoordinate_sp array = gctools::GC<ComplexVectorCoordinate_O>::allocate_container(false,1/*CRANK*/,dimension,fillPointer,gc::As_unsafe<core::Array_sp>(dataOrDisplacedTo),displacedToP,displacedIndexOffset);
      return array;
    }
    static ComplexVectorCoordinate_sp make_vector(size_t dimension, simple_element_type initialElement, core::T_sp fillPointer) {
      return make_vector(dimension,initialElement,fillPointer,nil<T_O>(),false,core::clasp_make_fixnum(0));
    }
    /*! Build a ComplexVectorCoordinate from a SimpleVectorCoordinate */
    static ComplexVectorCoordinate_sp make_vector(geom::SimpleVectorCoordinate_sp data) {
      return make_vector(data->length(),Vector3(),core::make_fixnum(data->length()),data,false,core::make_fixnum(0));
    }
    size_t vectorPush_Vector3(const Vector3& newElement);
    size_t vectorPushExtend_Vector3(const Vector3& newElement, size_t extension = 0);
  public:
//    virtual bool equalp(T_sp o) const final;
  };
};

namespace geom {

  typedef core::Array_O  ArrayCoordinate_O;
  typedef gc::smart_ptr<ArrayCoordinate_O> ArrayCoordinate_sp;
  
  void geom__in_place_transform(ArrayCoordinate_sp array, const Matrix& transform);
  void geom__in_place__TIMES_scalar(ArrayCoordinate_sp array, double mul);


};
#endif          /* ifndef SimpleVectorCoordinate_H */
