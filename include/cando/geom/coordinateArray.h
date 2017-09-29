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



#if 0
namespace geom {
  FORWARD(SimpleVectorCoordinate);
  extern core::Symbol_sp& arr_sym_vector3;
};

template <>
struct gctools::GCInfo<geom::SimpleVectorCoordinate_O> {
  static bool constexpr NeedsInitialization = false;
  static bool constexpr NeedsFinalization = false;
  static GCInfo_policy constexpr Policy = atomic;
};


namespace geom
{
  SMART(DisplayList);
  SMART(SimpleVectorCoordinate);
  typedef core::abstract_SimpleVector<SimpleVectorCoordinate_O,Vector3,core::AbstractSimpleVector_O> abstract_SimpleVectorCoordinate;
  class SimpleVectorCoordinate_O : public abstract_SimpleVectorCoordinate
  {
    LISP_CLASS(geom,GeomPkg,SimpleVectorCoordinate_O,"SimpleVectorCoordinate",core::AbstractSimpleVector_O);
  public:
#if 0
    bool fieldsp() const { return true; };
    void	fields(core::Record_sp node);
#endif
  public:
//    static SimpleVectorCoordinate_sp make(core::T_sp fnsize,core::List_sp vals);
  public:
    void initialize();
  public:
    virtual core::List_sp encode();
    virtual void decode(core::List_sp);
  private:
    void fillFromList(core::List_sp vals);
  public:
    static SimpleVectorCoordinate_sp create(uint size)
    {
      GC_ALLOCATE(SimpleVectorCoordinate_O,a);
      a->resize(size);
      return a;
    }
    static SimpleVectorCoordinate_sp create(core::List_sp vertexCons);

    typedef abstract_SimpleVectorCoordinate TemplatedBase;
    typedef typename TemplatedBase::leaf_type leaf_type;
    typedef typename TemplatedBase::value_type value_type;
    typedef typename TemplatedBase::vector_type vector_type;
    typedef typename TemplatedBase::iterator iterator;
    typedef typename TemplatedBase::const_iterator const_iterator;
    typedef value_type container_value_type;
  public:
    static value_type initial_element_from_object(core::T_sp obj, bool supplied) {
      IMPLEMENT_ME();
    }
    static value_type from_object(core::T_sp obj) {return gc::As<geom::OVector3_sp>(obj)->get(); }
    static core::T_sp to_object(const value_type& v) { return geom::OVector3_O::create(v); };
  public:
  SimpleVectorCoordinate_O(size_t length, value_type initialElement=value_type(), bool initialElementSupplied=false, size_t initialContentsSize=0, const value_type* initialContents=NULL) : TemplatedBase(length,initialElement,initialElementSupplied,initialContentsSize,initialContents) {};
    static SimpleVectorCoordinate_sp make(size_t length, value_type initialElement=value_type(), bool initialElementSupplied=false, size_t initialContentsSize=0, const value_type* initialContents=NULL) {
      auto bs = gctools::GC<SimpleVectorCoordinate_O>::allocate_container(length,length,initialElement,initialElementSupplied,initialContentsSize,initialContents);
      return bs;
    }
  public:
    // Specific to SimpleVectorCoordinate_O
//    virtual void __write__(T_sp stream) const final;
  public:
    virtual core::T_sp array_type() const final { return cl::_sym_simple_array; };
    virtual core::T_sp element_type() const override { return geom::_sym_vector3;};
  public:
    virtual core::clasp_elttype elttype() const { IMPLEMENT_MEF(BF("Implement elttype")); };
    virtual core::T_sp arrayElementType() const override { return geom::_sym_vector3; };
  public:
    // Implement these methods for simple vectors - some are implemented in parent classes
    // for convenience if not speed
    virtual bool equal(core::T_sp other) const override { return this->eq(other);};
  public:
    SimpleVectorCoordinate_O( const SimpleVectorCoordinate_O& orig );

    INHERIT_SEQUENCE virtual core::T_sp elt(int index) const { return this->getElementObject(index); };
    INHERIT_SEQUENCE virtual core::T_sp setf_elt(int index, core::T_sp value) { return this->setElementObject(index,value); };

    iterator begin() { return &this->_Data[0];};
    iterator end() { return &this->_Data[this->_Data._Length]; }
    const_iterator begin() const { return &this->_Data[0];};
    const_iterator end() const { return &this->_Data[this->_Data._Length]; }

    void zeroAllCoordinates();
    void add(SimpleVectorCoordinate_sp c);
//    void multiplyByScalar(double d);

    double distanceTo(SimpleVectorCoordinate_sp arr);

    Vector3& getElement(uint i);
    Vector3 value_getElement(uint i);
    core::T_sp getElementObject(uint i) const;
    void	setElement(uint i, Vector3 pos );
    core::T_sp setElementObject(uint i, core::T_sp val);
#if 0
    void	appendElement(const Vector3& pos) { this->_Points.push_back(pos); };
    CL_LISPIFY_NAME("coordinate-array-append");
    CL_DEFMETHOD     void	append(const Vector3& pos) { this->_Points.push_back(pos); };
    void	push_back(const Vector3& pos) { this->_Points.push_back(pos); };
#endif
    void	resize(int i);
//    void	clear() { this->_Points.clear();};
    void	transform(const Matrix& transMatrix);
    void	multiplyBy(double multiplier);

    string __repr__() const;
    string	asString() const;

#ifdef RENDER
    Render_sp asDisplayListOfPoints();
#endif

//	string parseFromStream(core::Stream_sp stream);
    void writeToStream(string const& info, core::T_sp stream);

    SimpleVectorCoordinate_sp	copy() const;
    virtual core::T_sp deepCopy() const { return this->copy();};
    virtual core::T_sp shallowCopy() const { return this->copy();};
  };
};


#endif





// ----------------------------------------------------------------------
// Arrays specialized for Coordinate
//

namespace geom {
  FORWARD(SimpleVectorCoordinate);
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
    static value_type initial_element_from_object(core::T_sp obj, bool supplied) {
      if (supplied) {
        return translate::from_object<Vector3>(obj)._v;
      }
      return Vector3(0.0,0.0,0.0);
    }
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
                                          const value_type* initialContents=NULL) {
      auto bs = gctools::GC<SimpleVectorCoordinate_O>::allocate_container(length,length,initialElement,initialElementSupplied,initialContentsSize,initialContents);
      return bs;
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
    virtual core::T_sp array_type() const final { return cl::_sym_simple_array; };
    virtual core::T_sp element_type() const override { return geom::_sym_OVector3_O;};
    virtual core::T_sp arrayElementType() const override { return geom::_sym_OVector3_O; };
    virtual core::clasp_elttype elttype() const { return core::clasp_aet_non_standard; };
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
  public: // make vector
  SimpleMDArrayCoordinate_O(size_t rank1, size_t dimension, core::Array_sp data) : TemplatedBase(dimension,data) {};
    static SimpleMDArrayCoordinate_sp make_vector(size_t dimension, simple_element_type initialElement/*=_Nil<T_O>()*/, core::T_sp data/*=_Nil<T_O>()*/) {
      LIKELY_if (data.nilp()) {
        data = SimpleVectorCoordinate_O::make(dimension,initialElement,true);
      }
      SimpleMDArrayCoordinate_sp array = gctools::GC<SimpleMDArrayCoordinate_O>::allocate_container(1,1,dimension,gc::As_unsafe<core::Array_sp>(data));
      return array;
    }
    static SimpleMDArrayCoordinate_sp make_vector(size_t dimension, simple_element_type initialElement) {
      return make_vector(dimension,initialElement,_Nil<T_O>());
    }
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
      SimpleMDArrayCoordinate_sp array = gctools::GC<SimpleMDArrayCoordinate_O>::allocate_container(rank,rank,dim_desig,gc::As<core::Array_sp>(data));
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
  public: // make vector
  MDArrayCoordinate_O(size_t dummy_rank_1,
                      size_t dimension,
                      core::T_sp fillPointer,
                      core::Array_sp data,
                      bool displacedToP,
                      core::Fixnum_sp displacedIndexOffset) : TemplatedBase(core::Rank1(),dimension,fillPointer,data,displacedToP,displacedIndexOffset) {};
    static MDArrayCoordinate_sp make_vector(size_t dimension, simple_element_type initialElement/*=simple_element_type()*/, core::T_sp fillPointer/*=_Nil<T_O>()*/, core::T_sp dataOrDisplacedTo/*=_Nil<T_O>()*/, bool displacedToP/*=false*/, core::Fixnum_sp displacedIndexOffset/*=clasp_make_fixnum(0)*/ ) {
      LIKELY_if (dataOrDisplacedTo.nilp()) {
        dataOrDisplacedTo = simple_type::make(dimension,initialElement,true);
      }
      MDArrayCoordinate_sp array = gctools::GC<MDArrayCoordinate_O>::allocate_container(1,1,dimension,fillPointer,gc::As_unsafe<core::Array_sp>(dataOrDisplacedTo),displacedToP,displacedIndexOffset);
      return array;
    }
    static MDArrayCoordinate_sp make_vector(size_t dimension, simple_element_type initialElement, core::T_sp fillPointer) {
      return make_vector(dimension,initialElement,fillPointer,_Nil<T_O>(),false,core::clasp_make_fixnum(0));
    }
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
      MDArrayCoordinate_sp array = gctools::GC<MDArrayCoordinate_O>::allocate_container(rank,rank,dim_desig,gc::As<core::Array_sp>(dataOrDisplacedTo),displacedToP,displacedIndexOffset);
      return array;
    }
    size_t vectorPush_Vector3(const Vector3& newElement);
    size_t vectorPushExtend_Vector3(const Vector3& newElement, size_t extension = 0);
  public:
//    virtual bool equalp(T_sp o) const final;
  };
};




#if 0
  SMART(SimpleVectorCoordinateWithHash);
  class SimpleVectorCoordinateWithHash_O : public SimpleVectorCoordinate_O
  {
    LISP_CLASS(geom,GeomPkg,SimpleVectorCoordinateWithHash_O,"SimpleVectorCoordinateWithHash",SimpleVectorCoordinate_O);
//	DECLARE_SERIALIZE();
  public:
    void initialize();
  public:
//	void	archiveBase(core::ArchiveP node);
  private:
    uint	_Hash;
  public:
    static SimpleVectorCoordinateWithHash_sp create(uint sz);
  public:
    SimpleVectorCoordinateWithHash_O( const SimpleVectorCoordinateWithHash_O& orig );


CL_LISPIFY_NAME("coordinate-array-with-hash-gethash");
CL_DEFMETHOD     uint getHash() { return this->_Hash;};
CL_LISPIFY_NAME("coordinate-array-with-hash-setHash");
CL_DEFMETHOD     void setHash(uint h) { this->_Hash = h;};

    DEFAULT_CTOR_DTOR(SimpleVectorCoordinateWithHash_O);
  };
#endif



namespace geom {

  typedef core::Array_O  ArrayCoordinate_O;
  typedef gc::smart_ptr<ArrayCoordinate_O> ArrayCoordinate_sp;
  
  void geom__in_place_transform(ArrayCoordinate_sp array, const Matrix& transform);
  void geom__in_place__TIMES_scalar(ArrayCoordinate_sp array, double mul);


};
#endif          /* ifndef SimpleVectorCoordinate_H */
