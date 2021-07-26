/*
    File: coordinateArray.cc
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

//
// (C) 2004 Christian E. Schafmeister
//
#include <clasp/core/common.h>
#include <clasp/core/evaluator.h>
#include <cando/geom/coordinateArray.h>
//#include "core/archiveNode.h"
//#include "core/serializerNode.h"
#include <cando/geom/ovector3.h>
#include <cando/geom/omatrix.h>
#include <clasp/core/lispStream.h>
#ifdef RENDER
#include <cando/geom/color.h>
#include <geom/render.h>
#endif
#include <clasp/core/wrappers.h>



namespace geom 
{

SYMBOL_EXPORT_SC_(GeomPkg,vectorPushExtend_Vector3);
size_t ComplexVectorCoordinate_O::vectorPushExtend_Vector3(const Vector3& newElement, size_t extension) {
  unlikely_if (!this->_Flags.fillPointerP()) noFillPointerError(_sym_vectorPushExtend_Vector3,this->asSmartPtr());
  size_t idx = this->_FillPointerOrLengthOrDummy;
  unlikely_if (idx >= this->_ArrayTotalSize) {
    if (extension <= 0) extension = 32;
    size_t new_size = this->_ArrayTotalSize+extension;
    unlikely_if (!cl::_sym_adjust_array || !cl::_sym_adjust_array->boundP()) {
      this->internalAdjustSize_(new_size);
    } else {
      core::eval::funcall(cl::_sym_adjust_array,this->asSmartPtr(),_Nil<core::T_O>(),cl::_sym_fill_pointer,core::clasp_make_fixnum(this->_FillPointerOrLengthOrDummy));
    }
  }
  (*this)[idx] = newElement;
  ++this->_FillPointerOrLengthOrDummy;
  return idx;
}



SYMBOL_EXPORT_SC_(GeomPkg,vector3);
#if 0
CL_LAMBDA(fnsize &optional vals);
CL_LISPIFY_NAME(make-coordinate-array);
CL_DEF_CLASS_METHOD SimpleVectorCoordinate_sp SimpleVectorCoordinate_O::make(core::T_sp fnsize, core::List_sp vals)
{
  if ( fnsize.notnilp() && vals.notnilp() )
  {
    SIMPLE_ERROR(BF("You can only pass size or vals"));
  }
  auto me = gctools::GC<SimpleVectorCoordinate_O>::allocate_with_default_constructor();
  if ( fnsize.notnilp() )
  {
    core::Fixnum size = core::clasp_to_fixnum(fnsize);
    me->resize(size);
  } else if ( vals.notnilp() ) {
    me->fillFromList(vals);
  }
  return me;
}
#endif




void ArrayCoordinate_fillFromList(ArrayCoordinate_sp array, core::List_sp vals)
{
  IMPLEMENT_MEF("Come up with a loop invariant way of filling array");
#if 0
  if ( this->size() != 0 ) {
    SIMPLE_ERROR(BF("Coordinate array must be empty to fillFromCons"));
  }
  if ( vals.notnilp() ) {
    for ( auto cur : vals ) {
      core::T_sp obj = core::oCar(cur);
      if ( core::cl__consp(obj) )
      {
        core::Cons_sp consVec = obj.as<core::Cons_O>();
        Vector3 vec;
        vec.fillFromCons(consVec);
        this->appendElement(vec);
      } else if ( OVector3_sp ovec = obj.asOrNull<OVector3_O>()) {
        this->appendElement(ovec->value());
      } else {
        SIMPLE_ERROR(BF("Cannot convert value[%s] to Vector3") % _rep_(obj) );
      }
    }
  }
#endif
}


void SimpleVectorCoordinate_O::zeroAllCoordinates()
{
  for (size_t i(0); i<this->length(); ++i ) {
    (*this)[i] = Vector3(0.0,0.0,0.0);
  }
}

void SimpleVectorCoordinate_O::addPoints(SimpleVectorCoordinate_sp y)
{_OF();
  int  sz;
//  ASSERT_eq(this->_Points.size(),y->size());
  sz = this->length();
  for ( size_t i=0; i<sz; i++ ) {
    (*this)[i] = (*this)[i].add((*y)[i]);
  }
}


double SimpleVectorCoordinate_O::distanceTo(SimpleVectorCoordinate_sp y)
{_OF();
  int  sz;
  double dist = 0.0;
  sz = this->length();
  ASSERT_eq(this->length(),y->length());
  ASSERT(this->length()>0);
  for ( size_t i=0; i<sz; i++ ) {
    Vector3 v = (*this)[i].sub((*y)[i]);
    double d = v.dotProduct(v);
    dist += d;
  }
  dist /= (double)(sz);
  return sqrt(dist);
}


#ifdef RENDER
Render_sp SimpleVectorCoordinate_O::asDisplayListOfPoints()
{
  int  sz;
  DisplayList_sp dl = DisplayList_O::create();
  GrPoints_sp pl = GrPoints_O::create();
  sz = this->_Points.size();
  for ( int i=0; i<sz; i++ ) 
  {
    pl->appendPoint(this->_Points[i],Color_O::white(_lisp));
  }
  dl->append(pl);
  return dl;
}
#endif

#if 0
void SimpleVectorCoordinate_O::multiplyBy( double mul )
{
  int  sz;
  sz = this->_Points.size();
  for ( int i=0; i<sz; i++ ) {
    this->_Points[i] = this->_Points[i].multiplyByScalar(mul);
  }
}
#endif
#if 0
void SimpleVectorCoordinate_O::fields(core::Record_sp node)
{
  if ( node->saving() ) {
    uint sz = this->size();
    node->attribute("size",sz);
    int i;
    core::ArchiveP vnode;
    vector<Vector3>::iterator vi;
    for ( i=0,vi=this->_Points.begin(); vi!=this->_Points.end(); i++,vi++ ) {
      vnode = node->createChildNode("V");
      vnode->attribute("i",i);
      (*vi).archive(vnode);
    }
  } else {
    uint sz=0, i=0;
    node->attribute("size",sz);
    this->resize(sz);
    vector<core::ArchiveP>::iterator	vi;
    i = 0;
    for ( vi=node->begin_Children(); vi!=node->end_Children(); vi++ ) {
      (*vi)->setRecognized(true);
      (*vi)->attribute("i",i);
      this->_Points[i].archive(*vi);
      i++;
    }
  }
}
#endif


#if 0
CL_LISPIFY_NAME("write-to-stream");
CL_DEFMETHOD     void SimpleVectorCoordinate_O::writeToStream(string const& info, core::T_sp sout)
{_OF();
  core::clasp_write_string((BF("+++SimpleVectorCoordinate %d\n") % this->_Points.size() ).str(),sout);
  core::clasp_write_string((BF("+ %s\n") % info).str(),sout);
  for (int i=0; i<(int)this->_Points.size(); i++ )
  {
    this->_Points[i].write(sout);
    core::clasp_terpri(sout);
  }
}
#endif

#if 0
core::List_sp SimpleVectorCoordinate_O::encode() {
  core::Vector_sp v = core::core__make_vector(cl::_sym_DoubleFloat_O,3*this->_Points.size());
  size_t cur(0);
  for ( size_t i(0); i<this->length(); ++i ) {
    v->rowMajorAset(cur++,core::clasp_make_double_float((*this)[i].getX()));
    v->rowMajorAset(cur++,core::clasp_make_double_float((*this)[i].getY()));
    v->rowMajorAset(cur++,core::clasp_make_double_float((*this)[i].getZ()));
  }
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,pnts),core::clasp_make_fixnum(this->_Points.size())),
                                  core::Cons_O::create(INTERN_(kw,pnts),v));
}

void SimpleVectorCoordinate_O::decode(core::List_sp c) {
  core::Fixnum_sp fn_size = gc::As<core::Fixnum_sp>(oCdr(oFirst(c)));
  core::Vector_sp v = gc::As<core::Vector_sp>(oCdr(oSecond(c)));
  this->resize(core::clasp_to_fixnum(fn_size));
  size_t cur(0);
  for ( size_t i(0); i<this->_Points.size(); ++i ) {
    // TODO: Optimize this if v is a specialized vector of doubles
    double x = core::clasp_to_double(v->rowMajorAref(cur++));
    double y = core::clasp_to_double(v->rowMajorAref(cur++));
    double z = core::clasp_to_double(v->rowMajorAref(cur++));
    Vector3 v3(x,y,z);
    this->_Points[i] = v3;
  }
}
#endif


#if 0
string SimpleVectorCoordinate_O::parseFromStream(core::Stream_sp sin)
{_OF();
  streampos spos;
  string line, info;
  int numLines;
  bool hitEof;
  spos = sin->tell();
  sin->readLine(line,hitEof);
  if ( !hitEof )
  {
    SIMPLE_ERROR(BF("Could not read header for SimpleVectorCoordinate - got[%s] eof[%d] ")
                 % line % hitEof );
  }
  if ( line.substr(0,18)=="+++SimpleVectorCoordinate")
  {
    stringstream snum;
    snum.str(line.substr(18,999));
    snum >> numLines;
    ASSERTF(numLines>=0,BF("Illegal number of entries[%d]") % numLines);
    this->resize(numLines);
    sin->readLine(line,hitEof);
//	    getline(sin->stream(),line);
    if ( !hitEof )
    {
      SIMPLE_ERROR(BF("Could not read info - got[%s]") % line );
    }
    info = line.substr(2,9999);
    for ( int i=0; i<numLines; i++ )
    {
      Vector3 pos;
      sin->readLine(line,hitEof);
//		getline(sin->stream(),line);
      if ( !hitEof )
      {
        SIMPLE_ERROR(BF("Error reading relative line[%d] contents[%s] for SimpleVectorCoordinate") % i % line );
      }
      pos.parseFromString(line);
      this->getElement(i)=pos;
    }
  } else
  {
    sin->seek(spos);
    SIMPLE_ERROR(BF("Could not read SimpleVectorCoordinate - header was[%s]") % line );
  }
  return info;
}
#endif





#if 0
SimpleVectorCoordinateWithHash_sp SimpleVectorCoordinateWithHash_O::create(uint sz)
{
  auto h = gctools::GC<SimpleVectorCoordinateWithHash_O>::allocate_with_default_constructor();
  h->resize(sz);
  return h;
}







void SimpleVectorCoordinateWithHash_O::initialize()
{
  this->Base::initialize();
  this->_Hash = UndefinedUnsignedInt;
}

#ifdef XML_ARCHIVE
void SimpleVectorCoordinateWithHash_O::archiveBase(core::ArchiveP node)
{
  this->Base::archiveBase(node);
  node->attribute("hash",this->_Hash);
}
#endif


#ifdef OLD_SERIALIZE
void SimpleVectorCoordinateWithHash_O::serialize(serialize::SNode node)
{
  this->Base::serialize(node);
  node->attribute("hash",this->_Hash);
}
#endif
#endif


};


namespace geom {



CL_LISPIFY_NAME(transform);
CL_DEFUN void geom__in_place_transform(ArrayCoordinate_sp array, const Matrix& transform) {
  core::AbstractSimpleVector_sp sv;
  size_t start;
  size_t end;
  array->asAbstractSimpleVectorRange(sv,start,end);
  SimpleVectorCoordinate_sp svc = gc::As<SimpleVectorCoordinate_sp>(sv);
  for ( size_t i(start); i<end; ++i ) {
    (*svc)[i] = transform.multiplyByVector3((*svc)[i]);
  }
};

CL_LISPIFY_NAME(transform);
CL_DEFUN void geom__in_place_mul_scalar(ArrayCoordinate_sp array, double mul) {
  core::AbstractSimpleVector_sp sv;
  size_t start;
  size_t end;
  array->asAbstractSimpleVectorRange(sv,start,end);
  SimpleVectorCoordinate_sp svc = gc::As<SimpleVectorCoordinate_sp>(sv);
  for ( size_t i(start); i<end; ++i ) {
    (*svc)[i] = (*svc)[i]*mul;
  }
};

CL_DEFUN core::T_sp geom__unit_vectors()
{
  Vector3 v;
  SimpleVectorCoordinate_sp a = SimpleVectorCoordinate_O::make(0,v,true,0,NULL);
  v.set(1.0,0.0,0.0);
  a->getElement(0) = v;
  v.set(0.0,1.0,0.0);
  a->getElement(1) = v;
  v.set(0.0,0.0,1.0);
  a->getElement(2) = v;
  return a;
}

};
