/*
    File: objectSet.h
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

#ifndef ObjectSet_H //[
#define ObjectSet_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/hashTableEq.h>
#include <cando/adapt/adaptPackage.fwd.h>
namespace adapt {
  using namespace core;
class ObjectSetCartesianProductWrapper {
public:
  virtual T_sp operator()(T_sp obj1, T_sp obj2) const = 0;
};

SMART(ObjectSet);

 class ObjectSet_O : public core::General_O {
  LISP_CLASS(adapt, AdaptPkg, ObjectSet_O, "ObjectSet",core::General_O);
  void initialize();

public:
#if defined(XML_ARCHIVE)
  void archive(ArchiveP node);
#endif // defined(XML_ARCHIVE)

GCPRIVATE:
  HashTableEq_sp _Set;

public:
CL_NAME("CONTAINS");
CL_DEFMETHOD   bool contains(T_sp obj) { return this->_Set->contains(obj); };
CL_NAME("INSERT");
CL_DEFMETHOD   void insert(T_sp obj) { this->_Set->setf_gethash(obj, _lisp->_true()); };
  void remove(T_sp obj) { this->_Set->remhash(obj); };
CL_LISPIFY_NAME("size");
CL_DEFMETHOD   int size() { return this->_Set->hashTableSize(); };

  List_sp asCons();

  void addObjectsInCons(List_sp objs);

  void addObjects(ObjectSet_sp b);

  // Set theory operations

  //! A setUnion B = (x:x E A || x E B)
  ObjectSet_sp setUnion(ObjectSet_sp b);

  //! A intersection B = (x:x E A && x E B)
  ObjectSet_sp intersection(ObjectSet_sp b);

  //! A-B = (x: x E A && not x E B )
  ObjectSet_sp relativeComplement(ObjectSet_sp b);

  ObjectSet_sp removeAll(ObjectSet_sp b) { return this->relativeComplement(b); };

  /*! AxB = ("x,y": x E A ; y E B )
	  Return x,y as ObjectPair
	*/
  ObjectSet_sp cartesianProduct(ObjectSet_sp b);

  /*! AxB = ("x,y": x E A ; y E B )
	  Return x,y wrapped by a callback functor
	*/
  ObjectSet_sp cartesianProductWrapped(ObjectSet_sp b, const ObjectSetCartesianProductWrapper &wrapper);

  string asString() const;

  void map(std::function<void(T_sp)> const &fn);
  void map(std::function<void(T_sp)> const &fn) const;

public:
  DEFAULT_CTOR_DTOR(ObjectSet_O);
};
};

TRANSLATE(adapt::ObjectSet_O);
#endif //]
