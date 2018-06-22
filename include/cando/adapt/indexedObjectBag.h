/*
    File: indexedObjectBag.h
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

#ifndef indexed_object_bag_H //[
#define indexed_object_bag_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/hashTableEq.h>
#include <cando/adapt/adaptPackage.fwd.h>
namespace adapt {

  SMART(IndexedObjectBag);

  class IndexedObjectBag_O : public core::CxxObject_O {
  LISP_CLASS(adapt, AdaptPkg, IndexedObjectBag_O, "IndexedObjectBag",core::CxxObject_O);
  void initialize();

public:
GCPRIVATE:
  gc::Nilable<core::HashTableEq_sp> _Bag;
  gc::Fixnum         _NextIndex;
public:

  gc::Fixnum set_next_available_entry(core::T_sp obj);
  core::T_sp entry(gc::Fixnum index);
  void remove_entry(gc::Fixnum index);

 public:
 IndexedObjectBag_O() : _NextIndex(0) {};
  virtual ~IndexedObjectBag_O() {};
};
};

TRANSLATE(adapt::IndexedObjectBag_O);
#endif //]
