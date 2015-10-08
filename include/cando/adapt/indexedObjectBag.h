/*
    File: objectSet.h
*/

/*
Copyright (c) 2014, Christian E. Schafmeister
 
CLASP is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.
 
See directory 'clasp/licenses' for full details.
 
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
/* -^- */

#ifndef indexed_object_bag_H //[
#define indexed_object_bag_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/hashTableEq.h>
#include <cando/adapt/adaptPackage.fwd.h>
namespace adapt {
  using namespace core;

  SMART(IndexedObjectBag);

class IndexedObjectBag_O : public T_O {
  LISP_BASE1(T_O);
  LISP_CLASS(adapt, AdaptPkg, IndexedObjectBag_O, "IndexedObjectBag");
  DECLARE_INIT();
  void initialize();

public:
GCPRIVATE:
  gc::Nilable<HashTableEq_sp> _Bag;
  Fixnum         _NextIndex;
public:

  Fixnum set_next_available_entry(core::T_sp obj);
  core::T_sp entry(Fixnum index);
  void remove_entry(Fixnum index);

 public:
 IndexedObjectBag_O() : _NextIndex(0) {};
  virtual ~IndexedObjectBag_O() {};
};
};

TRANSLATE(adapt::IndexedObjectBag_O);
#endif //]
