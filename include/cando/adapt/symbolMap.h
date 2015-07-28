/*
    File: symbolMap.h
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

//
// (C) 2004 Christian E. Schafmeister
//

#ifndef SymbolMap_H
#define SymbolMap_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lispVector.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/symbol.h>
#include <cando/adapt/adaptPackage.fwd.h>


namespace adapt {
template <class OT>
  class SymbolMap : public gctools::SmallMap<core::Symbol_sp, gctools::smart_ptr<OT>> {
  typedef gctools::SmallMap<core::Symbol_sp, gctools::smart_ptr<OT>> BaseType;
 public:
  typedef pair<core::Symbol_sp,gc::smart_ptr<OT>> key_value;
 public:
  void insert2(core::Symbol_sp key, gc::smart_ptr<OT> value) {
    this->insert(key_value(key,value));
  }
    
public:
  typedef typename BaseType::iterator iterator;
  typedef typename BaseType::const_iterator const_iterator;
};
};


namespace translate {
template <typename OT>
  struct from_object<adapt::SymbolMap<OT>> {
  typedef adapt::SymbolMap<OT> DeclareType;
  DeclareType _v;
  from_object(core::T_sp o) {
    if ( core::Cons_sp co = o.asOrNull<core::Cons_O>() ) {
      for ( auto cur : (core::List_sp)(co) ) {
        core::Cons_sp apair = gc::As<core::Cons_sp>(oCar(cur));
        core::Symbol_sp key = gc::As<core::Symbol_sp>(oCar(apair));
        core::T_sp value = oCdr(apair);
        this->_v.insert2(key,value);
      }
      return;
    } else if ( o.nilp() ) {
      this->_v.clear();
      return;
    }
    SIMPLE_ERROR(BF("Could not convert object to SymbolMap, object: %s") % _rep_(o));
  }
 };

 template <typename OT>
   struct to_object<adapt::SymbolMap<OT> > {
   typedef adapt::SymbolMap<OT> GivenType;
   static core::T_sp convert(const GivenType& v) {
     // Convert the SymbolMap into an alist
     core::List_sp res = _Nil<core::T_O>();
     for ( auto it : v ) {
       res = core::Cons_O::create(core::Cons_O::create(it.first,it.second),res);
     }
     return res;
   }
 };


};
#endif
