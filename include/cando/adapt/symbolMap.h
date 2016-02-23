/*
    File: symbolMap.h
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
