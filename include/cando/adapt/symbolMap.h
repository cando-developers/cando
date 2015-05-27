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
  void insert2(core::Symbol_sp key, gc::smart_ptr<OT> value) {
    this->insert(key_value(key,value));
  }
    
public:
  typedef typename BaseType::iterator iterator;
  typedef typename BaseType::const_iterator const_iterator;
};
};
#endif
