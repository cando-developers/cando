/*
    File: indexedObjectBag.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/environment.h>
#include <clasp/core/multipleValues.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/wrappers.h>

namespace adapt {

using namespace core;

void IndexedObjectBag_O::initialize() {
  this->_Bag = HashTableEq_O::create_default();
}

CL_LISPIFY_NAME("set_next_available_entry");
CL_DEFMETHOD Fixnum IndexedObjectBag_O::set_next_available_entry(core::T_sp obj) {
  Fixnum fi = this->_NextIndex;
  ++this->_NextIndex;
  this->_Bag->setf_gethash(core::clasp_make_fixnum(fi),obj);
  return fi;
}

CL_LISPIFY_NAME("get_entry");
CL_DEFMETHOD core::T_sp IndexedObjectBag_O::entry(Fixnum index)
{
  return this->_Bag->gethash(core::clasp_make_fixnum(index));
}

CL_LISPIFY_NAME("remove_entry");
CL_DEFMETHOD void IndexedObjectBag_O::remove_entry(Fixnum index)
{
  this->_Bag->remhash(core::clasp_make_fixnum(index));
}

#define ARGS_chem__makeIndexedObjectBag "()"
#define DECL_chem__makeIndexedObjectBag ""
#define DOCS_chem__makeIndexedObjectBag "makeIndexedObjectBag"
CL_DEFUN IndexedObjectBag_sp chem__makeIndexedObjectBag()
{
  IndexedObjectBag_sp obj = IndexedObjectBag_O::create();
  return obj;
}


};
