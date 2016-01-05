/*
    File: indexedObjectBag.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/environment.h>
#include <clasp/core/multipleValues.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/wrappers.h>

namespace adapt {

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

#define ARGS_chem_makeIndexedObjectBag "()"
#define DECL_chem_makeIndexedObjectBag ""
#define DOCS_chem_makeIndexedObjectBag "makeIndexedObjectBag"
IndexedObjectBag_sp chem_makeIndexedObjectBag()
{
  IndexedObjectBag_sp obj = IndexedObjectBag_O::create();
  return obj;
}

void IndexedObjectBag_O::exposeCando(Lisp_sp lisp) {
  class_<IndexedObjectBag_O>()
    .def("set_next_available_entry", &adapt::IndexedObjectBag_O::set_next_available_entry)
    .def("get_entry", &adapt::IndexedObjectBag_O::entry)
    .def("remove_entry", &adapt::IndexedObjectBag_O::remove_entry)
    ;
  core::af_def(AdaptPkg,"make-indexed-object-bag",&chem_makeIndexedObjectBag);
}

void IndexedObjectBag_O::exposePython(Lisp_sp lisp) {
  _G();
#ifdef USEBOOSTPYTHON //[
  PYTHON_CLASS(AdaptPkg, IndexedObjectBag, "", "", _lisp)
    ;
#endif
}
EXPOSE_CLASS(adapt, IndexedObjectBag_O);
};
