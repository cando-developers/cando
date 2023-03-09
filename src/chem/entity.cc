#ifdef USE_TOPOLOGY
/*
    File: entity.cc
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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <cando/chem/entity.h>
#include <cando/chem/candoDatabase.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/constitution.h>
#include <clasp/core/wrappers.h>

namespace chem {

//
// Constructor
//

void Entity_O::initialize() { this->Base::initialize(); }

//
// Copy Constructor
//
Entity_O::Entity_O(const Entity_O &emr) : core::CxxObject_O(emr) {}

//
// Destructor
//

#ifdef XML_ARCHIVE
void Entity_O::archiveBase(core::ArchiveP node) { this->Base::archiveBase(node); }
#endif

CL_LISPIFY_NAME("constitution");
CL_DEFMETHOD Constitution_sp Entity_O::constitution() {
  Constitution_sp x = nil<Constitution_O>();
  return x;
}

string Entity_O::__repr__() const {
  stringstream ss;
  ss << "[" << this->className() << "/name=" << this->__repr__() << "]";
  return ss.str();
}

CL_LISPIFY_NAME("minimalRepresentativeList");
CL_DEFMETHOD RepresentativeList_sp Entity_O::minimalRepresentativeList() const {
  
  RepresentativeList_sp expandedList = this->expandedRepresentativeList();
  core::HashTableEq_sp binder = core::HashTableEq_O::create_default();
  for (int i(0), iEnd(expandedList->length()); i < iEnd; ++i) {
    RepresentedEntityNameSet_sp curNameSet = expandedList->rowMajorAref(i).as<RepresentedEntityNameSet_O>();
    core::Symbol_sp representor = curNameSet->getRepresentative();
    RepresentedEntityNameSet_sp minimalRepresentedEntityNameSet;
    core::KeyValuePair* pair = binder->find(representor);
    if (pair) {
      minimalRepresentedEntityNameSet =
          pair->_Value.as<RepresentedEntityNameSet_O>(); // binder->indexed_value(bi).as<RepresentedEntityNameSet_O>();
      minimalRepresentedEntityNameSet->mergeEntityNames(curNameSet);
    } else {
      auto temp = gctools::GC<RepresentedEntityNameSet_O>::copy(*curNameSet); // = RP_Copy<RepresentedEntityNameSet_O>(curNameSet);
      minimalRepresentedEntityNameSet = temp;
      minimalRepresentedEntityNameSet->expandToTerminalEntityNames();
      binder->setf_gethash(representor, minimalRepresentedEntityNameSet);
    }
  }
  RepresentativeList_sp minimalList = RepresentativeList_O::create();
  binder->maphash([&minimalList](core::T_sp key, core::T_sp val) {
    minimalList->vectorPushExtend(val); // binder->indexed_value(val));
  });
  return minimalList;
}

}; // namespace chem
#endif
