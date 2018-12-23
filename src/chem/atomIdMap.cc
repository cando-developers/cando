/*
    File: atomIdMap.cc
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

#include <clasp/core/common.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/atomId.h>
#include <cando/chem/atomIdMap.h>
#include <cando/chem/atom.h>
#include <clasp/core/wrappers.h>
namespace chem
{

// ----------------------------------------------------------------------
//



class BuildFields : public AtomIdMapFunctor {
public:
  core::Vector_sp vector;
  BuildFields(core::Vector_sp v) : vector(v) {};
  virtual void operator()(const AtomId& atomId, Atom_sp atom) {
    vector->vectorPushExtend(core::Cons_O::create(translate::to_object<AtomId>::convert(atomId),atom));
  }
};

void AtomIdToAtomMap_O::fields(core::Record_sp node)
{
  // this->Base::fields(node);
  switch (node->stage()) {
  case core::Record_O::initializing: 
  case core::Record_O::loading: {
    core::Vector_sp valueVec;
    node->field(INTERN_(kw,data), valueVec );
    this->_AtomIdMap.resize(0);
    for ( size_t i(0), iEnd(cl__length(valueVec)); i<iEnd; ++i ) {
      core::Cons_sp keyValue = gc::As<core::Cons_sp>(valueVec->rowMajorAref(i));
      core::T_sp key = CONS_CAR(keyValue);
      core::T_sp value = CONS_CDR(keyValue);
      this->_AtomIdMap[translate::from_object<AtomId>(key)._v] = gc::As_unsafe<Atom_sp>(value);
    }
  }
      break;
  case core::Record_O::saving: {
    core::Vector_sp valueVec = core::core__make_vector(_lisp->_true(),
                                                       256,
                                                       true,
                                                       core::make_fixnum(0));
    size_t idx = 0;
    BuildFields buildFields(valueVec);
    this->_AtomIdMap.iterate(buildFields);
    node->field(INTERN_(kw,data), valueVec);
  }
      break;
  case core::Record_O::patching: {
    IMPLEMENT_MEF("Do patching");
  }
      break;
  }
}


CL_DEFMETHOD  Atom_sp AtomIdToAtomMap_O::lookupAtom(const AtomId& atomId) const
{
  return this->_AtomIdMap.safeLookup(atomId);
}


void AtomIdToAtomMap_O::initialize()
{_OF();
  this->Base::initialize();
}


void AtomIdToAtomMap_O::resize(int numMols)
{_OF();
  this->_AtomIdMap.resize(numMols);
}

void AtomIdToAtomMap_O::resize(int mol, int numRes)
{_OF();
  this->_AtomIdMap.resize(mol,numRes);
}

void AtomIdToAtomMap_O::resize(int mol, int res, int numAtoms)
{_OF();
  this->_AtomIdMap.resize(mol,res,numAtoms);
}
void AtomIdToAtomMap_O::set(AtomId const& atomId, Atom_sp atom)
{_OF();
  this->_AtomIdMap[atomId] = atom;
}


class WalkMap : public AtomIdMapFunctor {
public:
  core::Function_sp callback;
  WalkMap(core::Function_sp c) : callback(c) {};
  virtual void operator()(const AtomId& atomId, Atom_sp atom) {
    core::T_sp catomid = translate::to_object<AtomId>::convert(atomId);
    core::eval::funcall(callback,catomid,atom);
  }
};

CL_DEFMETHOD void AtomIdToAtomMap_O::walk(core::Function_sp callback)
{
  WalkMap walker(callback);
  this->_AtomIdMap.iterate(walker);
}

}; /* chem */
