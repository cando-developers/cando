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
  virtual void operator()(const AtomId& atomId, core::T_sp atom) {
    vector->vectorPushExtend(core::Cons_O::create(translate::to_object<AtomId>::convert(atomId),atom));
  }
};

void AtomIdMap_O::fields(core::Record_sp node)
{
  // this->Base::fields(node);
#if 1
  node->field(INTERN_(kw,data),this->_AtomIdMap);
#else
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
      this->_AtomIdMap[translate::from_object<AtomId>(key)._v] = value;
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
    // Nothing
  }
      break;
  }
#endif
}


CL_DEFUN chem::AtomIdMap_sp make_AtomIdMap() {
  return gctools::GC<AtomIdMap_O>::allocate_with_default_constructor();
}


CL_DEFMETHOD  core::T_sp AtomIdMap_O::AtomIdMap_get(const AtomId& atomId) const
{
  return this->safeLookup(atomId);
}

CL_DEFUN core::T_sp chem__lookup_atom(AtomIdMap_sp map, const chem::AtomId& atomId) {
  return map->AtomIdMap_get(atomId);
}

void AtomIdMap_O::initialize()
{_OF();
  this->_AtomIdMap = core::ComplexVector_T_O::make(0);
  this->Base::initialize();
}


CL_DEFMETHOD void AtomIdMap_O::resizeAggregate(int numMols)
{_OF();
  size_t oldNumMols = this->_AtomIdMap->length();
  this->_AtomIdMap->resize(numMols);
  for ( size_t it = oldNumMols; it<numMols; it++ ) {
    this->_AtomIdMap->rowMajorAset(it,core::ComplexVector_T_O::make(0));
  }
}

CL_DEFMETHOD void AtomIdMap_O::resizeMolecule(int mol, int numRes)
{_OF();
  if (mol<this->_AtomIdMap->length()) {
    core::ComplexVector_T_sp residues = gc::As<core::ComplexVector_T_sp>(this->_AtomIdMap->rowMajorAref(mol));
    size_t oldNumRes = residues->length();
    residues->resize(numRes);
    for ( size_t it = oldNumRes; it<numRes; it++ ) {
      residues->rowMajorAset(it,core::ComplexVector_T_O::make(0));
    }
    return;
  }
  SIMPLE_ERROR(BF("molecule index %lu is out of range %lu") % mol % this->_AtomIdMap->length());
}

CL_DEFMETHOD void AtomIdMap_O::resizeResidue(int mol, int res, int numAtoms)
{_OF();
  if (mol<this->_AtomIdMap->length()) {
    core::ComplexVector_T_sp residues = gc::As<core::ComplexVector_T_sp>(this->_AtomIdMap->rowMajorAref(mol));
    if (res<residues->length()) {
      core::ComplexVector_T_sp atoms = gc::As<core::ComplexVector_T_sp>(residues->rowMajorAref(res));
      size_t oldNumAtoms = atoms->length();
      atoms->resize(numAtoms);
      return;
    }
    SIMPLE_ERROR(BF("residue index %lu is out of range %lu") % res % residues->length() );
  }
  SIMPLE_ERROR(BF("molecule index %lu is out of range %lu") % mol % this->_AtomIdMap->length());
}

void AtomIdMap_O::set(AtomId const& atomId, core::T_sp atom)
{_OF();
  this->safeSet( atomId, atom );
}


CL_DEFMETHOD void AtomIdMap_O::AtomIdMap_set(const AtomId& atomId, core::T_sp atom)
{
  this->safeSet( atomId, atom );
}


class WalkMap : public AtomIdMapFunctor {
public:
  core::Function_sp callback;
  WalkMap(core::Function_sp c) : callback(c) {};
  virtual void operator()(const AtomId& atomId, core::T_sp atom) {
    core::T_sp catomid = translate::to_object<AtomId>::convert(atomId);
    core::eval::funcall(callback,catomid,atom);
  }
};

CL_DEFMETHOD void AtomIdMap_O::walk(core::Function_sp callback)
{
  WalkMap walker(callback);
  this->iterate(walker);
}

}; /* chem */
