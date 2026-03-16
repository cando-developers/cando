/*
    File: energyFlatBottomAnchorRestraint.h
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

#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>

// Kernel-generated term definition
#include "cando/chem/energyKernels/flat_bottom_anchor_restraint.h"

namespace chem {

FORWARD(EnergyFunction);
FORWARD(EnergyFlatBottomAnchorRestraint);
FORWARD(AtomTable);

/*! Wrapper that pairs a kernel term with its source Atom */
class EnergyFlatBottomAnchorRestraint : public EnergyTerm {
public:
  string className() { return "EnergyFlatBottomAnchorRestraint"; };

  flat_bottom_anchor_restraint_term term;
  Atom_sp _Atom1;

  Atom_sp getAtom() { return this->_Atom1; };

  core::List_sp encode() const;
  EnergyFlatBottomAnchorRestraint() : _Atom1(nil<Atom_O>()) {}
  EnergyFlatBottomAnchorRestraint(Atom_sp atom, double ka, double R0,
                                   Vector3 center, size_t I3)
    : _Atom1(atom),
      term(ka, R0, R0*R0, center.getX(), center.getY(), center.getZ(), I3)
  {}
  virtual ~EnergyFlatBottomAnchorRestraint() {}
};

} // namespace chem

namespace translate {
  template <>
  struct to_object<chem::EnergyFlatBottomAnchorRestraint> {
    typedef core::Cons_sp ExpectedType;
    typedef core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyFlatBottomAnchorRestraint& tt) {
      return tt.encode();
    }
  };

  template <>
  struct from_object<chem::EnergyFlatBottomAnchorRestraint> {
    typedef chem::EnergyFlatBottomAnchorRestraint ExpectedType;
    typedef ExpectedType DeclareType;
    DeclareType _v;
    from_object(core::T_sp o) { SIMPLE_ERROR("Implement me"); }
  };
} // namespace translate

namespace chem {

class EnergyFlatBottomAnchorRestraint_O : public EnergyComponent_O {
  LISP_CLASS(chem, ChemPkg, EnergyFlatBottomAnchorRestraint_O,
             "EnergyFlatBottomAnchorRestraint", EnergyComponent_O);

public:
  void initialize();
  virtual bool restraintp() const override { return true; };
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

public:
  typedef EnergyFlatBottomAnchorRestraint TermType;

public:
  gctools::Vec0<TermType> _Terms;

public:
  typedef gctools::Vec0<TermType>::iterator iterator;
  iterator begin() { return this->_Terms.begin(); };
  iterator end() { return this->_Terms.end(); };

  virtual std::string descriptionOfContents() const;
  virtual size_t numberOfTerms() { return this->_Terms.size(); };

public:
  void addTerm(const TermType& term);

  CL_DEFMETHOD
  size_t addFlatBottomAnchorRestraintTerm(
      EnergyFunction_sp energyFunction, Atom_sp a1,
      Vector3 center, double R0, double forceConstant);

  virtual void dumpTerms(core::HashTable_sp atomTypes);
  virtual bool is_restraint() const { return true; };

  virtual void setupHessianPreconditioner(
      NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m,
      core::T_sp activeAtomMask);

  virtual double evaluateAllComponent(
      ScoringFunction_sp scorer, NVector_sp pos,
      core::T_sp energyScale, core::T_sp componentEnergy,
      bool calcForce, gc::Nilable<NVector_sp> force,
      bool calcDiagonalHessian, bool calcOffDiagonalHessian,
      gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
      gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec,
      core::T_sp activeAtomMask, core::T_sp debugInteractions);

  EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory,
                                 SetupAccumulator& setupAcc);

public:
  EnergyFlatBottomAnchorRestraint_O(const EnergyFlatBottomAnchorRestraint_O& ss);
  DEFAULT_CTOR_DTOR(EnergyFlatBottomAnchorRestraint_O);
};

} // namespace chem