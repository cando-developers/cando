/*
    File: energyRosettaNonbond.h
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


/*
 *	energyRosettaNonbond.h
 *
 *	Rosetta-style nonbond component using generated kernels
 */

#pragma once

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
#include "clasp/core/ql.h"
// #include "geom/render.fwd.h" // energyRosettaNonbond.h wants DisplayList needs render.fwd.h

// Kernel parameter/term definitions
#include "cando/chem/energyKernels/rosetta_nonbond_dd_cutoff-params.h"
#include "cando/chem/energyKernels/rosetta_nonbond_dd_cutoff.h"

namespace chem {
  FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
  FORWARD(EnergyRosettaNonbond);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);
  FORWARD(Matter);

  /*! A Rosetta nonbond term
   */
  class EnergyRosettaNonbond : public EnergyTerm {
  public:
    Atom_sp              _Atom1_enb;
    Atom_sp              _Atom2_enb;
    rosetta_nonbond_term term;
  public:
    string className() { return "EnergyRosettaNonbond"; };
    Atom_sp getAtom1() { return this->_Atom1_enb; };
    Atom_sp getAtom2() { return this->_Atom2_enb; };
    bool defineForAtomPair(core::T_sp forceField,
                           Atom_sp a1,
                           Atom_sp a2,
                           size_t a1CoordinateIndexTimes3,
                           size_t a2CoordinateIndexTimes3,
                           EnergyRosettaNonbond_sp energyRosettaNonbond,
                           core::HashTable_sp atomTypes,
                           core::T_sp keepInteraction,
                           const rosetta_nonbond_parameters& params);

  public:
    core::List_sp encode() const;
    void decode(core::List_sp alist);
  };
};

namespace translate {

  template <>
  struct to_object<chem::EnergyRosettaNonbond>
  {
    typedef core::Cons_sp ExpectedType;
    typedef core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyRosettaNonbond& nonbond)
    {
      return nonbond.encode();
    }
  };

  template <>
  struct from_object<chem::EnergyRosettaNonbond>
  {
    typedef chem::EnergyRosettaNonbond ExpectedType;
    typedef ExpectedType DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
      SIMPLE_ERROR("Implement me");
    }
  };
};

namespace chem {

  class EnergyRosettaNonbond_O : public EnergyComponent_O
  {
    LISP_CLASS(chem, ChemPkg, EnergyRosettaNonbond_O, "EnergyRosettaNonbond", EnergyComponent_O);

  public:
    virtual bool restraintp() const override { return false; };
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);

  public: // virtual functions inherited from Object
    void initialize();

  public:
    typedef EnergyRosettaNonbond TermType;

  public: // instance variables
    gctools::Vec0<TermType> _Terms;
    AtomTable_sp            _AtomTable;
    core::T_sp              _NonbondForceField;
    core::HashTable_sp      _AtomTypes;
    core::T_sp              _KeepInteractionFactory;
    core::T_sp              _DisplacementBuffer;
    // If _Matter1 and _Matter2 are defined, build pair-list between them instead of using _AtomTable
    core::T_sp              _Matter1;
    core::T_sp              _Matter2;

    // Rosetta parameters (used to construct terms)
    rosetta_nonbond_parameters      _Parameters;

  public:
    virtual std::string implementation_details() const;
    virtual std::string descriptionOfContents() const;
    typedef gctools::Vec0<TermType>::iterator iterator;
    static EnergyRosettaNonbond_sp make(EnergyFunction_sp energyFunction,
                                           core::T_sp keepInteractionFactory,
                                           SetupAccumulator& setupAcc);
  public:
    CL_DEFMETHOD virtual size_t numberOfTerms() { return this->_Terms.size(); };

  public:

  public: // for building the pairList
    // In energyRosettaNonbond.h:
    double rpairlist() const { return _Parameters.rpairlist; }
    double rcut() const { return _Parameters.rcut; }
    AtomTable_sp atomTable() const { return _AtomTable; }
    CL_DEFMETHOD core::T_sp keepInteractionFactory() const { return this->_KeepInteractionFactory; };
    CL_DEFMETHOD core::T_sp matter1() const { return _Matter1; }
    CL_DEFMETHOD core::T_sp matter2() const { return _Matter2; }
    CL_DEFMETHOD void setMatter1(core::T_sp matter) { this->_Matter1 = matter; };
    CL_DEFMETHOD void setMatter2(core::T_sp matter) { this->_Matter2 = matter; };
    CL_DEFMETHOD void setMatters(core::T_sp matter1, core::T_sp matter2 ) {
      this->_Matter1 = matter1;
      this->_Matter2 = matter2;
    }
    void clearTerms() { _Terms.clear(); }
    void setDisplacementBuffer(NVector_sp buf) { _DisplacementBuffer = buf; }
    core::T_sp displacementBuffer() const { return _DisplacementBuffer; }

    bool tryAddTerm(Atom_sp a1, Atom_sp a2, size_t i3x1, size_t i3x2,
                    core::T_sp keepInteraction) {
      EnergyRosettaNonbond term;
      if (term.defineForAtomPair(_NonbondForceField, a1, a2, i3x1, i3x2,
                                 this->asSmartPtr(), _AtomTypes,
                                 keepInteraction, _Parameters)) {
        addTerm(term);
        return true;
      }
      return false;
    }

  public:
    void addTerm(const TermType& term);
    virtual void dumpTerms(core::HashTable_sp atomTypes);

    virtual double evaluateAllComponent(ScoringFunction_sp scorer,
                                        NVector_sp pos,
                                        core::T_sp energyScale,
                                        core::T_sp componentEnergy,
                                        bool calcForce,
                                        gc::Nilable<NVector_sp> force,
                                        bool calcDiagonalHessian,
                                        bool calcOffDiagonalHessian,
                                        gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                                        gc::Nilable<NVector_sp> hdvec,
                                        gc::Nilable<NVector_sp> dvec,
                                        core::T_sp activeAtomMask,
                                        core::T_sp debugInteractions);

    double debugAllComponent(ScoringFunction_sp scorer,
                             NVector_sp pos,
                             core::T_sp energyScale,
                             core::T_sp componentEnergy,
                             bool calcForce,
                             gc::Nilable<NVector_sp> force,
                             bool calcDiagonalHessian,
                             bool calcOffDiagonalHessian,
                             gc::Nilable<AbstractLargeSquareMatrix_sp> hessian,
                             gc::Nilable<NVector_sp> hdvec,
                             gc::Nilable<NVector_sp> dvec,
                             core::T_sp activeAtomMask,
                             core::T_sp debugInteractions);

    CL_DEFMETHOD void constructNonbondTermsBetweenMatters(Matter_sp mat1, Matter_sp mat2,
                                                           EnergyFunction_sp energyFunction,
                                                           core::T_sp keepInteractionFactory);

    core::T_mv maybeRebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairListBetweenMatters(core::T_sp tcoordinates);

    virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m,
                                          core::T_sp activeAtomMask );

    virtual void atomsForEachTerm(core::Function_sp callback);

    EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  public:
    EnergyRosettaNonbond_O(const EnergyRosettaNonbond_O& ss); //!< Copy constructor

    EnergyRosettaNonbond_O() :
        _KeepInteractionFactory(nil<core::T_O>()),
        _DisplacementBuffer(nil<core::T_O>()),
        _Matter1(nil<core::T_O>()),
        _Matter2(nil<core::T_O>())
    {};
  };

};
