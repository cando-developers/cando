/*
    File: energyRosettaLKSolvation.h
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
 *	energyRosettaLKSolvation.h
 *
 *	Rosetta-style LK solvation component using generated kernels
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
// #include "geom/render.fwd.h" // energyRosettaLKSolvation.h wants DisplayList needs render.fwd.h

// Kernel parameter/term definitions
#include "cando/chem/energyKernels/rosetta_lk_solvation-params.h"
#define pi 3.14159265358979
#include "cando/chem/energyKernels/rosetta_lk_solvation.h"
#undef pi

namespace chem {
  FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
  FORWARD(EnergyRosettaLKSolvation);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);

  /*! A Rosetta LK solvation term
   */
  class EnergyRosettaLKSolvation : public EnergyTerm {
  public:
    Atom_sp                  _Atom1_enb;
    Atom_sp                  _Atom2_enb;
    rosetta_lk_solvation_term term;
  public:
    string className() { return "EnergyRosettaLKSolvation"; };
    Atom_sp getAtom1() { return this->_Atom1_enb; };
    Atom_sp getAtom2() { return this->_Atom2_enb; };
    bool defineForAtomPair(core::T_sp forceField,
                           Atom_sp a1,
                           Atom_sp a2,
                           size_t a1CoordinateIndexTimes3,
                           size_t a2CoordinateIndexTimes3,
                           EnergyRosettaLKSolvation_sp energyRosettaLKSolvation,
                           core::HashTable_sp atomTypes,
                           core::T_sp keepInteraction,
                           const rosetta_lk_solvation_parameters& params);

  public:
    core::List_sp encode() const;
    void decode(core::List_sp alist);
  };
};

namespace translate {

  template <>
  struct to_object<chem::EnergyRosettaLKSolvation>
  {
    typedef core::Cons_sp ExpectedType;
    typedef core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyRosettaLKSolvation& lk_solv)
    {
      return lk_solv.encode();
    }
  };

  template <>
  struct from_object<chem::EnergyRosettaLKSolvation>
  {
    typedef chem::EnergyRosettaLKSolvation ExpectedType;
    typedef ExpectedType DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
      SIMPLE_ERROR("Implement me");
    }
  };
};

namespace chem {

  class EnergyRosettaLKSolvation_O : public EnergyComponent_O
  {
    LISP_CLASS(chem, ChemPkg, EnergyRosettaLKSolvation_O, "EnergyRosettaLKSolvation", EnergyComponent_O);

  public:
    virtual bool restraintp() const override { return false; };
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);

  public: // virtual functions inherited from Object
    void initialize();

  public:
    typedef EnergyRosettaLKSolvation TermType;

  public: // instance variables
    gctools::Vec0<TermType>     _Terms;
    AtomTable_sp                _AtomTable;
    core::T_sp                  _NonbondForceField;
    core::HashTable_sp          _AtomTypes;
    core::T_sp                  _KeepInteractionFactory;
    core::T_sp                  _DisplacementBuffer;

    // Rosetta parameters (used to construct terms)
    rosetta_lk_solvation_parameters _Parameters;

  public:
    virtual std::string implementation_details() const;
    virtual std::string descriptionOfContents() const;
    typedef gctools::Vec0<TermType>::iterator iterator;
    static EnergyRosettaLKSolvation_sp make(AtomTable_sp atomTable, core::T_sp nbForceField,
                                            core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory,
                                            SetupAccumulator& setupAcc,
                                            core::T_sp tcoordinates = nil<core::T_O>());
  public:
    CL_DEFMETHOD virtual size_t numberOfTerms() { return this->_Terms.size(); };
    void callForEachTerm(core::Function_sp callback);

  public:
    CL_DEFMETHOD core::T_sp keepInteractionFactory() const { return this->_KeepInteractionFactory; };

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

    core::T_mv maybeRebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairList(core::T_sp tcoordinates);

    EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  public:
    EnergyRosettaLKSolvation_O(const EnergyRosettaLKSolvation_O& ss); //!< Copy constructor

    EnergyRosettaLKSolvation_O() :
        _KeepInteractionFactory(nil<core::T_O>()),
        _DisplacementBuffer(nil<core::T_O>())
    {};
  };

};
