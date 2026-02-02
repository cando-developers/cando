/*
    File: energyRosettaElec.h
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
 *	energyRosettaElec.h
 *
 *	Rosetta-style electrostatics component using generated kernels
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
// #include "geom/render.fwd.h" // energyRosettaElec.h wants DisplayList needs render.fwd.h

// Kernel parameter/term definitions
#include "cando/chem/energyKernels/rosetta_elec_cutoff-params.h"
#include "cando/chem/energyKernels/rosetta_elec_cutoff.h"

namespace chem {
  FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
  FORWARD(EnergyRosettaElec);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);

  /*! A Rosetta electrostatics term
   */
  class EnergyRosettaElec : public EnergyTerm {
  public:
    Atom_sp            _Atom1_enb;
    Atom_sp            _Atom2_enb;
    rosetta_elec_term  term;
  public:
    string className() { return "EnergyRosettaElec"; };
    Atom_sp getAtom1() { return this->_Atom1_enb; };
    Atom_sp getAtom2() { return this->_Atom2_enb; };
    bool defineForAtomPair(core::T_sp forceField,
                           Atom_sp a1,
                           Atom_sp a2,
                           size_t a1CoordinateIndexTimes3,
                           size_t a2CoordinateIndexTimes3,
                           EnergyRosettaElec_sp energyRosettaElec,
                           core::HashTable_sp atomTypes,
                           core::T_sp keepInteraction,
                           const rosetta_elec_parameters& params);

  public:
    core::List_sp encode() const;
    void decode(core::List_sp alist);
  };
};

namespace translate {

  template <>
  struct to_object<chem::EnergyRosettaElec>
  {
    typedef core::Cons_sp ExpectedType;
    typedef core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyRosettaElec& elec)
    {
      return elec.encode();
    }
  };

  template <>
  struct from_object<chem::EnergyRosettaElec>
  {
    typedef chem::EnergyRosettaElec ExpectedType;
    typedef ExpectedType DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
      SIMPLE_ERROR("Implement me");
    }
  };
};

namespace chem {
#include <cando/chem/energyKernels/rosetta_elec_cutoff-params.h>
};

namespace chem {


  class EnergyRosettaElec_O : public EnergyComponent_O
  {
    LISP_CLASS(chem, ChemPkg, EnergyRosettaElec_O, "EnergyRosettaElec", EnergyComponent_O);

  public:
    virtual bool restraintp() const override { return false; };
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);

  public: // virtual functions inherited from Object
    void initialize();

  public:
    typedef EnergyRosettaElec TermType;

  public: // instance variables
    gctools::Vec0<TermType> _Terms;
    AtomTable_sp            _AtomTable;
    core::T_sp              _NonbondForceField;
    core::HashTable_sp      _AtomTypes;
    core::T_sp              _KeepInteractionFactory;
    core::T_sp              _DisplacementBuffer;
    rosetta_elec_parameters _Parameters;
    // Rosetta parameters (used to construct terms)
    double                  _Nonbond_r_pairlist;
    size_t                  _InteractionsKept;
    size_t                  _InteractionsDiscarded;

  public:
    virtual std::string implementation_details() const;
    virtual std::string descriptionOfContents() const;
    typedef gctools::Vec0<TermType>::iterator iterator;

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

    void constructNonbondTermsFromAtomTable(AtomTable_sp atomTable,
                                            core::T_sp nbforceField,
                                            core::HashTable_sp atomTypes,
                                            core::T_sp keepInteractionFactory,
                                            core::T_sp tcoordinates = nil<core::T_O>());

    core::T_mv maybeRebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairList(core::T_sp tcoordinates);

    EnergyRosettaElec_sp copyFilter(core::T_sp keepInteractionFactory);

  public:
    EnergyRosettaElec_O(const EnergyRosettaElec_O& ss); //!< Copy constructor

    EnergyRosettaElec_O() :
        _KeepInteractionFactory(nil<core::T_O>()),
        _DisplacementBuffer(nil<core::T_O>()),
        _InteractionsKept(0),
        _InteractionsDiscarded(0)
    {};
  };

};
