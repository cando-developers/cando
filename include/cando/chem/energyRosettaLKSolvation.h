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
#include <cstdint>
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
  FORWARD(RosettaLKTermCache);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);

  /*! The immutable coefficient table shared by every LK component over one AtomTable.
   *
   * There is deliberately no lock here.  A blueprint and all of its energy components are owned
   * by one worker thread; concurrently constructing caches on the same AtomTable is unsupported.
   */
  class RosettaLKTermCache_O : public core::CxxObject_O {
    LISP_CLASS(chem, ChemPkg, RosettaLKTermCache_O, "RosettaLKTermCache", core::CxxObject_O);
  public:
    size_t _Generation = 0;
    size_t _NTypeSlots = 0;
    double _C0 = 0.0;
    double _C1 = 0.0;
    double _RSolvLow = 0.0;
    double _RSolvHigh = 0.0;
    gctools::Vec0<rosetta_lk_solvation_term> _Terms;

    RosettaLKTermCache_O() : _Terms(true) {}

    bool matches(size_t generation, size_t nTypeSlots,
                 const rosetta_lk_solvation_parameters& params) const {
      return this->_Generation == generation
          && this->_NTypeSlots == nTypeSlots
          && this->_C0 == params.c0
          && this->_C1 == params.c1
          && this->_RSolvLow == params.r_solv_low
          && this->_RSolvHigh == params.r_solv_high;
    }
  };

  /*! The component-specific portion of an LK pair.
   *
   * CACHEINDEX selects the immutable coefficients in the AtomTable's shared cache.  The atom
   * pointers are recoverable as AtomTable->_Atoms[I3X/3].atom(), so retaining them here would be
   * redundant as well.
   */
  struct RosettaLKPair {
    uint32_t i3x1;
    uint32_t i3x2;
    uint32_t cacheIndex;
  };
  static_assert(sizeof(RosettaLKPair) == 12, "RosettaLKPair must remain a compact triple");

  /*! A standalone Rosetta LK solvation term retained for the legacy construction interface.
   */
  class EnergyRosettaLKSolvation : public EnergyTerm {
  public:
    Atom_sp                  _Atom1_enb;
    Atom_sp                  _Atom2_enb;
    uint32_t                 _I3x1;
    uint32_t                 _I3x2;
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

  class EnergyRosettaLKSolvation_O : public EnergyPairlistComponent_O
  {
    LISP_CLASS(chem, ChemPkg, EnergyRosettaLKSolvation_O, "EnergyRosettaLKSolvation", EnergyPairlistComponent_O);

  public:
    virtual bool restraintp() const override { return false; };
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);

  public: // virtual functions inherited from Object
    void initialize();

  public:
    typedef RosettaLKPair TermType;

  public: // instance variables
    gctools::Vec0<TermType>     _Terms;
    AtomTable_sp                _AtomTable;
    core::T_sp                  _LKSolvationForceField;
    core::HashTable_sp          _AtomTypes;
    // Rosetta parameters (used to construct terms)
    rosetta_lk_solvation_parameters _Parameters;

    //! Shared immutable type-pair coefficients; owned by the AtomTable cache bank.
    RosettaLKTermCache_sp _ParameterCache;
    //! The AtomTable slot generation used to construct the compact pair cache indices.
    size_t _PairCacheGeneration = (size_t)-1;

  public:
    void ensureParameterCache();   // defined in the .cc
    void invalidateParameterCache() {
      this->_ParameterCache = nil<RosettaLKTermCache_O>();
    }
    bool tryAddTermCached(Atom_sp a1, Atom_sp a2, size_t li, size_t lj,
                          size_t i3x1, size_t i3x2, core::T_sp keepInteraction);
  public:
    // pairList.h duck-typed interface
    CL_DEFMETHOD double rpairlist() const { return _Parameters.rpairlist; }
    double rcut() const { return _Parameters.r_solv_high; }  // <-- note: r_solv_high, not rcut
    AtomTable_sp atomTable() const { return _AtomTable; }
    // CL_DEFMETHOD void setMatter1(core::T_sp matter) { this->_Matter1 = matter; };
    //CL_DEFMETHOD void setMatter2(core::T_sp matter) { this->_Matter2 = matter; };
    // setMatters is inherited from EnergyPairlistComponent_O.
    void clearTerms() { _Terms.clear(); }

    bool tryAddTerm(Atom_sp a1, Atom_sp a2, size_t i3x1, size_t i3x2,
                    core::T_sp keepInteraction);
  public:
    virtual std::string implementation_details() const;
    virtual std::string descriptionOfContents() const;
    typedef gctools::Vec0<TermType>::iterator iterator;
    static EnergyRosettaLKSolvation_sp make(EnergyFunction_sp energyFunction,
                                            core::T_sp keepInteractionFactory,
                                            SetupAccumulator& setupAcc );
  public:
    CL_DEFMETHOD virtual size_t numberOfTerms() { return this->_Terms.size(); };
    void callForEachTerm(core::Function_sp callback);

  public:

  public:
    void addTerm(const TermType& term);
    virtual void dumpTerms(core::HashTable_sp atomTypes);
  virtual void atomsForEachTerm(core::Function_sp callback);

    virtual void setupHessianPreconditioner(NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m, core::T_sp activeAtomMask );

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

    /*! Point this component at a pair of matters and (re)initialize everything the
        between-matters pair-list build needs.  Mirrors
        EnergyRosettaNonbond_O::constructNonbondTermsBetweenMatters.  Invalidating the
        pair list is essential, not incidental: maybeRebuildPairList decides whether to
        rebuild from coordinate drift alone, and the MATTERS change without a single
        coordinate moving - so without the reset a second pair evaluated against an
        unchanged coordinate vector would silently reuse the previous pair's terms. */
    CL_DEFMETHOD void constructNonbondTermsBetweenMatters(Matter_sp mat1, Matter_sp mat2,
                                                   EnergyFunction_sp energyFunction,
                                                   core::T_sp keepInteractionFactory);

    core::T_mv maybeRebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairList(core::T_sp tcoordinates);

    EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  public:
    EnergyRosettaLKSolvation_O(const EnergyRosettaLKSolvation_O& ss); //!< Copy constructor

    EnergyRosettaLKSolvation_O() :
        _Terms(true),
        _ParameterCache(nil<RosettaLKTermCache_O>())
    {};
  };

};
