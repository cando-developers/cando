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
    typedef EnergyRosettaLKSolvation TermType;

  public: // instance variables
    gctools::Vec0<TermType>     _Terms;
    AtomTable_sp                _AtomTable;
    core::T_sp                  _LKSolvationForceField;
    core::HashTable_sp          _AtomTypes;
    // Rosetta parameters (used to construct terms)
    rosetta_lk_solvation_parameters _Parameters;

    // The term is a pure function of the two atoms' LK parameters (dGfree,
    // lambda, radius, volume), which are resolved from :lk-solvation-atom-type
    // and are therefore per-TYPE.  So collapse the distinct parameter tuples
    // into slots and precompute the term for every ordered slot pair - the term
    // is NOT symmetric in i/j (atom i desolvated by j differs from j by i), so
    // the table is full nt x nt, not triangular.
    gctools::Vec0<int>                        _TypeSlot;    // atom index -> slot, -1 = no params
    size_t                                    _NTypeSlots = 0;
    gctools::Vec0<rosetta_lk_solvation_term>  _TermCache;   // _NTypeSlots^2, row-major
    gctools::Vec0<char>                       _TermCacheValid;
    core::T_sp          _CachedForAtomTable;   // init nil in ctor
    //! The AtomTable_O::_LKGeneration this component's _TypeSlot and _TermCache were derived from.
    //! Not redundant with _CachedForAtomTable: the shared table can be rebuilt IN PLACE on the
    //! same atom table, and then the pointer still matches while the contents no longer do.
    size_t              _CachedLKGeneration = (size_t)-1;

  public:
    void ensureParameterCache();   // defined in the .cc
    void invalidateParameterCache() {
      this->_CachedForAtomTable = nil<core::T_O>();
      this->_TypeSlot.clear();
      this->_TermCache.clear();
      this->_TermCacheValid.clear();
      this->_NTypeSlots = 0;
      this->_CachedLKGeneration = (size_t)-1;
    }
    // Cached hot-path term add: no plist scan, no find-lksolvation-type funcall,
    // and no rosetta_lk_solvation_term construction - the term for this type
    // pair was precomputed in ensureParameterCache.
    bool tryAddTermCached(Atom_sp a1, Atom_sp a2, size_t li, size_t lj,
                          size_t i3x1, size_t i3x2, core::T_sp /*keepInteraction*/) {
      int s1 = this->_TypeSlot[li];
      int s2 = this->_TypeSlot[lj];
      if (s1 < 0)
        SIMPLE_ERROR("Could not find LKSolvation parameter for atom {} - property-list {}",
                     _rep_(a1), _rep_(a1->getProperties()));
      if (s2 < 0)
        SIMPLE_ERROR("Could not find LKSolvation parameter for atom {} - property-list {}",
                     _rep_(a2), _rep_(a2->getProperties()));
      size_t k = (size_t)s1 * this->_NTypeSlots + (size_t)s2;
      // ensureParameterCache fills every slot pair, so a miss here is a cache
      // construction bug - fail loudly rather than silently dropping the pair
      // (the caller has no else branch and would not count it as discarded).
      if (!this->_TermCacheValid[k])
        SIMPLE_ERROR("LKSolvation term cache miss for slot pair {},{} of {} - cache construction bug",
                     s1, s2, this->_NTypeSlots);
      TermType term;
      term._Atom1_enb = a1;
      term._Atom2_enb = a2;
      term.term = this->_TermCache[k];     // struct copy - no term arithmetic
      term.term.i3x1 = (int)i3x1;          // only the geometry-dependent fields
      term.term.i3x2 = (int)i3x2;
      this->addTerm(term);
      return true;
    }
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
                    core::T_sp keepInteraction) {
      EnergyRosettaLKSolvation term;
      if (term.defineForAtomPair(_LKSolvationForceField, a1, a2, i3x1, i3x2,
                                 this->asSmartPtr(), _AtomTypes,
                                 keepInteraction, _Parameters)) {
        addTerm(term);
        return true;
      }
      return false;
    }
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
        _CachedForAtomTable(nil<core::T_O>())
    {};
  };

};
