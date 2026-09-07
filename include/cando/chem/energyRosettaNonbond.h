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
#include <cstdint>
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
  FORWARD(RosettaNonbondTermCache);
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

/*! Shared Rosetta nonbond coefficients for one AtomTable NB slot generation.
 *
 * The table is keyed by the shared radius/epsilon slots on AtomTable plus the coefficient-shaping
 * parameters. rep_weight is intentionally not part of the key because it is applied at evaluation
 * time and does not enter rosetta_nonbond_term construction.
 */
class RosettaNonbondTermCache_O : public core::CxxObject_O {
  LISP_CLASS(chem, ChemPkg, RosettaNonbondTermCache_O, "RosettaNonbondTermCache", core::CxxObject_O);
public:
  size_t _Generation = 0;
  size_t _NTypeSlots = 0;
  double _RSwitch = 0.0;
  double _RCut = 0.0;
  gctools::Vec0<rosetta_nonbond_term> _Terms;
  gctools::Vec0<char> _TermValid;

  RosettaNonbondTermCache_O() : _Terms(true), _TermValid(true) {}

  bool matches(size_t generation, size_t nTypeSlots,
               const rosetta_nonbond_parameters& params) const {
    return this->_Generation == generation
        && this->_NTypeSlots == nTypeSlots
        && this->_RSwitch == params.rswitch
        && this->_RCut == params.rcut;
  }
};

/*! Component-specific Rosetta nonbond pair record.
 */
struct RosettaNonbondPair {
  uint32_t i3x1;
  uint32_t i3x2;
  uint32_t cacheIndex;
};
static_assert(sizeof(RosettaNonbondPair) == 12, "RosettaNonbondPair must remain a compact triple");

// vdW combining rule: two atoms' (radius,epsilon) -> Lennard-Jones A/C.
// Returns false when the pair contributes nothing (epsilon == 0, e.g. polar H).
inline bool combineNonbondParams(double r1, double e1, double r2, double e2,
                                 double& parmA, double& parmC) {
  double rstar     = r1 + r2;
  double epsilonij = std::sqrt(e1 * e2);
  if (epsilonij == 0.0) return false;
  double r6  = rstar*rstar*rstar*rstar*rstar*rstar;
  parmA = epsilonij * (r6*r6);   
  parmC = 2.0 * epsilonij * r6;
  return true;
} 

class EnergyRosettaNonbond_O : public EnergyPairlistComponent_O
{
  LISP_CLASS(chem, ChemPkg, EnergyRosettaNonbond_O, "EnergyRosettaNonbond", EnergyPairlistComponent_O);

public:
  virtual bool restraintp() const override { return false; };
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

public: // virtual functions inherited from Object
  void initialize();

public:
  typedef RosettaNonbondPair TermType;

public: // instance variables
  gctools::Vec0<TermType> _Terms;
  AtomTable_sp            _AtomTable;
  core::T_sp              _NonbondForceField;
  core::HashTable_sp      _AtomTypes;
  double                  _LastFaRep = 0.0;
  
  // Rosetta parameters (used to construct terms)
  rosetta_nonbond_parameters      _Parameters;
  core::T_sp                      _CachedForAtomTable;
  //! Which generation of the atom table's SHARED _NBTypeSlot this component uses.
  //! (size_t)-1 is "never cached" - see AtomTable_O::_NBGeneration for why the atom-table pointer
  //! alone cannot answer this.
  size_t                          _CachedNBGeneration = (size_t)-1;

  RosettaNonbondTermCache_sp            _ParameterCache;
  size_t                                _PairCacheGeneration = (size_t)-1;

public:
  CL_DEFMETHOD double getLastFaRep() const { return this->_LastFaRep; };
public:
  void ensureParameterCache();
  void invalidateParameterCache() {
    this->_CachedForAtomTable = nil<core::T_O>();
    this->_CachedNBGeneration = (size_t)-1;
    this->_ParameterCache = nil<RosettaNonbondTermCache_O>();
    this->_PairCacheGeneration = (size_t)-1;
  }

  // Cached hot-path term add: no getType, no find-type, no gethash,
  // and no rosetta_nonbond_term construction — the term for this type
  // pair was precomputed in ensureParameterCache.
  bool tryAddTermCached(Atom_sp a1, Atom_sp a2, size_t li, size_t lj,
                        size_t i3x1, size_t i3x2, core::T_sp keepInteraction);
  
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
  CL_DEFMETHOD void setRepWeight(double rep_weight) { this->_Parameters.rep_weight = rep_weight; };
  CL_DEFMETHOD double getRepWeight() const { return this->_Parameters.rep_weight; };
  
public: // for building the pairList
  // In energyRosettaNonbond.h:
  CL_DEFMETHOD double rpairlist() const { return _Parameters.rpairlist; }
  double rcut() const { return _Parameters.rcut; }
  AtomTable_sp atomTable() const { return _AtomTable; }
  // setMatters is inherited from EnergyPairlistComponent_O.
  void clearTerms() { _Terms.clear(); }

  bool tryAddTerm(Atom_sp a1, Atom_sp a2, size_t i3x1, size_t i3x2,
                  core::T_sp keepInteraction) {
    this->ensureParameterCache();
    return this->tryAddTermCached(a1, a2, i3x1 / 3, i3x2 / 3,
                                  i3x1, i3x2, keepInteraction);
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
      _CachedForAtomTable(nil<core::T_O>()),
      _ParameterCache(nil<RosettaNonbondTermCache_O>())
  {};
};

};
