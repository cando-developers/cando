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
#include <cstdint>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
// #include "geom/render.fwd.h" // energyRosettaElec.h wants DisplayList needs render.fwd.h

#include "clasp/core/ql.h"
// Kernel parameter/term definitions
#include "cando/chem/energyKernels/rosetta_elec_cutoff-params.h"
#include "cando/chem/energyKernels/rosetta_elec_cutoff.h"

namespace chem {
FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
FORWARD(EnergyRosettaElec);
FORWARD(RosettaElecTermCache);
FORWARD(FFNonbondDb);
FORWARD(AtomTable);
FORWARD(Matter);

inline double calculate_dQ1Q2(double electrostaticScale, double electrostaticModifier, double charge1, double charge2) {
  return electrostaticScale * electrostaticModifier * charge1 * charge2;
}

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

  /*! Shared Rosetta electrostatic coefficients for kqq == 1.0.
   *
   * A pair's charge product is still per pair. Everything else in rosetta_elec_term is a
   * parameter-only linear coefficient, so one prototype can be scaled on the stack during
   * evaluation instead of storing ten doubles in every pair record.
   */
  class RosettaElecTermCache_O : public core::CxxObject_O {
    LISP_CLASS(chem, ChemPkg, RosettaElecTermCache_O, "RosettaElecTermCache", core::CxxObject_O);
  public:
    double _EpsCore = 0.0;
    double _EpsSolvent = 0.0;
    double _Rmin = 0.0;
    double _Rlow = 0.0;
    double _Rhi = 0.0;
    double _Rcut = 0.0;
    rosetta_elec_term _Prototype;

    bool matches(const rosetta_elec_parameters& params) const {
      return this->_EpsCore == params.eps_core
          && this->_EpsSolvent == params.eps_solvent
          && this->_Rmin == params.rmin
          && this->_Rlow == params.rlow
          && this->_Rhi == params.rhi
          && this->_Rcut == params.rcut;
    }
  };

  /*! Component-specific electrostatic pair record.
   *
   * kqq is the only coefficient that is genuinely per atom pair. Atom pointers are recovered from
   * the AtomTable when a caller iterates/debugs terms.
   */
  struct RosettaElecPair {
    uint32_t i3x1;
    uint32_t i3x2;
    double kqq;
  };
  static_assert(sizeof(RosettaElecPair) == 16, "RosettaElecPair must remain compact");


  class EnergyRosettaElec_O : public EnergyPairlistComponent_O
  {
    LISP_CLASS(chem, ChemPkg, EnergyRosettaElec_O, "EnergyRosettaElec", EnergyPairlistComponent_O);

  public:
    virtual bool restraintp() const override { return false; };
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);

  public: // virtual functions inherited from Object
    void initialize();

  public:
    typedef RosettaElecPair TermType;

  public: // instance variables
    gctools::Vec0<TermType> _Terms;
    AtomTable_sp            _AtomTable;
    core::T_sp              _NonbondForceField;
    core::HashTable_sp      _AtomTypes;
    rosetta_elec_parameters _Parameters;
    // Rosetta parameters (used to construct terms)

    // ---- hot-path cache -------------------------------------------------
    // Electrostatics keep kqq per pair and share the parameter-only spline coefficients.
    double                  _DQ1Q2Scale = 0.0;
    RosettaElecTermCache_sp _ParameterCache;

  public:
      void ensureParameterCache();
      void invalidateParameterCache() {
        this->_ParameterCache = nil<RosettaElecTermCache_O>();
      }
      bool tryAddTermCached(Atom_sp a1, Atom_sp a2, size_t li, size_t lj,
                            size_t i3x1, size_t i3x2, core::T_sp keepInteraction);
        
  public:
    // In energyRosettaElec.h:
    CL_DEFMETHOD double rpairlist() const { return _Parameters.rpairlist; }
    double rcut() const { return _Parameters.rcut; }
    AtomTable_sp atomTable() const { return _AtomTable; }
    // keepInteractionFactory() already exists (line 159)
    //    CL_DEFMETHOD void setMatter1(core::T_sp matter) { this->_Matter1 = matter; };
    //    CL_DEFMETHOD void setMatter2(core::T_sp matter) { this->_Matter2 = matter; };
  // setMatters is inherited from EnergyPairlistComponent_O.
    void clearTerms() { _Terms.clear(); }

    bool tryAddTerm(Atom_sp a1, Atom_sp a2, size_t i3x1, size_t i3x2,
                    core::T_sp keepInteraction) {
      this->ensureParameterCache();
      return this->tryAddTermCached(a1, a2, i3x1 / 3, i3x2 / 3,
                                    i3x1, i3x2, keepInteraction);
    }

  public:
    virtual std::string implementation_details() const;
    virtual std::string descriptionOfContents() const;
    typedef gctools::Vec0<TermType>::iterator iterator;

    static EnergyRosettaElec_sp make(EnergyFunction_sp energyFunction,
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

    CL_DEFMETHOD void constructNonbondTermsBetweenMatters(Matter_sp mat1, Matter_sp mat2,
                                                           EnergyFunction_sp energyFunction,
                                                           core::T_sp keepInteractionFactory);

    core::T_mv maybeRebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairList(core::T_sp tcoordinates);
    core::T_mv rebuildPairListBetweenMatters(core::T_sp tcoordinates);

    EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc);

  public:
    EnergyRosettaElec_O(const EnergyRosettaElec_O& ss); //!< Copy constructor

    // No initializer list left - every member it used to initialize now belongs to
    // EnergyPairlistComponent_O, whose own constructor nils them.
    EnergyRosettaElec_O() :
        _ParameterCache(nil<RosettaElecTermCache_O>())
    {};
  };

};
