/*
    File: energyNonbond.h
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
 *	energyNonbond.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyNonbond_H  //[
#define	EnergyNonbond_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyComponent.h>
//#include "geom/render.fwd.h"// energyNonbond.h wants DisplayList needs render.fwd.h

namespace       chem {
  FORWARD(EnergyFunction); // Declares class EnergyFunction_O {} and EnergyFunction_sp
  FORWARD(EnergyNonbond);
  FORWARD(FFNonbondDb);
  FORWARD(AtomTable);


struct TermNonBond
{
  INT  I1; //!< i*3 index into coordinate vector, must match Mathematica code!
  INT  I2; //!< i*3 index into coordinate vector, must match Mathematica code!
  double dQ1Q2;
  // Anonymous union/struct keeps term.dA/term.dC access while allowing
  // alternate parameterization.
  union {
    struct { // These are Amber parameters
      double dA;
      double dC;
    };
    struct { // These are Rosetta parameters
      double dEpsilon;
      double dSigma;
      double dM;
      double dB;
    };
  };
};


inline	string	XmlTag_Nonbond() { return "Nonbond"; };
inline	string	XmlTag_EnergyNonbond() { return "EnergyNonbond"; };

/*! A Nonbond term, this is subclassed in AmberFunction to deal with Amber specific stuff.
 * I should combine energyFunction.cc and amberFunction.cc together.
 */
class EnergyNonbond : public EnergyTerm {
public:
  Atom_sp      	_Atom1_enb;
  Atom_sp      	_Atom2_enb;
  TermNonBond	term;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  bool		_calcForce;
  bool		_calcDiagonalHessian;
  bool		_calcOffDiagonalHessian;
#include <cando/chem/energy_functions/_Nonbond_debugEvalDeclares.cc>
#endif
public:
  string	className()	{ return "EnergyNonbond"; };
  Atom_sp	getAtom1() { return this->_Atom1_enb; };
  Atom_sp	getAtom2() { return this->_Atom2_enb; };
#if 0
  bool	defineFrom(core::T_sp	forceField,
                   bool		is14,
                   EnergyAtom	*iea1,
                   EnergyAtom	*iea2,
                   EnergyNonbond_sp nb,
                   core::HashTable_sp atomTypes );
#endif
  bool	defineForAtomPair(core::T_sp	forceField,
                          bool		is14,
                          Atom_sp        a1,
                          Atom_sp        a2,
                          size_t         a1CoordinateIndexTimes3,
                          size_t         a2CoordinateIndexTimes3,
                          EnergyNonbond_sp energyNonbond,
                          core::HashTable_sp atomTypes,
                          core::T_sp     keepInteraction
                          );

public:
public:
  adapt::QDomNode_sp	asXml();
  void	parseFromXmlRelativeToContainer( adapt::QDomNode_sp xml, Matter_sp parent );
  void	parseFromXmlUsingAtomTable(adapt::QDomNode_sp xml, AtomTable_sp atomTable );

public:
  core::List_sp encode() const;
  void decode(core::List_sp alist);
};
};

namespace translate {

  template <>
struct	to_object<chem::EnergyNonbond>
  {
    typedef	core::Cons_sp ExpectedType;
    typedef	core::Cons_sp DeclareType;
    static core::T_sp convert(const chem::EnergyNonbond& nonbond)
    {
      return nonbond.encode();
    }
  };

  template <>
    struct	from_object<chem::EnergyNonbond>
  {
    typedef	chem::EnergyNonbond	ExpectedType;
    typedef	ExpectedType 		DeclareType;
    DeclareType _v;
    from_object(core::T_sp o)
    {
      SIMPLE_ERROR("Implement me");
    }
  };
};

namespace chem {

double	_evaluateEnergyOnly_Nonbond(ScoringFunction_sp score, core::T_sp energyScale,
                                    int I1, int I2, core::T_sp activeAtomMask,
                                    num_real x1, num_real y1, num_real z1,
                                    num_real x2, num_real y2, num_real z2,
                                    num_real dA, num_real dC, num_real dQ1Q2 );


class EnergyNonbond_O : public EnergyComponent_O
{
  LISP_CLASS(chem,ChemPkg,EnergyNonbond_O,"EnergyNonbond",EnergyComponent_O);

 public:
  virtual bool restraintp() const override {return false;};
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);

 public: // virtual functions inherited from Object
  void	initialize();
//	string	__repr__() const;
 public:
  typedef EnergyNonbond TermType;
 public: // instance variables
  double	          	_EnergyElectrostatic;
  bool                          _UsesExcludedAtoms;
  double                        _Nonbond_r_switch;
  double                        _Nonbond_r_cut;
  double                        _Nonbond_r_pairlist;
  double                        _Nonbond_invdd; // Distance dependent dielectric reciprocal
  size_t                        _NonbondsKept;
  size_t                        _NonbondsDiscarded;
  gctools::Vec0<TermType>	_Terms14;
  // Pairlist
  AtomTable_sp                  _AtomTable;
  core::T_sp                    _NonbondForceField;
  core::HashTable_sp            _AtomTypes;
  core::T_sp                    _KeepInteractionFactory;
  core::T_sp                    _DisplacementBuffer; // keep track of atom positions used to generate the pair-list in _Terms
  // If _Matter1 and _Matter2 is defined then _Matter1,_Matter2 are used to build the pair-list - otherwise the _AtomTable is used
  core::T_sp                    _Matter1;
  core::T_sp                    _Matter2;
  core::T_sp                    _NonbondForm; // :AMBER or :ROSETTA form for _Terms
  gctools::Vec0<TermType>	_Terms;
  double                        _RepulsionWeight; // 0.0 - 1.0 weight for ROSETTA NonbondForm
    // Correct way of defining nonbonds using excluded atom indexes
  // FIXME:  Get rid of _FFNonbondDb - and the old code that uses it.
  core::T_sp                      _FFNonbondDb;
  //  Tables for nonbonded calculation using excluded atoms and the Amber way
  size_t                          _ntypes;          // ntypes
  core::SimpleVector_sp           _atom_name_vector;  // atom-name-vector
  NVector_sp                      _charge_vector;          // charge-vector
  NVector_sp                      _mass_vector;            // masses
  core::SimpleVector_int32_t_sp   _atomic_number_vector;    // vec

  // Nonbond parameters generated by generate-nonbond-parameters
  core::SimpleVector_int32_t_sp   _ico_vec;             // ico-vec
  core::SimpleVector_int32_t_sp   _iac_vec;             // iac-vec
  core::SimpleVector_sp      _atom_type_vector;    // Amber atom type names
  core::SimpleVector_int32_t_sp   _local_typej_vec;      // local-typej-vec
  NVector_sp                      _cn1_vec;
  NVector_sp                      _cn2_vec;
  // Excluded atom table
  core::SimpleVector_int32_t_sp   _NumberOfExcludedAtomIndexes;
  core::SimpleVector_int32_t_sp   _ExcludedAtomIndexes;
  size_t                          _InteractionsKept;
  size_t                          _InteractionsDiscarded;

 public:
  virtual std::string implementation_details() const;
  virtual std::string descriptionOfContents() const;
  typedef gctools::Vec0<TermType>::iterator iterator;
  //  iterator begin() { return this->_Terms.begin(); };
  //  iterator end() { return this->_Terms.end(); };
//added by G 7.19.2011
 public:
//  core::List_sp termAtIndex(size_t index) const;
  CL_DEFMETHOD virtual size_t numberOfTerms() { return this->_Terms.size(); };
  CL_DEFMETHOD virtual size_t numberOfTerms14() { return this->_Terms14.size(); };
  void callForEachTerm(core::Function_sp callback);
  void callForEachTerm14(core::Function_sp callback);
 public:

  size_t runTestCalls(core::T_sp stream, chem::NVector_sp coords) const;
  void set_nonbond_pairlist_parameters(double r_switch, double r_cut, double r_pairlist, double distance_dielectric);

  CL_DEFMETHOD core::T_sp keepInteractionFactory() const { return this->_KeepInteractionFactory; };
  CL_DEFMETHOD core::SimpleVector_int32_t_sp number_excluded_atoms() const { return this->_NumberOfExcludedAtomIndexes;}
  CL_DEFMETHOD core::SimpleVector_int32_t_sp excluded_atom_list() const { return this->_ExcludedAtomIndexes;}
 public:
  void constructNonbondTermsBetweenResidues(Residue_sp res1, Residue_sp res2, core::T_sp nbForceField, core::HashTable_sp atomTypes );
  void addTerm14(const TermType& term);
  void addTerm(const TermType& term);
  virtual void dumpTerms(core::HashTable_sp atomTypes);

  virtual core::List_sp extract_vectors_as_alist() const;

  virtual void setupHessianPreconditioner(NVector_sp nvPosition,
                                          AbstractLargeSquareMatrix_sp m,
                                          core::T_sp activeAtomMask );

  void verifyExcludedAtoms(Matter_sp matter, ScoringFunction_sp score);

  virtual double evaluateAllComponent( ScoringFunction_sp scorer,
                                       NVector_sp 	pos,
                                       core::T_sp energyScale,
                                       core::T_sp componentEnergy,
                                       bool 		calcForce,
                                       gc::Nilable<NVector_sp> 	force,
                                       bool		calcDiagonalHessian,
                                       bool		calcOffDiagonalHessian,
                                       gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                       gc::Nilable<NVector_sp>	hdvec,
                                       gc::Nilable<NVector_sp> dvec,
                                       core::T_sp activeAtomMask,
                                       core::T_sp debugInteractions );

  double debugAllComponent( ScoringFunction_sp scorer,
                            NVector_sp 	pos,
                            core::T_sp energyScale,
                            core::T_sp componentEnergy,
                            bool 		calcForce,
                            gc::Nilable<NVector_sp> 	force,
                            bool		calcDiagonalHessian,
                            bool		calcOffDiagonalHessian,
                            gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                            gc::Nilable<NVector_sp>	hdvec,
                            gc::Nilable<NVector_sp> dvec,
                            core::T_sp activeAtomMask,
                            core::T_sp debugInteractions );
  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm( ScoringFunction_sp score,
                                                                                NVector_sp pos,
                                                                                core::T_sp energyScale,
                                                                                core::T_sp activeAtomMask );

  void expandExcludedAtomsToTerms(ScoringFunction_sp score, core::T_sp energyScale );

//    int countBadVdwOverlaps(double scaleSumOfVdwRadii, NVector_sp pos, geom::DisplayList_sp displayIn, core::LispPtr );

  core::T_sp getFFNonbondDb();

  core::T_mv maybeRebuildPairList(core::T_sp tcoordinates);
  core::T_mv rebuildPairList(core::T_sp tcoordinates);
  core::T_mv rebuildPairListBetweenMatters(core::T_sp tcoordinates);
  void constructNonbondTermsFromAtomTable( AtomTable_sp atomTable, core::T_sp nbforceField, core::HashTable_sp atomTypes, core::T_sp keepInteractionFactory, core::T_sp coordinates);
  void constructNonbondTermsBetweenMatters( Matter_sp matter1, Matter_sp matter2, EnergyFunction_sp energyFunction, core::T_sp keepInteractionFactory );
  void construct14InteractionTerms(AtomTable_sp atomTable, Matter_sp matter, core::T_sp nbforceField, core::T_sp keepInteractionFactory, core::HashTable_sp atomTypes );
  void constructExcludedAtomListFromAtomTable(AtomTable_sp atomTable, core::T_sp nbforceField, core::T_sp keepInteractionFactory );

  void constructNonbondTermsFromAtomTableUsingExcludedAtoms(EnergyFunction_sp energyFunction,
                                                            core::T_sp prepareAmberEnergyNonbond,
                                                            core::T_sp keepInteractionFactory );

  void constructNonbondTermsFromAList(core::List_sp values);
  core::List_sp nonbondTermsAsAList();

  void setNonbondExcludedAtomInfo(AtomTable_sp atom_table, core::SimpleVector_int32_t_sp excluded_atoms_list, core::SimpleVector_int32_t_sp number_excluded_atoms);
  void checkEnergyNonbond();
  EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc );

  virtual void emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const;

 public:
  EnergyNonbond_O( const EnergyNonbond_O& ss ); //!< Copy constructor

  EnergyNonbond_O() :
      _UsesExcludedAtoms(false),
      _Nonbond_r_switch(10.0),
      _Nonbond_r_cut(12.0),
      _Nonbond_r_pairlist(14.0),
      _Nonbond_invdd(1.0),
      _KeepInteractionFactory(nil<core::T_O>()),
      _DisplacementBuffer(nil<core::T_O>()),
      _Matter1(nil<core::T_O>()),
      _Matter2(nil<core::T_O>()),
      _FFNonbondDb(nil<core::T_O>()),
      _InteractionsKept(0),
      _InteractionsDiscarded(0)
  {};
};

EnergyFunction_sp energyFunctionNonbondParameters(ScoringFunction_sp score,
                                                  core::T_sp energyScale,
                                                  double& dielectricConstant,
                                                  double& dQ1Q2Scale,
                                                  double& cutoff );

};

#endif //]
