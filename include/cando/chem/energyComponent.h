/*
    File: energyComponent.h
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
#define	DEBUG_LEVEL_FULL

//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	energyComponent.h
 *
 *	Maintains a database of stretch types
 */

#ifndef EnergyComponent_H 
#define	EnergyComponent_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/quickDom.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/nVector.h>
#include <cando/chem/scoringFunction.h>
#include <cando/chem/residue.h>


#define USE_MATHEMATICA_KERNELS 0

#define SIZE_T size_t

		//
		// SCALE FACTOR FOR ELECTROSTIC INTERACTIONS???
		//
#define	DECLARE_FLOAT(x) num_real x = 0.0

	// If we declare mathematica generated term variables
	// on the fly then set ON_THE_FLY_DECLARES
	// and variable declarations won't be included
	//
	// This should avoid a lot of variables being declared but
	// not used and eliminate a lot of compiler warnings
	//
#define	USE_EXPLICIT_DECLARES 1



#ifndef	delta2
#define		delta2	0.00001
#endif
#include <cando/adapt/quickDom.fwd.h>// energyComponent.h wants QDomNode needs quickDom.fwd.h

namespace chem {


typedef enum { energyEval, gradientEval, hessianEval } EnergyComponentEvalType;

inline EnergyComponentEvalType determineEnergyComponentEvalType(gc::Nilable<NVector_sp> force, gc::Nilable<NVector_sp> hdvec, gc::Nilable<NVector_sp> dvec ) {
  if (hdvec.notnilp()) {
    if (force.nilp()) {
      SIMPLE_ERROR("You cannot have hdvec not nil but force nil");
    } else if (dvec.nilp()) {
      SIMPLE_ERROR("You cannot have hdvec not nil but dvec nil");
    }
    return hessianEval;
  }
  if (force.notnilp()) return gradientEval;
  return energyEval;
}


class KahanSummation {
private:
  double sum;
  double c; // A running compensation for lost low-order bits.

public:
      // Constructor initializes sum and compensation to zero.
  KahanSummation() : sum(0.0), c(0.0) {}

      // Adds a new value to the running total using Kahan summation algorithm.
  void add(double input) {
    double y = input - c;
    double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }

      // Returns the current value of the sum.
  double getSum() const {
    return sum;
  }
};
};
namespace       chem {


#define MAYBE_SETUP_ACTIVE_ATOM_MASK() \
  bool hasActiveAtomMask = (activeAtomMask.notnilp()); \
  auto bitvectorActiveAtomMask = gc::As_unsafe<core::SimpleBitVector_sp>(activeAtomMask);
#define MAYBE_SETUP_DEBUG_INTERACTIONS(dbgint) bool doDebugInteractions = dbgint;

// Active-atom masks are Cartesian-coordinate masks: three bits per atom, at
// the same i3x offsets used by energy terms. A term contributes whenever at
// least one coordinate of at least one participating atom is movable.
inline bool activeAtomMaskAtomIsActive(core::SimpleBitVector_sp mask, size_t i3x) {
  return mask->testBit(i3x) || mask->testBit(i3x + 1) || mask->testBit(i3x + 2);
}

inline bool activeAtomMaskAnyAtomIsActive(core::SimpleBitVector_sp mask, size_t i3x1) {
  return activeAtomMaskAtomIsActive(mask, i3x1);
}

inline bool activeAtomMaskAnyAtomIsActive(core::SimpleBitVector_sp mask,
                                          size_t i3x1, size_t i3x2) {
  return activeAtomMaskAtomIsActive(mask, i3x1) ||
         activeAtomMaskAtomIsActive(mask, i3x2);
}

inline bool activeAtomMaskAnyAtomIsActive(core::SimpleBitVector_sp mask,
                                          size_t i3x1, size_t i3x2, size_t i3x3) {
  return activeAtomMaskAtomIsActive(mask, i3x1) ||
         activeAtomMaskAtomIsActive(mask, i3x2) ||
         activeAtomMaskAtomIsActive(mask, i3x3);
}

inline bool activeAtomMaskAnyAtomIsActive(core::SimpleBitVector_sp mask,
                                          size_t i3x1, size_t i3x2,
                                          size_t i3x3, size_t i3x4) {
  return activeAtomMaskAtomIsActive(mask, i3x1) ||
         activeAtomMaskAtomIsActive(mask, i3x2) ||
         activeAtomMaskAtomIsActive(mask, i3x3) ||
         activeAtomMaskAtomIsActive(mask, i3x4);
}

SMART(QDomNode);
SMART(AbstractLargeSquareMatrix);
FORWARD(EnergyFunction);

#if TURN_ENERGY_FUNCTION_DEBUG_ON == 1
#define	LOG_ENERGY_CLEAR() {this->_DebugLog.str("");};
#define LOG_ENERGY(...) {fmt::print("{}:{}:{} ", __FILE__, __LINE__, __FUNCTION__ ); fmt::print(__VA_ARGS__);}
#else
#define	LOG_ENERGY_CLEAR() {}
#define LOG_ENERGY(...) {}
#endif


class EnergyTerm 
{
 private:
 public:
  virtual string	className() {THROW_HARD_ERROR("Subclass must implement");};

  EnergyTerm() {};
};


inline string	atomLabel(Atom_sp a)
{
  stringstream	sstr;
  sstr << _rep_(a->getName());
  return sstr.str();
}

#ifndef VERYSMALL
#define VERYSMALL       0.000000000001
#define verysmall       VERYSMALL
#endif
#define VERYSMALLSQUARED       (VERYSMALL*VERYSMALL)
#define	TENM3		0.001


#ifndef	SIGN
#define	SIGN(x)	((x)<0.0?-1.0:1.0)
#endif
#define	myPow(x,y)	(pow(x,y))
#define	ArcCos(x)	(acos(x))
#define	ArcTan(x,y)	(atan2(y,x))
#define myatan2(y,x)    (atan2(y,x))
#define	Cos(x)		(cos(x))
#define	Sin(x)		(sin(x))

#define	Pi		(3.14159265)
#define	TWOPI		(3.14159265*2.0)
#define	reciprocal(x)	(1.0/(x))
#define	mysqrt(x)	sqrt(x)
#define	mymod(x,y)	fmod(x,y)
#define	power1(x)	(x)
#define	power2(x)	((x)*(x))
#define	powern2(x)	(1.0/(x)*(x))
#define	c_two		2.0
#define	c_six		6.0
#define	oneHalf		0.5
#define	c_twelve	12.0
#define	c_True		true
#define	c_False		false

	// ChiralOffset is a term added in the chiral restraint energy function
	// to try and prevent the four atoms in a chiral restraint from
	// settling into the same plane
#define	ChiralOffset	0.0



#define	TEST_FORCE( func, delta, argLow, argHigh, term, idx ) {	\
    num_real eLow = func argLow;\
    num_real eHigh = func argHigh;\
    num_real numForce = -(eHigh-eLow)/(delta);\
    LOG(("eHigh = %le eLow = %le delta = %le") , (eHigh) , (eLow) , (delta) ); \
    LOG(("numForce = %le") , numForce );			\
    if ( !_areValuesClose( numForce, term, #func, #term, idx ) ) { fails++; };\
  }

#define	TEST_DIAGONAL_HESSIAN( func, delta, argLow, argMiddle, argHigh, term, idx ) {\
    num_real eLow = func argLow;\
    num_real eMiddle = func argMiddle;\
    num_real eHigh = func argHigh;\
    num_real numHessian = (eHigh+eLow-2.0*eMiddle)/(((delta)/2.0)*((delta)/2.0));\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }

#define	TEST_OFF_DIAGONAL_HESSIAN( func, delta, argMM, argPM, argMP, argPP, term, idx ) {\
    num_real eMM = func argMM;\
    num_real eMP = func argMP;\
    num_real ePM = func argPM;\
    num_real ePP = func argPP;\
    num_real numHessian = ((ePP-ePM)-(eMP-eMM))/(delta*delta);\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }


// ----------------------------------------------------------------------
// New kernel function accumulators
// ----------------------------------------------------------------------


#define KernelGradientAcc(ii1,oo1,vv1) force[ii1+oo1] += (-vv1)

//
// Accumulate an off diagonal Hessian element
//
#define	KernelOffDiagHessAcc(i1,o1,i2,o2,v) {\
    auto v22 = v*dvec[i2+o2];\
    auto v11 = v*dvec[i1+o1];\
    hdvec[i1+o1] += v22; \
    hdvec[i2+o2] += v11; \
  }

//
// Accumulate a diagonal Hessian element
//
#define	KernelDiagHessAcc(i1,o1,i2,o2,v) {\
      auto vd = v*dvec[i1+o1];\
      hdvec[i1+o1] += vd; \
  }



struct NoHessian {};

template <typename HESSIAN=NoHessian>
inline void KernelHessDiagAcc(size_t positionSize, HESSIAN hessian, double* dvec, double* hdvec, size_t i3x1, size_t of1, size_t i3x2, size_t of2, double hval) {
  KernelDiagHessAcc( i3x1, of1, i3x2, of2, hval );
}

template <typename HESSIAN=NoHessian>
inline void KernelHessOffDiagAcc(size_t positionSize, HESSIAN hessian, double* dvec, double* hdvec, size_t i3x1, size_t of1, size_t i3x2, size_t of2, double hval) {
  KernelOffDiagHessAcc( i3x1, of1, i3x2, of2, hval );
}

template <>
inline void KernelHessDiagAcc(size_t positionSize, double* hessian, double* dvec, double* hdvec, size_t i3x1, size_t of1, size_t i3x2, size_t of2, double hval) {
  hessian[(i3x1+of1)*positionSize+(i3x2+of2)] += hval;
  KernelDiagHessAcc( i3x1, of1, i3x2, of2, hval );
}

template <>
inline void KernelHessOffDiagAcc(size_t positionSize, double* hessian, double* dvec, double* hdvec, size_t i3x1, size_t of1, size_t i3x2, size_t of2, double hval) {
  hessian[(i3x1+of1)*positionSize+(i3x2+of2)] += hval;
  hessian[(i3x2+of2)*positionSize+(i3x1+of1)] += hval;
  KernelOffDiagHessAcc( i3x1, of1, i3x2, of2, hval );
}


// ----------------------------------------------------------------------
// Old accumulators
// ----------------------------------------------------------------------

#define	ForceAcc(i,o,v) {\
      if ( hasForce ) {\
        force->setElement((i)+(o),(v)+force->getElement((i)+(o)));\
      }\
  }

//
// Accumulate an off diagonal Hessian element
//
#define	OffDiagHessAcc(i1,o1,i2,o2,v) {\
    if ( hasHessian ) {\
      hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
    }\
    if ( hasHdAndD ) {\
      auto v22 = v*dvec->element((i2)+(o2));\
      auto v11 = v*dvec->element((i1)+(o1));\
      hdvec->addToElement((i1)+(o1),v22); \
      hdvec->addToElement((i2)+(o2),v11); \
    }\
  }

//
// Accumulate a diagonal Hessian element
//
#define	DiagHessAcc(i1,o1,i2,o2,v) {\
    if ( hasHessian ) {\
      hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
    }\
    if ( hasHdAndD ) {\
      auto vd = v*dvec->element((i1)+(o1));\
      hdvec->addToElement((i1)+(o1),vd); \
    }\
  }



#define	ALL_ENERGY_COMPONENTS(msg) {		\
    this->_Stretch->msg; 			\
      hdvec->addToElement((i1)+(o1),v*dvec->element((i2)+(o2)));\
      hdvec->addToElement((i2)+(o2),v*dvec->element((i1)+(o1)));\
    }\
  }

#if 0
// This is a duplicate of the one above - why is it here?  Delete it
//
// Accumulate a diagonal Hessian element
//
#define	DiagHessAcc(i1,o1,i2,o2,v) {\
      if ( hasHessian ) {\
          hessian->addToElement((i1)+(o1),(i2)+(o2),v);\
      }\
      if ( hasHdAndD ) {\
          hdvec->addToElement((i1)+(o1),v*dvec->element((i1)+(o1)));\
      }\
  }
#endif


#define	ALL_ENERGY_COMPONENTS(msg) {		\
    this->_Stretch->msg; 			\
    this->_Angle->msg; 			\
    this->_Dihedral->msg; 			\
    this->_Nonbond->msg; 			\
    this->_DihedralRestraint->msg; 		\
    this->_ChiralRestraint->msg; 		\
    this->_AnchorRestraint->msg; 		\
    this->_FixedNonbondRestraint->msg; 	\
  }


SMART(EnergyComponent );
SMART(EnergyPairlistComponent );
class EnergyComponent_O : public core::CxxObject_O
{
  LISP_CLASS(chem,ChemPkg,EnergyComponent_O,"EnergyComponent",core::CxxObject_O);
public:
public: // virtual functions inherited from Object
  void initialize();
public: // instance variables
  bool		        _Enabled;
  double		_Scale;
public:
  size_t                _Evaluations;
public:
  virtual bool restraintp() const = 0;
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
  CL_DEFMETHOD virtual std::string implementation_details() const { return "energyComponent.h: no-details-available-subclass-me"; };
  CL_DEFMETHOD virtual size_t numberOfTerms() {_OF(); SUBCLASS_MUST_IMPLEMENT();};
  void setScale(double s) {this->_Scale = s; };
  double getScale() { return this->_Scale ; };
  CL_DEFMETHOD bool isEnabled() const { return this->_Enabled; };
  CL_LISPIFY_NAME("enable");
  CL_DEFMETHOD 	void enable() {this->_Enabled = true; };
  CL_LISPIFY_NAME("disable");
  CL_DEFMETHOD 	void disable() {this->_Enabled = false; };

  CL_DOCSTRING("Some energy-components are restraints - and should be disabled to calculate conformational energies");
  CL_DEFMETHOD
  virtual bool is_restraint() const { return false; };
  
  CL_LISPIFY_NAME("energy-component-evaluations");
  CL_DEFMETHOD 	size_t	evaluations() const { return this->_Evaluations; };

  virtual void atomsForEachTerm(core::Function_sp callback);

  CL_DEFMETHOD virtual void emitTestCalls(core::T_sp stream, chem::NVector_sp pos) const {SUBCLASS_MUST_IMPLEMENT(); };
  CL_DEFMETHOD virtual size_t runTestCalls(core::T_sp stream, chem::NVector_sp pos) const {SUBCLASS_MUST_IMPLEMENT(); };
  
  CL_DEFMETHOD virtual core::List_sp extract_vectors_as_alist() const { SUBCLASS_MUST_IMPLEMENT(); };
 
  string enabledAsString();
  string debugLogAsString();

public:	// Virtual methods

  CL_DEFMETHOD virtual void dumpTerms(core::HashTable_sp atomTypes) {_OF();SUBCLASS_MUST_IMPLEMENT();};

  virtual EnergyComponent_sp copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc ) {SUBCLASS_MUST_IMPLEMENT();};
  
  CL_DEFMETHOD virtual	double evaluateAllComponent( ScoringFunction_sp scorer,
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
                                                     core::T_sp debugInteractions ) = 0;


  virtual	void	compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp pos ) {_OF();SUBCLASS_MUST_IMPLEMENT();};
public:
  EnergyComponent_O( const EnergyComponent_O& ss ); //!< Copy constructor

  EnergyComponent_O() : _Enabled(true), _Scale(1.0) {};

  virtual void setupHessianPreconditioner(NVector_sp nvPosition, AbstractLargeSquareMatrix_sp m, core::T_sp activeAtomMask) = 0;

};


/*! Base for every component that builds a PAIR LIST - EnergyNonbond_O, EnergyRosettaNonbond_O,
    EnergyRosettaElec_O and EnergyRosettaLKSolvation_O.

    These five members were declared identically in all four classes.  rebuildPairListImpl
    (pairList.h) is a template that reaches them by duck typing, with the requirements written out
    as a comment at the top of that file; collecting them here makes that contract a declared base
    class instead, so a new pair-list component inherits the interface rather than reproducing it.

    Nothing here is virtual.  The template still binds by static type, so there is no dispatch cost
    and the derived classes stay layout-compatible with what they were. */
class EnergyPairlistComponent_O : public EnergyComponent_O
{
  LISP_ABSTRACT_CLASS(chem,ChemPkg,EnergyPairlistComponent_O,"EnergyPairlistComponent",EnergyComponent_O);
public:
  core::T_sp            _KeepInteractionFactory;
  /*! When the current pair list was built, and how far the coordinates had already drifted then.
   *
   * These two scalars replace _DisplacementBuffer, which was a full COPY of the coordinates - 162 KB
   * per component - walked in its entirety on every evaluation.  A blueprint holds ~2829 pair-list
   * components over 6762 atoms, so that was ~460 MB of copies and ~57 million coordinate comparisons
   * per cycle, all of them recomputing the same answer from the same numbers.
   *
   * The snapshot now lives once on the AtomTable (AtomTable_O::_RefCoords) and each component keeps
   * only its position relative to it: the epoch the reference was in, and the drift it had reached.
   * maybeRebuildPairListImpl turns those into the same decision with two loads and a compare - see
   * pairList.h and the _RefCoords comment for the triangle-inequality bound that makes it sound.
   *
   * EPOCH 0 MEANS NEVER BUILT, which is why AtomTable_O::_RefEpoch starts at 1 and skips 0 on wrap.
   * Everything that used to nil out _DisplacementBuffer now zeroes the epoch, and it means the same
   * thing: the next evaluation must rebuild rather than trust what is cached. */
  size_t                _PairListEpoch = 0;
  double                _PairListDrift = 0.0;
  //! False when no bonded path can span the two matters - two molecules, or a sidechain against an
  //! intermolecular partner - so the 1-2/1-3/1-4 probes would always miss and are pure cost.
  bool                  _ExclusionsPossible = true;
  //! When both are set the pair list is built BETWEEN them rather than over the whole atom table.
  core::T_sp            _Matter1;
  core::T_sp            _Matter2;

  /*! QUERY SET - which atoms this component enumerates candidate pairs FROM.
   *
   * Unset (the default) means "every atom": the full-table grid enumeration rebuildPairListImpl
   * has always done.  Every existing caller leaves it unset and is unaffected.
   *
   * Set, it means: iterate THESE atoms, read each one's neighbour list off the AtomTable, and form
   * candidates from that.  Scoring one rotamer the old way examined ~6.5 million candidate pairs -
   * gridding all 6762 atoms to find the ~4500 involving 15 of them - and did it once per component
   * per evaluation, with the coordinates unchanged throughout.
   *
   * Two forms, because the two callers have different shapes.  A fan-out slot is a CONTIGUOUS i3
   * range, so it is two integers and no allocation.  A backbone is not contiguous and needs the
   * mask, which its caller already holds.  Forcing both through a mask would spend 470 masks of
   * 2.5 KB to say what two fixnums say.
   *
   * INVARIANT: the query set must be a SUPERSET of the atoms keepInteraction would accept on the
   * query side.  It is an enumeration hint, not a second filter - the predicate still runs and
   * still decides.  Too wide costs only speed; TOO NARROW SILENTLY DROPS TERMS.
   */
  size_t                _QueryLo3 = 0;
  size_t                _QueryHi3 = 0;      //!< 0 means no range set
  core::T_sp            _QueryMask;         //!< NIL means no mask set; else an i3-indexed bit-vector
  /*! Which atoms may be FOUND as neighbours - the other end of the query set.
   *
   * The query set says what to enumerate FROM; this says what may be enumerated TO, and it is what
   * the AtomTable's neighbour list gets built against.  For a blueprint it is the backbone, so a
   * slot finds backbone partners and never another rotamer's atoms.
   *
   * NIL means every atom is a target.  Only consulted when a query set is set. */
  core::T_sp            _QueryTargetMask;
  /*! Cutoff to BUILD the shared neighbour list at, when it differs from this component's own.
   *
   * 0 means "use my own rpairlist".  Set it to the MAX rpairlist across the group's components and
   * every one of them builds the same list: validity is `_NeighborCutoff >= requested`, so a wider
   * list answers a narrower request, and each component still applies its own exact distance test
   * per pair.
   *
   * Without it the FIRST component to evaluate sets the width.  Today that is rosetta nonbond at
   * 9.0 and the others reuse it, but only by the order components happen to be added - override
   * :rpairlist to make nonbond narrower than LK and the first build would be too small, forcing a
   * second full pass.  PROTEIN-BINDER.LISP:2328 takes the same max for the same reason. */
  double                _QueryCutoff = 0.0;

public:
  void fields(core::Record_sp node);

public:
  CL_DEFMETHOD core::T_sp keepInteractionFactory() const { return this->_KeepInteractionFactory; };
  core::T_sp matter1() const { return this->_Matter1; };
  core::T_sp matter2() const { return this->_Matter2; };
  bool exclusionsPossible() const { return this->_ExclusionsPossible; };
  size_t pairListEpoch() const { return this->_PairListEpoch; };
  double pairListDrift() const { return this->_PairListDrift; };
  /*! Record that the pair list was just built against EPOCH with the reference at DRIFT. */
  void notePairListBuilt(size_t epoch, double drift) {
    this->_PairListEpoch = epoch;
    this->_PairListDrift = drift;
  };
  CL_LISPIFY_NAME("invalidate-pair-list");
  CL_DOCSTRING(R"dx(Force the next evaluation to rebuild this component's pair list.)dx");
  CL_DEFMETHOD void invalidatePairList() { this->_PairListEpoch = 0; };

  bool hasQuerySet() const {
    return this->_QueryHi3 != 0 || this->_QueryMask.notnilp();
  };
  size_t queryLo3() const { return this->_QueryLo3; };
  size_t queryHi3() const { return this->_QueryHi3; };
  core::T_sp queryMask() const { return this->_QueryMask; };

  CL_LISPIFY_NAME("set-query-range");
  CL_DOCSTRING(R"dx(Enumerate candidate pairs from the atoms with i3 from LO3 up to but not including HI3.)dx");
  CL_DEFMETHOD void setQueryRange(size_t lo3, size_t hi3) {
    if (hi3 <= lo3)
      SIMPLE_ERROR("Empty query range lo3={} hi3={} - use clear-query-set for 'every atom'", lo3, hi3);
    this->_QueryLo3 = lo3;
    this->_QueryHi3 = hi3;
    this->invalidatePairList();
  };
  CL_LISPIFY_NAME("set-query-mask");
  CL_DOCSTRING(R"dx(Enumerate candidate pairs from the atoms whose i3 bit is set in MASK.)dx");
  CL_DEFMETHOD void setQueryMask(core::T_sp mask) {
    this->_QueryMask = mask;
    this->invalidatePairList();
  };
  /*! Build the pair list BETWEEN matter1 and matter2 instead of over the whole atom table.
   *
   * EXCLUSIONS-POSSIBLE NIL says no bonded path can span the two matters - two molecules, or a
   * sidechain against an intermolecular partner - so the 1-2/1-3/1-4 probes would always miss and
   * are pure cost per pair.
   *
   * Invalidating the pair list is not incidental: maybeRebuildPairList decides whether to rebuild
   * from how far the coordinates have drifted, and the coordinates do not move when the MATTERS
   * change - so leaving the stamp in place would let the next evaluation happily reuse a pair list
   * built for the PREVIOUS pair of matters.
   *
   * The three rosetta components each had this verbatim and byte-identical, including a
   * commented-out invalidateParameterCache() call.  Leaving that out is correct: their parameter
   * caches are keyed on the atom table and the force field, neither of which setMatters touches. */
  CL_LAMBDA((self chem:energy-pairlist-component) matter1 matter2 &optional (exclusions-possible t));
  CL_DEFMETHOD void setMatters(core::T_sp matter1, core::T_sp matter2, bool exclusionsPossible = true) {
    this->_Matter1 = matter1;
    this->_Matter2 = matter2;
    this->_ExclusionsPossible = exclusionsPossible;
    this->invalidatePairList();
  };
  core::T_sp queryTargetMask() const { return this->_QueryTargetMask; };
  double queryCutoff() const { return this->_QueryCutoff; };
  CL_LISPIFY_NAME("set-query-cutoff");
  CL_DOCSTRING(R"dx(Build the shared neighbour list at CUTOFF rather than at this component's own rpairlist.  Pass the max across the group so every component shares one list.  Zero reverts to using rpairlist.)dx");
  CL_DEFMETHOD void setQueryCutoff(double cutoff) { this->_QueryCutoff = cutoff; };
  CL_LISPIFY_NAME("set-query-targets");
  CL_DOCSTRING(R"dx(Restrict the atoms that may be FOUND as neighbours to those whose i3 bit is set in MASK.  NIL means every atom.)dx");
  CL_DEFMETHOD void setQueryTargets(core::T_sp mask) {
    this->_QueryTargetMask = mask;
    this->invalidatePairList();
  };
  CL_LISPIFY_NAME("clear-query-set");
  CL_DOCSTRING(R"dx(Revert to enumerating candidate pairs from every atom.)dx");
  CL_DEFMETHOD void clearQuerySet() {
    this->_QueryLo3 = 0;
    this->_QueryHi3 = 0;
    this->_QueryMask = nil<core::T_O>();
    this->_QueryTargetMask = nil<core::T_O>();
    this->invalidatePairList();
  };

public:
  EnergyPairlistComponent_O()
    : _KeepInteractionFactory(nil<core::T_O>())
    , _Matter1(nil<core::T_O>())
    , _Matter2(nil<core::T_O>())
    , _QueryMask(nil<core::T_O>())
    , _QueryTargetMask(nil<core::T_O>())
  {};
};

template <typename SP>
SP safe_alist_lookup(core::List_sp list, core::T_sp key) {
  if (list.consp()) {
    core::T_sp entry = list.unsafe_cons()->assoc(key,nil<core::T_O>(),nil<core::T_O>(),nil<core::T_O>());
    if (entry.consp()) {
      return gctools::As<SP>(CONS_CDR(entry));
    }
  }
  SIMPLE_ERROR("Could not find {}" , _rep_(key));
}

void copyEnergyComponent(EnergyComponent_sp newComponent, EnergyComponent_sp originalComponent );
;
};

namespace chem {

void test_zero( size_t num,
                double* force_new, double* force_ground,
                double* hessian_new, double* hessian_ground,
                double* dvec_new, double* dvec_ground,
                double* hdvec_new, double* hdvec_ground );

bool test_match( core::T_sp stream, const char* label, size_t num,
                double* force_new, double* force_ground,
                double* hessian_new, double* hessian_ground,
                double* hdvec_new, double* hdvec_ground );

void test_position(core::T_sp stream, size_t pos_size, double* position );

};

#endif
