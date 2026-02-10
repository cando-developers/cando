/*
    File: energyAnchorRestraint.cc
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
       
#define	DEBUG_LEVEL_NONE

#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <clasp/core/wrappers.h>

namespace chem {

#define ANCHOR_RESTRAINT_APPLY_ATOM_MASK(I1) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         ) \
    ) goto SKIP_term;
#define ANCHOR_RESTRAINT_APPLY_DEBUG_INTERACTIONS(I1) \
    if (doDebugInteractions) { \
      core::eval::funcall(debugInteractions,EnergyAnchorRestraint_O::static_classSymbol(), \
                          mk_double_float(Energy), \
                          core::make_fixnum(I1); \
    }

EnergyAnchorRestraint::EnergyAnchorRestraint()
{
  this->_Atom1 = nil<Atom_O>();
}

EnergyAnchorRestraint::~EnergyAnchorRestraint()
{
}

core::List_sp EnergyAnchorRestraint::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,ka),core::clasp_make_double_float(this->term.ka)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,xa), core::clasp_make_double_float(this->term.xa)),
                                  core::Cons_O::create(INTERN_(kw,ya), core::clasp_make_double_float(this->term.ya)),
                                  core::Cons_O::create(INTERN_(kw,za), core::clasp_make_double_float(this->term.za)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1));
}


void EnergyAnchorRestraint_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}

#ifdef XML_ARCHIVE
void	EnergyAnchorRestraint::archive(core::ArchiveP node)
{
  node->attribute("ka",this->term.ka);
  node->attribute("xa",this->term.xa);
  node->attribute("ya",this->term.ya);
  node->attribute("za",this->term.za);
  node->attribute("I1",this->term.I1);
  node->attribute("a1",this->_Atom1);
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
  node->attributeIfDefined("calcForce",this->_calcForce,this->_calcForce);
  node->attributeIfDefined("calcDiagonalHessian",this->_calcDiagonalHessian,this->_calcDiagonalHessian);
  node->attributeIfDefined("calcOffDiagonalHessian",this->_calcOffDiagonalHessian,this->_calcOffDiagonalHessian);
#include <cando/chem/energy_functions/_AnchorRestraint_debugEvalSerialize.cc>
#endif //]
}
#endif

#if 0
adapt::QDomNode_sp	EnergyAnchorRestraint::asXml()
{
  adapt::QDomNode_sp	node,child;
  Vector3	vdiff;

  node = adapt::QDomNode_O::create(env,"EnergyAnchorRestraint");
  node->addAttributeString("atomName",this->_Atom1->getName());
  node->addAttributeInt("I1",this->term.I1);
  node->addAttributeDoubleScientific("ka",this->term.ka);
  node->addAttributeDoubleScientific("xa",this->term.xa);
  node->addAttributeDoubleScientific("ya",this->term.ya);
  node->addAttributeDoubleScientific("za",this->term.za);
#if TURN_ENERGY_FUNCTION_DEBUG_ON
  adapt::QDomNode_sp xml = adapt::QDomNode_O::create(env,"Evaluated");
  xml->addAttributeBool("calcForce",this->_calcForce );
  xml->addAttributeBool("calcDiagonalHessian",this->_calcDiagonalHessian );
  xml->addAttributeBool("calcOffDiagonalHessian",this->_calcOffDiagonalHessian );
#include <_AnchorRestraint_debugEvalXml.cc>
  node->addChild(xml);
#endif
  return node;
}

void	EnergyAnchorRestraint::parseFromXmlUsingAtomTable(adapt::QDomNode_sp	xml,
                                                          AtomTable_sp at )
{
  this->term.ka = xml->getAttributeDouble("ka");
  this->term.I1 = xml->getAttributeInt("I1");
  this->term.xa = xml->getAttributeDouble("xa");
  this->term.ya = xml->getAttributeDouble("ya");
  this->term.za = xml->getAttributeDouble("za");
  this->_Atom1 = at->findEnergyAtomWithCoordinateIndex(this->term.I1)->atom();
}
#endif


//
// Copy this from implementAmberFunction.cc
//
double	_evaluateEnergyOnly_AnchorRestraint(
                                                    num_real x1, num_real y1, num_real z1,
                                                    num_real xa, num_real ya, num_real za,
                                                    num_real ka, core::T_sp activeAtomMask )
{
  IMPLEMENT_ME();
  #if 0
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_CALC_FORCE	// Don't calculate FORCE or HESSIAN


#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_AnchorRestraint_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_AnchorRestraint_termCode.cc>

  return Energy;
  #endif
}

std::string EnergyAnchorRestraint_O::descriptionOfContents() const {
  stringstream ss;
  ss << ":enabled " << ((this->_Enabled) ? "T" : "NIL");
  ss << " number-of-terms " << this->_Terms.size();
  return ss.str();
}


void	EnergyAnchorRestraint_O::addTerm(const EnergyAnchorRestraint& r)
{
  this->_Terms.push_back(r);
}

CL_DEFMETHOD
size_t EnergyAnchorRestraint_O::addAnchorRestraintTerm(EnergyFunction_sp energyFunction, Atom_sp a1, Vector3 position, double forceConstant ) {
  AtomTable_sp atomTable = energyFunction->atomTable();
  EnergyAtom* ea = atomTable->getEnergyAtomPointer(a1);
  //EnergyAnchorRestraint restraint(a1,position,ea->coordinateIndexTimes3(),forceConstant);
  this->_Terms.emplace_back(a1,position,ea->coordinateIndexTimes3(),forceConstant);
  return this->_Terms.size()-1;
}



void	EnergyAnchorRestraint_O::dumpTerms(core::HashTable_sp atomTypes)
{
}



void	EnergyAnchorRestraint_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                            AbstractLargeSquareMatrix_sp m,
                                                            core::T_sp activeAtomMask )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  bool		calcForce = true;
  bool		calcDiagonalHessian = true;
  bool		calcOffDiagonalHessian = true;

//
// Copy from implementAmberFunction::setupHessianPreconditioner
//
// -----------------------

#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x=cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
    m->addToElement((i1)+(o1),(i2)+(o2),v);\
  }
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
    m->addToElement((i1)+(o1),(i2)+(o2),v);\
  }
#define ANCHOR_RESTRAINT_CALC_FORCE
#define ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN

  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_AnchorRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    ;
    num_real x1,y1,z1,xa,ya,za,ka;
    int	I1;
    for ( gctools::Vec0<EnergyAnchorRestraint>::iterator cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++ ) {
#include <cando/chem/energy_functions/_AnchorRestraint_termCode.cc>
    }
  }

}





double EnergyAnchorRestraint_O::evaluateAllComponent( ScoringFunction_sp score,
                                                      NVector_sp 	pos,
                                                      core::T_sp energyScale,
                                                      core::T_sp energyComponents,
                                                      bool 		calcForce,
                                                      gc::Nilable<NVector_sp> 	force,
                                                      bool		calcDiagonalHessian,
                                                      bool		calcOffDiagonalHessian,
                                                      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                                      gc::Nilable<NVector_sp>	hdvec,
                                                      gc::Nilable<NVector_sp> dvec,
                                                      core::T_sp activeAtomMask,
                                                      core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  double totalEnergy = 0.0;
  this->_Evaluations++;
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());
//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------

#define ANCHOR_RESTRAINT_CALC_FORCE
#define ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) totalEnergy += (e);
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE 		ForceAcc
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc


  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_AnchorRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,xa,ya,za,ka;
    int	I1, i;
    gctools::Vec0<EnergyAnchorRestraint>::iterator cri;
    for ( i=0,cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++,i++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfAnchorRestraintTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfAnchorRestraintTermsToCalculate ) {
          break;
        }
      }
#endif

			    /* Obtain all the parameters necessary to calculate */
			    /* the amber and forces */
#include <cando/chem/energy_functions/_AnchorRestraint_termCode.cc>

#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      cri->_calcForce = calcForce;
      cri->_calcDiagonalHessian = calcDiagonalHessian;
      cri ->_calcOffDiagonalHessian = calcOffDiagonalHessian;
			// Now write all of the calculated values into the eval structure
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ cri->eval.var=val;};
#include <cando/chem/energy_functions/_AnchorRestraint_debugEvalSet.cc>
#endif //]

			    /* Add the forces */

      if ( calcForce ) {
#if !USE_EXPLICIT_DECLARES
        num_real fx1 = 0.0;
        num_real fy1 = 0.0;
        num_real fz1 = 0.0;
#endif
//		    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		    _lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
      }
    }
  }
  maybeSetEnergy( energyComponents, EnergyAnchorRestraint_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}







void	EnergyAnchorRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(NVector_sp 	pos,
                                                                                        core::T_sp activeAtomMask )
{
#if 0
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  core::T_sp debugInteractions = nil<core::T_O>();
  MAYBE_SETUP_DEBUG_INTERACTIONS(false);
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;

//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------

#define ANCHOR_RESTRAINT_CALC_FORCE
#define ANCHOR_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define ANCHOR_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	ANCHOR_RESTRAINT_SET_PARAMETER
#define	ANCHOR_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	ANCHOR_RESTRAINT_SET_POSITION
#define	ANCHOR_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	ANCHOR_RESTRAINT_PHI_SET
#define	ANCHOR_RESTRAINT_PHI_SET(x) {}
#undef	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE
#define	ANCHOR_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	ANCHOR_RESTRAINT_FORCE_ACCUMULATE
#define	ANCHOR_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	ANCHOR_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


  {
   
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_AnchorRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,xa,ya,za,ka;
    int	I1, i;
    gctools::Vec0<EnergyAnchorRestraint>::iterator cri;
    for ( i=0,cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_AnchorRestraint_termCode.cc>
      LOG("fx1 = {}" , fx1 );
      LOG("fy1 = {}" , fy1 );
      LOG("fz1 = {}" , fz1 );
      int index = i;
#if !USE_EXPLICIT_DECLARES
      num_real fx1 = 0.0;
      num_real fy1 = 0.0;
      num_real fz1 = 0.0;
      num_real dhx1x1 = 0.0;
      num_real ohx1y1 = 0.0;
      num_real dhy1y1 = 0.0;
      num_real ohy1z1 = 0.0;
      num_real dhz1z1 = 0.0;
      num_real ohx1z1 = 0.0;
#endif
#include <cando/chem/energy_functions/_AnchorRestraint_debugFiniteDifference.cc>

    }
  }
  #endif
  IMPLEMENT_ME(); // must return some sort of integer value
}



void EnergyAnchorRestraint_O::initialize()
{
  this->Base::initialize();
}

EnergyComponent_sp EnergyAnchorRestraint_O::copyFilter(core::T_sp keepInteractionFactory, SetupAccumulator& setupAcc) {
  EnergyAnchorRestraint_sp copy = EnergyAnchorRestraint_O::create();
  copyEnergyComponent( copy, this->asSmartPtr() );
  core::T_sp keepInteraction = specializeKeepInteractionFactory(keepInteractionFactory,EnergyAnchorRestraint_O::staticClass());
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    if ( skipInteraction( keepInteraction, a1 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}

core::T_sp test_allocate_EnergyAnchorRestraint_O() {
#if 1
  // succeeds
  auto thing = EnergyAnchorRestraint_O::create();
#elif 1
  // succeeds ---- I got it to fail once I think but it may have been that I didn't save 
  auto thing = gctools::GC<EnergyAnchorRestraint_O>::allocate();
#elif 0
  // FAILS
  auto kind = gctools::Header_s::StampWtagMtag::make_StampWtagMtag(EnergyAnchorRestraint_O::static_ValueStampWtagMtag);
  auto thing = gctools::GCObjectAllocator<EnergyAnchorRestraint_O>::template allocate_kind<gctools::RuntimeStage>(kind,gctools::sizeof_with_header<EnergyAnchorRestraint_O>());
#else
  // FAILS
  gctools::Header_s::BadgeStampWtagMtag the_header;
  auto thing = gctools::GCObjectAppropriatePoolAllocator<EnergyAnchorRestraint_O,gctools::normal>::allocate_in_appropriate_pool_kind<gctools::RuntimeStage>(the_header,9999);
#endif
  return thing;
}

};
