/*
    File: energyChiralRestraint.cc
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


#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/bond.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/largeSquareMatrix.h>
#include <clasp/core/wrappers.h>



//#define	CHIRAL_RESTRAINT_WEIGHT 1000000.0	// 10000.0 is used in MOE


namespace chem {


#define CHIRAL_RESTRAINT_APPLY_ATOM_MASK(I1,I2,I3,I4) \
if (hasActiveAtomMask \
    && !(bitvectorActiveAtomMask->testBit(I1/3) \
         && bitvectorActiveAtomMask->testBit(I2/3) \
         && bitvectorActiveAtomMask->testBit(I3/3) \
         && bitvectorActiveAtomMask->testBit(I4/3) \
         ) \
    ) goto SKIP_term;

string	EnergyChiralRestraint::description()
{
    stringstream ss;
    ss << "EnergyChiralRestraint[";
    ss << "atoms: ";
    ss << this->_Atom1->description() << "-";
    ss << this->_Atom2->description() << "-";
    ss << this->_Atom3->description() << "-";
    ss << this->_Atom4->description() ;
#if TURN_ENERGY_FUNCTION_DEBUG_ON
    ss << " eval.Energy=" << this->eval.Energy;
#endif
    return ss.str();
}





//
// Copy this from implementAmberFunction.cc
//
num_real	_evaluateEnergyOnly_ChiralRestraint(
		num_real x1, num_real y1, num_real z1,
		num_real x2, num_real y2, num_real z2,
		num_real x3, num_real y3, num_real z3,
		num_real x4, num_real y4, num_real z4,
		num_real K, num_real CO )
{
  IMPLEMENT_ME();
  #if 0
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_CALC_FORCE	// Don't calculate FORCE or HESSIAN

#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>

    return Energy;
    #endif
}

SYMBOL_EXPORT_SC_(KeywordPkg,central_atom3);

core::List_sp EnergyChiralRestraint::encode() const {
  return core::Cons_O::createList(core::Cons_O::create(INTERN_(kw,k),core::clasp_make_double_float(this->term.K)),
                                  core::Cons_O::create(INTERN_(kw,CO),core::clasp_make_double_float(this->term.CO)),
                                  core::Cons_O::create(INTERN_(kw,i1), core::make_fixnum(this->term.I1)),
                                  core::Cons_O::create(INTERN_(kw,i2), core::make_fixnum(this->term.I2)),
                                  core::Cons_O::create(INTERN_(kw,i3), core::make_fixnum(this->term.I3)),
                                  core::Cons_O::create(INTERN_(kw,i4), core::make_fixnum(this->term.I4)),
                                  core::Cons_O::create(INTERN_(kw,atom1), this->_Atom1),
                                  core::Cons_O::create(INTERN_(kw,atom2), this->_Atom2),
                                  core::Cons_O::create(INTERN_(kw,central_atom3), this->_Atom3),
                                  core::Cons_O::create(INTERN_(kw,atom4), this->_Atom4));
}

EnergyChiralRestraint_sp EnergyChiralRestraint_O::copyFilter(core::T_sp keepInteraction) {
  EnergyChiralRestraint_sp copy = EnergyChiralRestraint_O::create();
  for ( auto edi=this->_Terms.begin(); edi!=this->_Terms.end(); edi++ ) {
    Atom_sp a1 = edi->_Atom1;
    Atom_sp a2 = edi->_Atom2;
    Atom_sp a3 = edi->_Atom3;
    Atom_sp a4 = edi->_Atom4;
    if ( skipInteraction( keepInteraction, EnergyChiralRestraint_O::staticClass(), a1, a2, a3, a4 ) ) continue;
    copy->_Terms.push_back(*edi);
  }
  return copy;
}


void EnergyChiralRestraint_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,terms), this->_Terms );
  this->Base::fields(node);
}


void EnergyChiralRestraint_O::addTerm(const EnergyChiralRestraint& e)
{
    this->_Terms.push_back(e);
}

void EnergyChiralRestraint_O::dumpTerms(core::HashTable_sp atomTypes)
{
    gctools::Vec0<EnergyChiralRestraint>::iterator cri;
    int i;
    string as1, as2, as3, as4;
    for ( i=0,cri=this->_Terms.begin();
	  cri!=this->_Terms.end(); cri++,i++ ) {
	as1 = atomLabel(cri->_Atom1);
	as2 = atomLabel(cri->_Atom2);
	as3 = atomLabel(cri->_Atom3);
	as4 = atomLabel(cri->_Atom4);
        core::clasp_write_string(fmt::format("TERM 7CHIRAL {:<9} {:<9} {:<9} {:<9} {:8.2f} {:8.2f} ; I1={} I2={} I3={} I4={}\n" , as1 , as2 , as3 , as4 , cri->term.K , cri->term.CO , cri->term.I1 , cri->term.I2 , cri->term.I3 , cri->term.I4 ));
    }
}


string EnergyChiralRestraint_O::beyondThresholdInteractionsAsString()
{
    return component_beyondThresholdInteractionsAsString<EnergyChiralRestraint_O,EnergyChiralRestraint>(*this);
}




void	EnergyChiralRestraint_O::setupHessianPreconditioner(chem::NVector_sp nvPosition,
                                                            chem::AbstractLargeSquareMatrix_sp m,
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


#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x=cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of) {x=nvPosition->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {\
	  m->addToElement((i1)+(o1),(i2)+(o2),v);\
}
#define CHIRAL_RESTRAINT_CALC_FORCE
#define CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN

 {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include	<cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	double x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
	int	I1, I2, I3, I4;
	for ( gctools::Vec0<EnergyChiralRestraint>::iterator cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++ ) {
#include	<cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
	}
    }
}





num_real EnergyChiralRestraint_O::evaluateAllComponent( ScoringFunction_sp score,
                                                        chem::NVector_sp 	pos,
                                                        core::T_sp componentEnergy,
                                                        bool 		calcForce,
                                                        gc::Nilable<chem::NVector_sp> 	force,
                                                        bool		calcDiagonalHessian,
                                                        bool		calcOffDiagonalHessian,
                                                        gc::Nilable<chem::AbstractLargeSquareMatrix_sp>	hessian,
                                                        gc::Nilable<chem::NVector_sp>	hdvec,
                                                        gc::Nilable<chem::NVector_sp> dvec,
                                                        core::T_sp activeAtomMask,
                                                        core::T_sp debugInteractions )
{
  MAYBE_SETUP_ACTIVE_ATOM_MASK();
  MAYBE_SETUP_DEBUG_INTERACTIONS(debugInteractions.notnilp());
  this->_Evaluations++;
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY_CLEAR();
    LOG_ENERGY(("%s {\n") , this->className());
  }
  num_real totalEnergy = 0.0;
  ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool	hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());


//
// Copy from implementAmberFunction::evaluateAll
//
// -----------------------


#define CHIRAL_RESTRAINT_CALC_FORCE
#define CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) totalEnergy += (e);
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE 		ForceAcc
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE 	DiagHessAcc
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE OffDiagHessAcc


  {
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
    int	I1, I2, I3, I4, i;
    gctools::Vec0<EnergyChiralRestraint>::iterator cri;
    for ( i=0,cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++,i++ ) {
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
      if ( this->_Debug_NumberOfChiralRestraintTermsToCalculate > 0 ) {
        if ( i>= this->_Debug_NumberOfChiralRestraintTermsToCalculate ) {
          break;
        }
      }
#endif

			/* Obtain all the parameters necessary to calculate */
			/* the amber and forces */
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
      if (Energy>1.0) {
	auto as1 = atomLabel(cri->_Atom1);
        auto as2 = atomLabel(cri->_Atom2);
	auto as3 = atomLabel(cri->_Atom3);
	auto as4 = atomLabel(cri->_Atom4);
        core::clasp_write_string(fmt::format("TERM 7CHIRAL {} {} {:<9} {:<9} {:<9} {:<9} {:8.2f} {:8.2f} ; I1={} I2={} I3={} I4={}\n" , Energy, _rep_(cri->_Atom3), as1 , as2 , as3 , as4 , cri->term.K , cri->term.CO , cri->term.I1 , cri->term.I2 , cri->term.I3 , cri->term.I4 ));
      }
#if TURN_ENERGY_FUNCTION_DEBUG_ON //[
      cri->_calcForce = calcForce;
      cri->_calcDiagonalHessian = calcDiagonalHessian;
      cri ->_calcOffDiagonalHessian = calcOffDiagonalHessian;
		    // Now write all of the calculated values into the eval structure
#undef EVAL_SET
#define	EVAL_SET(var,val)	{ cri->eval.var=val;};
#include	<cando/chem/energy_functions/_ChiralRestraint_debugEvalSet.cc>
#endif //]

      if ( this->_DebugEnergy ) {
        LOG_ENERGY(( "MEISTER chiralRestraint %d args cando\n") , (i+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d K %lf\n") , (i+1) , K );
        LOG_ENERGY(( "MEISTER chiralRestraint %d x1 %5.3lf %d\n") , (i+1) , x1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d y1 %5.3lf %d\n") , (i+1) , y1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d z1 %5.3lf %d\n") , (i+1) , z1 , (I1/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d x2 %5.3lf %d\n") , (i+1) , x2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d y2 %5.3lf %d\n") , (i+1) , y2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d z2 %5.3lf %d\n") , (i+1) , z2 , (I2/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d x3 %5.3lf %d\n") , (i+1) , x3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d y3 %5.3lf %d\n") , (i+1) , y3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d z3 %5.3lf %d\n") , (i+1) , z3 , (I3/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d x4 %5.3lf %d\n") , (i+1) , x4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d y4 %5.3lf %d\n") , (i+1) , y4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d z4 %5.3lf %d\n") , (i+1) , z4 , (I4/3+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d results\n") , (i+1) );
        LOG_ENERGY(( "MEISTER chiralRestraint %d Energy %lf\n") , (i+1) , Energy);
        if ( calcForce ) {
//			LOG_ENERGY(( "MEISTER chiralRestraint %d DePhi %lf\n") , (i+1) , DePhi);
          LOG_ENERGY(( "MEISTER chiralRestraint %d fx1 %8.5lf %d\n") , (i+1) , fx1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fy1 %8.5lf %d\n") , (i+1) , fy1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fz1 %8.5lf %d\n") , (i+1) , fz1 , (I1/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fx2 %8.5lf %d\n") , (i+1) , fx2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fy2 %8.5lf %d\n") , (i+1) , fy2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fz2 %8.5lf %d\n") , (i+1) , fz2 , (I2/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fx3 %8.5lf %d\n") , (i+1) , fx3 , (I3/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fy3 %8.5lf %d\n") , (i+1) , fy3 , (I3/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fz3 %8.5lf %d\n") , (i+1) , fz3 , (I3/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fx4 %8.5lf %d\n") , (i+1) , fx4 , (I4/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fy4 %8.5lf %d\n") , (i+1) , fy4 , (I4/3+1) );
          LOG_ENERGY(( "MEISTER chiralRestraint %d fz4 %8.5lf %d\n") , (i+1) , fz4 , (I4/3+1) );
        }
        LOG_ENERGY(( "MEISTER chiralRestraint %d stop\n") , (i+1) );
      }
			/* Add the forces */

      if ( calcForce ) {
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz1>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz2>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz3>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fx4>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fy4>10000.0);
//		_lisp->profiler().eventCounter(core::forcesGreaterThan10000).recordCallAndProblem(fz4>10000.0);
      }
    }
  }
  if ( this->_DebugEnergy ) 
  {
    LOG_ENERGY(("%s }\n") , this->className());
  }
  maybeSetEnergy( componentEnergy, EnergyChiralRestraint_O::static_classSymbol(), totalEnergy );
  return totalEnergy;
}







void	EnergyChiralRestraint_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
		chem::NVector_sp 	pos)
{
  IMPLEMENT_ME();
#if 0
  int	fails = 0;
  bool	calcForce = true;
  bool	calcDiagonalHessian = true;
  bool	calcOffDiagonalHessian = true;


//
// copy from implementAmberFunction::compareAnalyticalAndNumericalForceAndHessianTermByTerm(
//
//------------------

#define CHIRAL_RESTRAINT_CALC_FORCE
#define CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#define CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


  {
		
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
    num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
    int	I1, I2, I3, I4, i;
    gctools::Vec0<EnergyChiralRestraint>::iterator cri;
    for ( i=0,cri=this->_Terms.begin();
          cri!=this->_Terms.end(); cri++,i++ ) {
			  /* Obtain all the parameters necessary to calculate */
			  /* the amber and forces */
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
      LOG("fx1 = {}" , fx1 );
      LOG("fy1 = {}" , fy1 );
      LOG("fz1 = {}" , fz1 );
      LOG("fx2 = {}" , fx2 );
      LOG("fy2 = {}" , fy2 );
      LOG("fz2 = {}" , fz2 );
      LOG("fx3 = {}" , fx3 );
      LOG("fy3 = {}" , fy3 );
      LOG("fz3 = {}" , fz3 );
      LOG("fx4 = {}" , fx4 );
      LOG("fy4 = {}" , fy4 );
      LOG("fz4 = {}" , fz4 );
      int index = i;
#include <cando/chem/energy_functions/_ChiralRestraint_debugFiniteDifference.cc>

    }
  }
#endif
}



#if 0
int	EnergyChiralRestraint_O::checkForBeyondThresholdInteractions(
			stringstream& info, chem::NVector_sp pos )
{
int	fails = 0;

    this->_BeyondThresholdTerms.clear();

//
// Copy from implementAmberFunction::checkForBeyondThresholdInteractions
//
//------------------


#undef CHIRAL_RESTRAINT_CALC_FORCE
#undef CHIRAL_RESTRAINT_CALC_DIAGONAL_HESSIAN
#undef CHIRAL_RESTRAINT_CALC_OFF_DIAGONAL_HESSIAN
#undef	CHIRAL_RESTRAINT_SET_PARAMETER
#define	CHIRAL_RESTRAINT_SET_PARAMETER(x)	{x = cri->term.x;}
#undef	CHIRAL_RESTRAINT_SET_POSITION
#define	CHIRAL_RESTRAINT_SET_POSITION(x,ii,of)	{x = pos->element(ii+of);}
#undef	CHIRAL_RESTRAINT_PHI_SET
#define	CHIRAL_RESTRAINT_PHI_SET(x) {}
#undef	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE
#define	CHIRAL_RESTRAINT_ENERGY_ACCUMULATE(e) {}
#undef	CHIRAL_RESTRAINT_FORCE_ACCUMULATE
#define	CHIRAL_RESTRAINT_FORCE_ACCUMULATE(i,o,v) {}
#undef	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}
#undef	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE
#define	CHIRAL_RESTRAINT_OFF_DIAGONAL_HESSIAN_ACCUMULATE(i1,o1,i2,o2,v) {}


    {
	    
#pragma clang diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#include <cando/chem/energy_functions/_ChiralRestraint_termDeclares.cc>
#pragma clang diagnostic pop
	num_real x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,K, CO;
	int	I1, I2, I3, I4, i;
        gctools::Vec0<EnergyChiralRestraint>::iterator cri;
	LOG("Entering checking loop, there are {} terms" , this->_Terms.end()-this->_Terms.begin() );
	for ( i=0,cri=this->_Terms.begin();
		    cri!=this->_Terms.end(); cri++,i++ ) {
		      /* Obtain all the parameters necessary to calculate */
		      /* the amber and forces */
	    LOG("Checking term# {}" , i  );
#include <cando/chem/energy_functions/_ChiralRestraint_termCode.cc>
	    LOG("Status" );
	    if ( ChiralTest>0.0 ) {
		chem::Atom_sp a1, a2, a3, a4;
		a1 = (*cri)._Atom1;
		a2 = (*cri)._Atom2;
		a3 = (*cri)._Atom3;
		a4 = (*cri)._Atom4;
	        LOG("Status" );
		info<< "ChiralRestraintDeviation ";
		info << "value " << ChiralTest << " Atoms(";
		info << a1->getName() << " ";
		info << a2->getName() << " ";
		info << a3->getName() << " ";
		info << a4->getName() << ")";
		info << std::endl;
	        LOG("Info: {}" , info.str().c_str()  );
		this->_BeyondThresholdTerms.push_back(*cri);
	        LOG("Status" );
		fails++;
	    }
	    LOG("Status" );
	}
	LOG("Status" );
    }

    return fails;
}

#endif









void EnergyChiralRestraint_O::initialize()
{
    this->Base::initialize();
    this->setErrorThreshold(0.2);
}


#ifdef XML_ARCHIVE
void EnergyChiralRestraint_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    IMPLEMENT_ME();
}
#endif


};
