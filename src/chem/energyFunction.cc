/*
    File: energyFunction.cc
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
__BEGIN_DOC(classes.EnergyFunction,section,EnergyFunction)
EnergyFunction objects are used to calculate Amber or Gaff molecular 
mechanics energies and first and second derivatives.
__END_DOC
*/

/*
 *	energyFunction.cc
 *
 */

#include <clasp/core/common.h>
#include <clasp/core/bformat.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/loop.h>
#include <cando/adapt/indexedObjectBag.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/environment.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/restraint.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/adapt/iterateCons.h>
#include <cando/chem/ringFinder.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyStretch.h>
#if USE_ALL_ENERGY_COMPONENTS
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyImproperRestraint.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyFixedNonbond.h>
#endif
#include <clasp/core/symbolTable.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffVdwDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/chemInfo.h>
#include <clasp/core/wrappers.h>


#if 0
#define	TURN_DEBUG_OFF
#define	MATHEMATICA_DEBUG_ON	1
#define	PARMDEBUG	1
#endif


#ifdef	darwin
#include <xmmintrin.h>
#endif



namespace chem
{

SYMBOL_EXPORT_SC_(ChemPkg,nonbondTerm);
SYMBOL_EXPORT_SC_(ChemPkg,restraintAnchor);

//Maybe use this in place of InteractionError?
SYMBOL_EXPORT_SC_(ChemPkg,InteractionError);
/*! LinearAngleError has a argument :atoms (atom1 atom2 atom3) */
SYMBOL_EXPORT_SC_(ChemPkg,LinearAngleError);

SYMBOL_EXPORT_SC_(ChemPkg,LinearDihedralError);
SYMBOL_EXPORT_SC_(ChemPkg,LinearImproperRestraintError);
SYMBOL_EXPORT_SC_(ChemPkg,OverlappingNonbondError);
SYMBOL_EXPORT_SC_(KeywordPkg,atoms);

bool energyFunctionInitialized = false;

int	_areValuesClose( double numVal, double analVal, const char* funcName, const char* termName, int index )
{
  double	rel = 0.0;
  if ( fabs(numVal)<FINITE_DIFFERENCE_TOLERANCE && fabs(analVal)<FINITE_DIFFERENCE_TOLERANCE ) {
    return true; // both close to zero
  }
  if ( numVal*analVal > 0.0 ) {
    double avg = (fabs(numVal)+fabs(analVal));
    rel = (fabs(numVal-analVal)/avg);
    if ( rel<FINITE_DIFFERENCE_TOLERANCE ) {
      return true;	// they are within a few percent of each other
    }
  }
  return false;
}

SYMBOL_EXPORT_SC_(ChemPkg,STARsecondaryAmideSmartsSTAR);
void energyFunction_initializeSmarts()
{
  if ( energyFunctionInitialized ) return;
  _sym_STARparserNodeHolderSTAR->defparameter(adapt::IndexedObjectBag_O::create());
  energyFunctionInitialized = true;
  SmartsRoot_sp secondaryAmide = chem__compile_smarts("[$([N:1]([C:2])(~[#1:3])~[C:4](=[O:5])[C:6])]",_Nil<core::T_O>());
  _sym_STARsecondaryAmideSmartsSTAR->defparameter(secondaryAmide);
}

CL_LAMBDA(matter &key (use-excluded-atoms t) active-atoms (assign-types t));
CL_LISPIFY_NAME(make_energy_function);
CL_DEF_CLASS_METHOD EnergyFunction_sp EnergyFunction_O::make(Matter_sp matter, bool useExcludedAtoms, core::T_sp activeAtoms, bool assign_types)
{
  GC_ALLOCATE(EnergyFunction_O, me );
  if ( matter.notnilp() ) {
    me->defineForMatter(matter,useExcludedAtoms,activeAtoms,assign_types);
  }
  return me;
};

bool inAtomSet(core::T_sp activeSet, Atom_sp a)
{
  if ( activeSet.nilp() ) return true;
  if ( core::Cons_sp cset = activeSet.asOrNull<core::Cons_O>() ) {
    for ( auto it : core::List_sp(cset) ) {
      if ( oCar(it) == a ) return true;
    }
    return false;
  } else if ( core::HashTable_sp htset = activeSet.asOrNull<core::HashTable_O>() ) {
    return htset->gethash(a).notnilp();
  } else if ( core::Vector_sp vset = activeSet.asOrNull<core::Vector_O>() ) {
    for ( size_t i(0), iEnd(core::cl__length(vset)); i<iEnd; ++i ) {
      if ( vset->rowMajorAref(i) == a ) return true;
    }
    return false;
  }
  SIMPLE_ERROR(BF("Add support to search %s as a set\n") % _rep_(activeSet));
}



void	EnergyFunction_O::initialize()
{
  this->Base::initialize();
  this->_AtomTable = AtomTable_O::create();
  this->_Stretch = EnergyStretch_O::create();
#if USE_ALL_ENERGY_COMPONENTS
  this->_Angle = EnergyAngle_O::create();
  this->_Dihedral = EnergyDihedral_O::create();
  this->_Nonbond = EnergyNonbond_O::create();
  this->_ChiralRestraint = EnergyChiralRestraint_O::create();
  this->_AnchorRestraint = EnergyAnchorRestraint_O::create();
  this->_ImproperRestraint = EnergyImproperRestraint_O::create();
  this->_FixedNonbondRestraint = EnergyFixedNonbondRestraint_O::create();
#endif
  this->setScoringFunctionName("");
  this->_Message = "";
  this->useDefaultSettings();
}

void	EnergyFunction_O::useDefaultSettings()
{
  ASSERTNOTNULL(this->_Stretch);
  ASSERTNOTNULL(this->_Stretch);
  ASSERT(this->_Stretch.notnilp());
  this->_Stretch->initialize();
#if USE_ALL_ENERGY_COMPONENTS
  ASSERTNOTNULL(this->_Angle);
  ASSERT(this->_Angle.notnilp());
  ASSERTNOTNULL(this->_Dihedral);
  ASSERT(this->_Dihedral.notnilp());
  ASSERTNOTNULL(this->_Nonbond);
  ASSERT(this->_Nonbond.notnilp());
  ASSERTNOTNULL(this->_ChiralRestraint);
  ASSERT(this->_ChiralRestraint.notnilp());
  ASSERTNOTNULL(this->_AnchorRestraint);
  ASSERT(this->_AnchorRestraint.notnilp());
  ASSERTNOTNULL(this->_ImproperRestraint);
  ASSERT(this->_ImproperRestraint.notnilp());
  ASSERTNOTNULL(this->_FixedNonbondRestraint);
  ASSERT(this->_FixedNonbondRestraint.notnilp());
  this->_Angle->initialize();
  this->_Dihedral->initialize();
  this->_Nonbond->initialize();
  this->_ChiralRestraint->initialize();
  this->_AnchorRestraint->initialize();
  this->_ImproperRestraint->initialize();
  this->_FixedNonbondRestraint->initialize();
#endif
  this->_ChiralRestraintWeight = DefaultChiralRestraintWeight;
  this->_ChiralRestraintOffset = DefaultChiralRestraintOffset;
  this->_AnchorRestraintWeight = DefaultAnchorRestraintWeight;
  this->_RestrainSecondaryAmides = true;
}



void EnergyFunction_O::fields(core::Record_sp node)
{
  node->field_if_not_unbound(INTERN_(kw,AtomTable),this->_AtomTable);
  node->field_if_not_unbound(INTERN_(kw,Stretch),this->_Stretch);
#if USE_ALL_ENERGY_COMPONENTS
  node->field_if_not_unbound(INTERN_(kw,Angle),this->_Angle);
  node->field_if_not_unbound(INTERN_(kw,Dihedral),this->_Dihedral);
  node->field_if_not_unbound(INTERN_(kw,Nonbond),this->_Nonbond);
  node->field_if_not_unbound(INTERN_(kw,ImproperRestraint),this->_ImproperRestraint);
  node->field_if_not_unbound(INTERN_(kw,ChiralRestraint),this->_ChiralRestraint);
  node->field_if_not_unbound(INTERN_(kw,AnchorRestraint),this->_AnchorRestraint);
  node->field_if_not_unbound(INTERN_(kw,FixedNonbondRestraint),this->_FixedNonbondRestraint);
#endif
}



#ifdef XML_ARCHIVE
void	EnergyFunction_O::archive(core::ArchiveP node)
{
  node->attribute("_Name",this->_Name);
  node->attribute("_Matter",this->_Matter);
  node->attribute("_AtomTable",this->_AtomTable );
//    node->archiveSharedObject<Dumb_StretchComponent>("_Stretch","StretchTerms",this->_Stretch);
//    node->archiveSharedObject<Dumb_AngleComponent>("_Angle","AngleTerms",this->_Angle);
//    node->archiveSharedObject<Dumb_DihedralComponent>("_Dihedral","DihedralTerms",this->_Dihedral);
//    node->archiveSharedObject<Dumb_ImproperRestraintComponent>("_ImproperRestraint","ImproperRestraintTerms",this->_ImproperRestraint);
//    node->archiveSharedObject<Dumb_ChiralRestraintComponent>("_ChiralRestraint","ChiralRestraintTerms",this->_ChiralRestraint);
//    node->archiveSharedObject<Dumb_AnchorRestraintComponent>("_AnchorRestraint","AnchorRestraintTerms",this->_AnchorRestraint);
//    node->archiveSharedObject<Dumb_NonbondComponent>("_Nonbond","NonbondTerms",this->_Nonbond);
  node->attribute("_Stretch",this->_Stretch);
#if USE_ALL_ENERGY_COMPONENTS
  node->attribute("_Angle",this->_Angle);
  node->attribute("_Dihedral",this->_Dihedral);
  node->attribute("_Nonbond",this->_Nonbond);
  node->attribute("_ImproperRestraint",this->_ImproperRestraint);
  node->attribute("_ChiralRestraint",this->_ChiralRestraint);
  node->attribute("_AnchorRestraint",this->_AnchorRestraint);
  node->attribute("_FixedNonbondRestraint",this->_FixedNonbondRestraint);
#endif
  node->attribute("_ChiralRestraintWeight",this->_ChiralRestraintWeight);
  node->attribute("_ChiralRestraintOffset",this->_ChiralRestraintOffset);
  node->attribute("_AnchorRestraintWeight",this->_AnchorRestraintWeight);
  node->attribute("_TotalEnergy",this->_TotalEnergy);
  node->archiveString("_Message",this->_Message);
  node->attribute("_DielectricConstant",this->_DielectricConstant);

}
#endif



size_t EnergyFunction_O::getNVectorSize() 
{ 
  return this->_AtomTable->getNVectorSize();
};

#if 0
double EnergyFunction_O::getDihedralComponentEnergy()
{
  return this->_Dihedral->getEnergy();
}

double EnergyFunction_O::getImproperRestraintComponentEnergy()
{
  return this->_ImproperRestraint->getEnergy();
}
#endif








void	EnergyFunction_O::setOption( core::Symbol_sp option, core::T_sp val)
{_OF();
  if ( option == _sym_nonbondTerm)
  {
    bool onval = val.isTrue();
#if USE_ALL_ENERGY_COMPONENTS
    if ( onval ) this->getNonbondComponent()->enable();
    else	this->getNonbondComponent()->disable();
#endif
  }
  else if ( option == _sym_restraintAnchor )
  {
    bool onval = val.isTrue();
#if USE_ALL_ENERGY_COMPONENTS
    if ( onval ) this->getAnchorRestraintComponent()->enable();
    else	this->getAnchorRestraintComponent()->disable();
#endif
  } else
  {
    SIMPLE_ERROR(BF("Unknown EnergyFunction setOption keyword[%s]") % _rep_(option) );
  }
}



/*
__BEGIN_DOC(classes.EnergyFunction.setOptions,subsection,setOptions)
\scriptCmd{setOptions}{energyFunction options}

Set the \scriptArg{options} for the energy function. The options are represented by a list of list pairs
of the form [: 'symbol value ]
  They include...\par
[ 'nonbond true/false ]\par
__END_DOC
*/
void	EnergyFunction_O::setOptions( core::List_sp options )
{
  while ( options.notnilp() )
  {
    core::Symbol_sp option = core::oCar(options).as<core::Symbol_O>();
    options = options.asCons()->cdr();
    core::T_sp val = core::oCar(options);
    options = options.asCons()->cdr();
    this->setOption(option,val);
  }
}




// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
//
//	This is where we include the code that interfaces to Mathematica
//	generated code.
//
//

class	FiniteDifferenceMismatch {
public:
  string	functionName;
  string	termName;
  int	index;
  string	getError() { return "FiniteDifferenceMismatch"; } //  @"+this->functionName()+" term:("+this->termName+")";
};



#define	delta2	0.00001



#undef	TEST_FORCE
#define	TEST_FORCE( func, delta, argLow, argHigh, term, idx ) {\
    double eLow = this->func argLow;\
    double eHigh = this->func argHigh;\
    double numForce = -(eHigh-eLow)/(delta);\
    if ( !_areValuesClose( numForce, term, #func, #term, idx ) ) { fails++; };\
  }

#undef	TEST_DIAGONAL_HESSIAN
#define	TEST_DIAGONAL_HESSIAN( func, delta, argLow, argMiddle, argHigh, term, idx ) {\
    double eLow = this->func argLow;\
    double eMiddle = this->func argMiddle;\
    double eHigh = this->func argHigh;\
    double numHessian = (eHigh+eLow-2.0*eMiddle)/(((delta)/2.0)*((delta)/2.0));\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }

#undef	TEST_OFF_DIAGONAL_HESSIAN
#define	TEST_OFF_DIAGONAL_HESSIAN( func, delta, argMM, argPM, argMP, argPP, term, idx ) {\
    double eMM = this->func argMM;\
    double eMP = this->func argMP;\
    double ePM = this->func argPM;\
    double ePP = this->func argPP;\
    double numHessian = ((ePP-ePM)-(eMP-eMM))/(delta*delta);\
    if ( !_areValuesClose( numHessian, term, #func, #term, idx ) ) { fails++; };\
  }




//
// Interface to amberFunction.cc
//

void EnergyFunction_O::setupHessianPreconditioner(NVector_sp nvPosition,
                                                               AbstractLargeSquareMatrix_sp m )
{
#ifdef ENERGY_FUNCTION_TIMER
//    _lisp->profiler().timer(core::timerPreconditioner).start();
//    _lisp->profiler().timer(core::timerPreconditionerSetup).start();
#endif
  m->fill(0.0);

  this->_Stretch->setupHessianPreconditioner(nvPosition, m );
#if USE_ALL_ENERGY_COMPONENTS
  this->_Angle->setupHessianPreconditioner(nvPosition, m );
  this->_Dihedral->setupHessianPreconditioner(nvPosition, m );
		// Nonbond doesn't contribute to hessian preconditioner
//    this->_Nonbond->setupHessianPreconditioner(nvPosition, m );
  this->_ChiralRestraint->setupHessianPreconditioner(nvPosition, m );
  this->_AnchorRestraint->setupHessianPreconditioner(nvPosition, m );
  this->_ImproperRestraint->setupHessianPreconditioner(nvPosition, m );
//    this->_FixedNonbondRestraint->setupHessianPreconditioner(nvPosition, m );
#endif

#ifdef ENERGY_FUNCTION_TIMER
//    _lisp->profiler().timer(core::timerPreconditionerSetup).stop();
//    _lisp->profiler().timer(core::timerPreconditioner).stop();
#endif
}







uint	EnergyFunction_O::countTermsBeyondThreshold()
{_OF();
  int		terms;
  terms = 0;
  SIMPLE_ERROR(BF("Should there be something here?"));
//    node = rawAccumulateTermsBeyondThresholdAsXml(terms);
  return terms;
}






//

double	EnergyFunction_O::evaluateAll(
                                      NVector_sp 	pos,
                                      bool 		calcForce,
                                      gc::Nilable<NVector_sp> 	force,
                                      bool		calcDiagonalHessian,
                                      bool		calcOffDiagonalHessian,
                                      gc::Nilable<AbstractLargeSquareMatrix_sp>	hessian,
                                      gc::Nilable<NVector_sp>	hdvec,
                                      gc::Nilable<NVector_sp> dvec)
{_G()
#ifdef DEBUG_ENERGY_FUNCTION
    printf("%s:%d:%s Entered\n", __FILE__, __LINE__, __FUNCTION__ );
#endif
#if 0 //[
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
      hdvec->addToElement((i1)+(o1),v*dvec->element((i2)+(o2)));\
      hdvec->addToElement((i2)+(o2),v*dvec->element((i1)+(o1)));\
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
      hdvec->addToElement((i1)+(o1),v*dvec->element((i1)+(o1)));\
    }\
  }
#endif //]

    ANN(force);
  ANN(hessian);
  ANN(hdvec);
  ANN(dvec);
  bool	hasForce = force.notnilp();
  bool   hasHessian = hessian.notnilp();
  bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

#ifdef	DEBUG_ON //[
	// Summarize entry state for debugging
  LOG(BF("calcForce = %d") % calcForce  );
  LOG(BF("calcDiagonalHessian = %d") % calcDiagonalHessian  );
  LOG(BF("calcOffDiagonalHessian = %d") % calcOffDiagonalHessian  );
  LOG(BF("hasForce = %d") % hasForce  );
  LOG(BF("hasHdAndD = %d") % hasHdAndD  );
  if ( hasForce && force->size() < pos->size() ) {
    SIMPLE_ERROR(BF("Force does not have the necessary dimensions"));
  }
#endif //]

  if ( !calcForce && ( calcDiagonalHessian || calcOffDiagonalHessian ) ) {
    SIMPLE_ERROR(BF("Inconsistant arguments: if you want to calcDiagonalHessian or calcOffDiagonalHessian you must calcForce"));
  }
  if ( !calcDiagonalHessian & calcOffDiagonalHessian ) {
    SIMPLE_ERROR(BF("Inconsistant arguments: if you want to calcOffDiagonalHessian you must calcDiagonalHessian"));
  }

////    _lisp->profiler().pushTimerStates();
  ALL_ENERGY_COMPONENTS(zeroEnergy());

    if ( hasForce ) force->zero();
    if ( hasHessian ) hessian->zero();
    if ( hasHdAndD ) {
      LOG(BF("Zeroing hdvec") );
      hdvec->zero();	// Zero the result
    }

    LOG(BF("Starting evaluation of energy") );
////	_lisp->profiler().timer(core::timerEnergy).start();

////	_lisp->profiler().timer(core::timerBondAngleDihedral).start();

	//
	// Evaluate the stretch term
	//
////	_lisp->profiler().timer(core::timerBond).start();
    this->_Stretch->evaluateAll( pos, calcForce, force,
                                 calcDiagonalHessian,
                                 calcOffDiagonalHessian,
                                 hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerBond).stop();
#if USE_ALL_ENERGY_COMPONENTS
////	_lisp->profiler().timer(core::timerAngle).start();
    this->_Angle->evaluateAll( pos, calcForce, force,
                               calcDiagonalHessian,
                               calcOffDiagonalHessian,
                               hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerAngle).stop();

////	_lisp->profiler().timer(core::timerDihedral).start();
    this->_Dihedral->evaluateAll( pos, calcForce, force,
                                  calcDiagonalHessian,
                                  calcOffDiagonalHessian,
                                  hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerDihedral).stop();
////	_lisp->profiler().timer(core::timerBondAngleDihedral).stop();

////	_lisp->profiler().timer(core::timerNonbond).start();
    this->_Nonbond->evaluateAll( pos, calcForce, force,
                                 calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerNonbond).stop();

////	_lisp->profiler().timer(core::timerImproperRestraint).start();
    this->_ImproperRestraint->evaluateAll( pos, calcForce, force,
                                           calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerImproperRestraint).stop();

////	_lisp->profiler().timer(core::timerChiralRestraint).start();
    this->_ChiralRestraint->evaluateAll( pos, calcForce, force,
                                         calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerChiralRestraint).stop();

//	_lisp->profiler().timer(core::timerAnchorRestraint).start();
    this->_AnchorRestraint->evaluateAll( pos, calcForce, force,
                                         calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerAnchorRestraint).stop();

////	_lisp->profiler().timer(core::timerFixedNonbondRestraint).start();
    this->_FixedNonbondRestraint->evaluateAll( pos, calcForce, force,
                                               calcDiagonalHessian, calcOffDiagonalHessian, hessian, hdvec, dvec );
////	_lisp->profiler().timer(core::timerFixedNonbondRestraint).stop();
#endif
    this->_TotalEnergy = this->_Stretch->getEnergy();
#if USE_ALL_ENERGY_COMPONENTS
    this->_TotalEnergy += this->_Angle->getEnergy();
    this->_TotalEnergy += this->_Dihedral->getEnergy();
    this->_TotalEnergy += this->_Nonbond->getEnergy();
    this->_TotalEnergy += this->_ImproperRestraint->getEnergy();
    this->_TotalEnergy += this->_ChiralRestraint->getEnergy();
    this->_TotalEnergy += this->_AnchorRestraint->getEnergy();
    this->_TotalEnergy += this->_FixedNonbondRestraint->getEnergy();
#endif

////	_lisp->profiler().timer(core::timerEnergy).stop();

    	// More energy terms
  return this->_TotalEnergy;
}


string EnergyFunction_O::energyComponentsAsString()
{
  stringstream ss;
  ss << boost::format("Stretch(%lf)") % this->_Stretch->getEnergy() << std::endl;
#if USE_ALL_ENERGY_COMPONENTS
  ss << boost::format("Angle(%lf)") % this->_Angle->getEnergy() << std::endl;
  ss << boost::format("Dihedral(%lf)") % this->_Dihedral->getEnergy() << std::endl;
  ss << boost::format("Nonbond(%lf)") % this->_Nonbond->getEnergy() << std::endl;
  ss << boost::format("ChiralRestraint(%lf)") % this->_ChiralRestraint->getEnergy() << std::endl;
  ss << boost::format("ImproperRestraint(%lf)") % this->_ImproperRestraint->getEnergy() << std::endl;
  ss << boost::format("AnchorRestraint(%lf)") % this->_AnchorRestraint->getEnergy() << std::endl;
  ss << boost::format("FixedNonbondRestraint(%lf)") % this->_FixedNonbondRestraint->getEnergy() << std::endl;
#endif
  return ss.str();
}


/*!
 * Compare the analytical force and hessian components term by term with
 * numerical ones.  Print a message for every mismatch
 * Return the number of mismatches
 */
int	EnergyFunction_O::compareAnalyticalAndNumericalForceAndHessianTermByTerm( NVector_sp 	pos)
{
  int	fails = 0;

  {
    this->_Stretch->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
#if USE_ALL_ENERGY_COMPONENTS
    this->_Angle->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_Dihedral->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_Nonbond->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_ImproperRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_ChiralRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_AnchorRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
    this->_FixedNonbondRestraint->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
#endif
  }
  return fails;
}


/*!
 * Compare the analytical force and hessian components term by term with
 * numerical ones at the current position.  Print a message for every mismatch
 * Return the number of mismatches
 */
int	EnergyFunction_O::compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition( )
{
  NVector_sp pos = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  return this->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
}




SYMBOL_EXPORT_SC_(KeywordPkg,stretch_deviations);
SYMBOL_EXPORT_SC_(KeywordPkg,angle_deviations);

/*!
 * Look for bad geometry
 *
 * Extract the coordinates of the atoms
 */
CL_DEFMETHOD
core::List_sp EnergyFunction_O::checkForBeyondThresholdInteractions(double threshold)
{
  NVector_sp	pos;
  stringstream	info;
  int	fails = 0;

  info.str("");
  pos = NVector_O::create(this->getNVectorSize());
  this->loadCoordinatesIntoVector(pos);
  fails = 0;
  ql::list result;
  result & this->_Stretch->checkForBeyondThresholdInteractionsWithPosition(pos,threshold);
  result & this->_Angle->checkForBeyondThresholdInteractionsWithPosition(pos,threshold);
  result & this->_Nonbond->checkForBeyondThresholdInteractionsWithPosition(pos,threshold);
  return result.cons();
}



//
//
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------
// ------------------------------------------------------------------------





string	EnergyFunction_O::energyTermsEnabled()
{
  stringstream ss;
  ss.str("");
  ss << this->_Stretch->enabledAsString();
#if USE_ALL_ENERGY_COMPONENTS
  ss << this->_Angle->enabledAsString();
  ss << this->_Dihedral->enabledAsString();
  ss << this->_Nonbond->enabledAsString();
  ss << this->_ImproperRestraint->enabledAsString();
  ss << this->_ChiralRestraint->enabledAsString();
  ss << this->_AnchorRestraint->enabledAsString();
  ss << this->_FixedNonbondRestraint->enabledAsString();
#endif
  return ss.str();
}




#define	mathematicaSinNPhi(n,s,c) 0.0
#define	mathematicaCosNPhi(n,s,c) 0.0













#define	DELTA	0.00000001

double	EnergyFunction_O::calculateNumericalDerivative(NVector_sp pos, double delta, uint i )
{
  double x, ylow, yhigh, fval;
  double	deltaDiv2 = delta/2.0;
  x = pos->element(i);
  pos->setElement(i,x-deltaDiv2);
  ylow = this->evaluateEnergy(pos);
  pos->setElement(i,x+deltaDiv2);
  yhigh = this->evaluateEnergy(pos);
  pos->setElement(i,x);
  fval = (yhigh-ylow)/delta;
  return fval;
}


double	EnergyFunction_O::calculateNumericalSecondDerivative(NVector_sp pos, double delta, uint i, uint j )
{
  double	x, fxmh, fx, fxph, f2;
  double	y, fpipj, fpimj, fmipj, fmimj, fp, fm;
  if ( i==j ) {
    x = pos->element(i);
    pos->setElement(i,x-delta);
    fxmh = this->evaluateEnergy(pos);
    pos->setElement(i,x+delta);
    fxph = this->evaluateEnergy(pos);
    pos->setElement(i,x);
    fx = this->evaluateEnergy(pos);
    f2 = (fxph+fxmh-2.0*(fx))/(delta*delta);
  } else {
    double	deltaDiv2 = delta/2.0;
    x = pos->element(i);
    y = pos->element(j);
    pos->setElement(i,x+deltaDiv2);
    pos->setElement(j,y+deltaDiv2);
    fpipj = this->evaluateEnergy(pos);
    pos->setElement(i,x+deltaDiv2);
    pos->setElement(j,y-deltaDiv2);
    fpimj = this->evaluateEnergy(pos);
    pos->setElement(i,x-deltaDiv2);
    pos->setElement(j,y+deltaDiv2);
    fmipj = this->evaluateEnergy(pos);
    pos->setElement(i,x-deltaDiv2);
    pos->setElement(j,y-deltaDiv2);
    fmimj = this->evaluateEnergy(pos);
    pos->setElement(i,x);
    pos->setElement(j,y);
    LOG(BF("fpipj = %le") % fpipj  );
    LOG(BF("fpimj = %le") % fpimj  );
    LOG(BF("fmipj = %le") % fmipj  );
    LOG(BF("fmimj = %le") % fmimj  );
    fp = (fpipj-fpimj)/delta;
    fm = (fmipj-fmimj)/delta;
    LOG(BF("fp = %le") % fp  );
    LOG(BF("fm = %le") % fm  );
    f2 = (fp-fm)/delta;
    LOG(BF("f2 = %le") % f2  );
  }
  return f2;
}




/*! Calculate the force numerically
 */
void	EnergyFunction_O::evaluateNumericalForce(NVector_sp pos, NVector_sp numForce, double delta)
{
  double		fval;
  uint		i;

  for (i=0; i<pos->size(); i++ ) {
    fval = -this->calculateNumericalDerivative(pos,delta,i);
    numForce->setElement(i,fval);
  }
}


/*! Calculate the hessian numerically
 */
void	EnergyFunction_O::evaluateNumericalHessian(NVector_sp pos, AbstractLargeSquareMatrix_sp hessian, bool calcOffDiagonal, double delta )
{_OF();
  double		fval;
  uint		c, r;

  if ( hessian->columns() != pos->size() || hessian->rows()!=pos->size() ) {
    SIMPLE_ERROR(BF("evaluateNumericalHessian must have the right size"));
  }
  hessian->zero();
  for ( c=0; c<pos->size(); c++ ) {
    fval = this->calculateNumericalSecondDerivative(pos,delta,c,c);
    hessian->setElement(c,c,fval);
  }
  if ( !calcOffDiagonal ) return;
  for ( c=0; c<pos->size(); c++ ) {
    for ( r=0; r<pos->size(); r++ ) {
      if ( c!=r) {
        fval = this->calculateNumericalSecondDerivative(pos,delta,c,r);
        hessian->setElement(c,r,fval);
      }
    }
  }
}


/*!
 * Check if the force matches the numerical force
 *
 * Return an XML node that indicates what is going on.
 * If there is a mis-match then dump the EnergyFunction into the result.
 *
 */
ForceMatchReport_sp EnergyFunction_O::checkIfAnalyticalForceMatchesNumericalForce(NVector_sp pos, NVector_sp analyticalForce )
{
  ForceMatchReport_sp report;
  NVector_sp	numForce, tempForce;
  stringstream	result;
  double	analyticalMag, numericalMag, dot;
  double	avg;

  report = ForceMatchReport_O::create();

  numForce = NVector_O::create(pos->size());
  this->evaluateNumericalForce(pos,numForce,DELTA);
  dot = dotProduct(numForce,analyticalForce);
  numericalMag = magnitude(numForce);
  analyticalMag = magnitude(analyticalForce);
  tempForce = NVector_O::create(pos->size());
    	// Evaluate the force at pos again
  this->evaluateEnergyForce(pos,true,tempForce);
  avg = (analyticalMag+numericalMag)/2.0;
  if ( analyticalMag < VERYSMALL && numericalMag < VERYSMALL ) {
    result.str("");
    result << "analyticalForce and numericalForces are both < VERYSMALL"<<std::endl;
    report->_Message = result.str();
    goto DONE;
  }
  if ( analyticalMag < VERYSMALL ) {
    result.str("");
    result << "analyticalForce < VERYSMALL and numericalForces is not"<<std::endl;
    report->_Message = result.str();
    goto DONE;
  }
  if ( numericalMag < VERYSMALL ) {
    result.str("");
    result << "numericalForce < VERYSMALL and analyticalForce is not"<<std::endl;
    report->_Message = result.str();
    goto DONE;
  }
  if ( avg < VERYSMALL ) {
    report->_Message = "average of Analytical & Numerical Forces is VERY small";
    goto DONE;
  }
  dot /= ( analyticalMag*numericalMag );

    	// Dump the analytical and the numerical force vectors to the log
	//
  report->_AnalyticalForce = analyticalForce;
  report->_NumericalForce = numForce;
  if ( fabs(analyticalMag-numericalMag)/avg >0.1 ) {
    result.str("");
    result << "Lengths of Analytical and Numerical forces differ by more than 10%" << std::endl;
    result << "|analyticalForce| == " << analyticalMag << "  |numericalForce| == "<< numericalMag << std::endl;
    result << "(analyticalForce/|analyticalForce|).(numericalForce/|numericalForce|) = "<< dot << std::endl;
    report->_Message = result.str();
    this->saveCoordinatesAndForcesFromVectors(pos,analyticalForce);
    goto DONE;
  }
  if ( dot < 0.98 ) {
    result.str("");
    result << "The Lengths of Analytical and Numerical forces are very similar but they are not parallel%" << std::endl;
    result << "|analyticalForce| == " << analyticalMag << "  |numericalForce| == "<< numericalMag << std::endl;
    result << "(analyticalForce/|analyticalForce|).(numericalForce/|numericalForce|) = "<< dot << std::endl;
    report->_Message = result.str();
    this->saveCoordinatesAndForcesFromVectors(pos,analyticalForce);
    goto DONE;
  }
  report->_Message = "Analytical and Numerical forces are virtually identical";
 DONE:
  return report;
}



void	EnergyFunction_O::summarizeTerms()
{
  core::write_bf_stream(BF("Number of atom terms: %d\n") % this->_AtomTable->getNumberOfAtoms() );
  core::write_bf_stream(BF("Number of Stretch terms: %d\n") % this->_Stretch->numberOfTerms() );
  core::write_bf_stream(BF("Number of Angle terms: %d\n") % this->_Angle->numberOfTerms() );
  core::write_bf_stream(BF("Number of Dihedral terms: %d\n") % this->_Dihedral->numberOfTerms() );
};



void	EnergyFunction_O::dumpTerms()
{
  this->_AtomTable->dumpTerms();
  this->_Stretch->dumpTerms();
#if USE_ALL_ENERGY_COMPONENTS
  this->_Angle->dumpTerms();
  this->_Dihedral->dumpTerms();
  this->_Nonbond->dumpTerms();
  this->_ImproperRestraint->dumpTerms();
  this->_ChiralRestraint->dumpTerms();
  this->_AnchorRestraint->dumpTerms();
  this->_FixedNonbondRestraint->dumpTerms();
#endif
}





int EnergyFunction_O::_applyRestraints(core::T_sp nonbondDb, core::Iterator_sp restraintIterator, core::T_sp activeAtoms )
{
  int terms = 0;
#if USE_ALL_ENERGY_COMPONENTS
  restraintIterator->first();
  while ( restraintIterator->notDone() ) {
    Restraint_sp restraint = restraintIterator->current<Restraint_O>();
    if ( !restraint->isActive() ) goto CONT;
    if ( restraint.isA<RestraintDihedral_O>() )
    {
      RestraintDihedral_sp dih = (restraint).as<RestraintDihedral_O>();
      EnergyImproperRestraint energyTerm;
      if ( activeAtoms.notnilp() &&
           (!inAtomSet(activeAtoms,dih->getAtomA())
            || !inAtomSet(activeAtoms,dih->getAtomB())
            || !inAtomSet(activeAtoms,dih->getAtomC())
            || !inAtomSet(activeAtoms,dih->getAtomD()) ) ) goto CONT;
      energyTerm._Atom1 = dih->getAtomA();
      energyTerm._Atom2 = dih->getAtomB();
      energyTerm._Atom3 = dih->getAtomC();
      energyTerm._Atom4 = dih->getAtomD();
      EnergyAtom* ea1 = this->getEnergyAtomPointer(dih->getAtomA());
      EnergyAtom* ea2 = this->getEnergyAtomPointer(dih->getAtomB());
      EnergyAtom* ea3 = this->getEnergyAtomPointer(dih->getAtomC());
      EnergyAtom* ea4 = this->getEnergyAtomPointer(dih->getAtomD());
      energyTerm.term.I1 = ea1->coordinateIndexTimes3();
      energyTerm.term.I2 = ea2->coordinateIndexTimes3();
      energyTerm.term.I3 = ea3->coordinateIndexTimes3();
      energyTerm.term.I4 = ea4->coordinateIndexTimes3();
      energyTerm.term.U = dih->getMaxDegrees()*0.0174533;
      energyTerm.term.L = dih->getMinDegrees()*0.0174533;
      energyTerm.term.K = dih->getWeight();
      this->_ImproperRestraint->addTerm(energyTerm);
      ++terms;
    } else if ( restraint.isA<RestraintAnchor_O>() )
    {
      RestraintAnchor_sp anchor = (restraint).as<RestraintAnchor_O>();
      EnergyAnchorRestraint	iterm;
      Vector3		anchorPos;
      Atom_sp a1 = anchor->getAtom();
      if ( activeAtoms.notnilp() && !inAtomSet(activeAtoms,a1) ) goto CONT;
      EnergyAtom* ea1 = this->getEnergyAtomPointer(a1);
      LOG(BF("Create an anchor restraint for %s") % a1->description()  );
      anchorPos = anchor->getAnchorPos();
      iterm.term.xa = anchorPos.getX();
      iterm.term.ya = anchorPos.getY();
      iterm.term.za = anchorPos.getZ();
      iterm.term.ka = this->_AnchorRestraintWeight;
      iterm.term.I1 = ea1->coordinateIndexTimes3();
      this->_AnchorRestraint->addTerm(iterm);
      ++terms;
    } else if ( restraint.isA<RestraintFixedNonbond_O>() )
    {
      this->_FixedNonbondRestraint->setupForEvaluation(this->_AtomTable,this->_NonbondCrossTermTable);
      RestraintFixedNonbond_sp fixedNonbond = restraint.as<RestraintFixedNonbond_O>();
      Matter_sp matter = fixedNonbond->getMatter();
//	    EnergyAtom	energyAtom(_lisp);
      Loop loop;
      {_BLOCK_TRACE("Defining FixedNonbond ATOMS");
        loop.loopTopGoal(matter,ATOMS);
        while ( loop.advanceLoopAndProcess() ) 
        {
          Atom_sp a1 = loop.getAtom();
          if ( activeAtoms.notnilp() && !inAtomSet(activeAtoms,a1))  goto CONT;
          if ( a1.isA<VirtualAtom_O>() ) continue; // skip virtuals
          this->_FixedNonbondRestraint->addFixedAtom(nonbondDb,a1);
          ++terms;
        }
      }
    } else if ( restraint.isA<RestraintDistance_O>() ) {
      RestraintDistance_sp rd = gc::As_unsafe<RestraintDistance_sp>(restraint);
      EnergyAtom* ea1 = this->getEnergyAtomPointer(rd->_A);
      EnergyAtom* ea2 = this->getEnergyAtomPointer(rd->_B);
      EnergyStretch   energyStretch(rd->_A,rd->_B,ea1->coordinateIndexTimes3(),ea2->coordinateIndexTimes3(),rd->_K,rd->_R0);
      this->_Stretch->addTerm(energyStretch);
    } else {
      SIMPLE_ERROR(BF("Handle restraint: %s") % _rep_(restraint));
    }
  CONT:
    restraintIterator->next();
  }
#endif
  return terms;
}

void EnergyFunction_O::_applyDihedralRestraint(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double minDegrees, double maxDegrees, double weight, core::T_sp activeAtoms)
{_OF();
  if ( activeAtoms.notnilp() &&
       (!inAtomSet(activeAtoms,a1) || !inAtomSet(activeAtoms,a2) || !inAtomSet(activeAtoms,a3) || !inAtomSet(activeAtoms,a4)) ) return;
#if USE_ALL_ENERGY_COMPONENTS
  EnergyImproperRestraint energyTerm;
  energyTerm._Atom1 = a1;
  energyTerm._Atom2 = a2;
  energyTerm._Atom3 = a3;
  energyTerm._Atom4 = a4;
  EnergyAtom* ea1 = this->getEnergyAtomPointer(energyTerm._Atom1);
  EnergyAtom* ea2 = this->getEnergyAtomPointer(energyTerm._Atom2);
  EnergyAtom* ea3 = this->getEnergyAtomPointer(energyTerm._Atom3);
  EnergyAtom* ea4 = this->getEnergyAtomPointer(energyTerm._Atom4);
  energyTerm.term.I1 = ea1->coordinateIndexTimes3();
  energyTerm.term.I2 = ea2->coordinateIndexTimes3();
  energyTerm.term.I3 = ea3->coordinateIndexTimes3();
  energyTerm.term.I4 = ea4->coordinateIndexTimes3();
  energyTerm.term.U = maxDegrees*0.0174533;
  energyTerm.term.L = minDegrees*0.0174533;
  energyTerm.term.K = weight;
  this->_ImproperRestraint->addTerm(energyTerm);
#endif
}

void EnergyFunction_O::__createSecondaryAmideRestraints(gctools::Vec0<Atom_sp>& nitrogens, core::T_sp activeAtoms )
{
  gctools::Vec0<Atom_sp>::iterator ni;
  double transMin = -160.0;
  double transMax = 160.0;
  double cisMin = 20.0;
  double cisMax = -20.0;
  double weight = 1.0;
  for ( ni=nitrogens.begin(); ni!=nitrogens.end(); ni++ )
  {
    SmartsRoot_sp secondaryAmide = gctools::As<SmartsRoot_sp>(chem::_sym_STARsecondaryAmideSmartsSTAR->symbolValue());
    core::T_mv match_mv = chem__chem_info_match(secondaryAmide,*ni);
    if ( match_mv.notnilp() ) {
      ChemInfoMatch_sp match = gc::As<ChemInfoMatch_sp>(match_mv.second());
      Atom_sp ax = match->tag(core::make_fixnum(1));
      Atom_sp ax1 = match->tag(core::make_fixnum(2));
      Atom_sp ax2 = match->tag(core::make_fixnum(3));
      Atom_sp ay = match->tag(core::make_fixnum(4));	// Carbonyl carbon
      Atom_sp ay1 = match->tag(core::make_fixnum(5));
      Atom_sp ay2 = match->tag(core::make_fixnum(6));
		    //
		    // If amide carbonyl carbon is in ring
		    // then we want a cis amide bond
		    // otherwise trans
		    //
		    // Achieve this by swapping ay1 and ay2 if its
		    // supposed to be a cis amide bond
		    //
      if ( ay->getRingMembershipCount() > 0 )
      {
        Atom_sp azz = ay1;
        ay1 = ay2;
        ay2 = azz;
      }
      LOG(BF("Applying a secondary amide restraint between %s and %s") % ax->description() % ay->description()  );
		    //
		    // H3(ax2) and O5(ay1) should be trans
      this->_applyDihedralRestraint(ax1,ax,ay,ay1,cisMin,cisMax,weight,activeAtoms);
      LOG(BF("Restrain cis %s - %s - %s -%s") % ax1->description() % ax->description() % ay->description() % ay1->description()  );
      this->_applyDihedralRestraint(ax1,ax,ay,ay2,transMin,transMax,weight,activeAtoms);
      LOG(BF("Restrain trans %s - %s - %s -%s") % ax1->description() % ax->description() % ay->description() % ay2->description()  );
      this->_applyDihedralRestraint(ax2,ax,ay,ay1,transMin,transMax,weight,activeAtoms);
      LOG(BF("Restrain trans %s - %s - %s -%s") % ax2->description() % ax->description() % ay->description() % ay1->description()  );
      this->_applyDihedralRestraint(ax2,ax,ay,ay2,cisMin,cisMax,weight,activeAtoms);
      LOG(BF("Restrain cis %s - %s - %s -%s") % ax2->description() % ax->description() % ay->description() % ay2->description()  );
    }
  }
}





SYMBOL_EXPORT_SC_(ChemPkg,find_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,force_fields_as_list);

SYMBOL_EXPORT_SC_(ChemPkg,add_shadowing_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,assign_force_field_types);
SYMBOL_EXPORT_SC_(ChemPkg,compute_merged_nonbond_force_field_for_aggregate);
SYMBOL_EXPORT_SC_(ChemPkg,force_field_component_merge);
SYMBOL_EXPORT_SC_(ChemPkg,clear_combined_force_field);
SYMBOL_EXPORT_SC_(ChemPkg,generate_molecule_energy_function_tables);
SYMBOL_EXPORT_SC_(ChemPkg,STARparameter_warningsSTAR);
SYMBOL_EXPORT_SC_(ChemPkg,report_parameter_warnings);

SYMBOL_EXPORT_SC_(ChemPkg,identify_aromatic_rings);
SYMBOL_EXPORT_SC_(ChemPkg,STARcurrent_aromaticity_informationSTAR);

SYMBOL_EXPORT_SC_(KeywordPkg,bounding_box);
CL_LISPIFY_NAME("defineForMatter");
CL_LAMBDA((energy-function !) matter &key use-excluded-atoms active-atoms (assign-types t));
CL_DEFMETHOD void EnergyFunction_O::defineForMatter(Matter_sp matter, bool useExcludedAtoms, core::T_sp activeAtoms, bool assign_types )
{_OF();

  if ( !(matter.isA<Aggregate_O>() || matter.isA<Molecule_O>() ) )
  {
    SIMPLE_ERROR(BF("You can only define energy functions for Aggregates or Molecules"));
  }

  core::DynamicScopeManager scope(_sym_STARparameter_warningsSTAR,_Nil<core::T_O>());
  
	//
	// Identify rings
	//
  if (chem__verbose(0)) core::write_bf_stream(BF("Searching for rings.\n"));

  core::List_sp rings = RingFinder_O::identifyRings(matter);
  core::DynamicScopeManager ring_scope(_sym_STARcurrent_ringsSTAR, rings );

  //
  // Calculate aromaticity using the rings we just calculated
  //
  core::T_sp aromaticity_info = core::eval::funcall(_sym_identify_aromatic_rings,rings);
  core::DynamicScopeManager aromaticity_scope(_sym_STARcurrent_aromaticity_informationSTAR,aromaticity_info);
  
  //
  // Assign relative Cahn-Ingold-Preylog priorities
  //
  if (chem__verbose(0)) core::write_bf_stream(BF("Assigning CIP priorities.\n"));
  CipPrioritizer_O::assignPriorities(matter);

    //
    // Assign atom types for each molecule
    //
  
  core::HashTableEq_sp force_fields = core::HashTableEq_O::create_default();
  Loop moleculeLoop;
  moleculeLoop.loopTopGoal(matter,MOLECULES);
  while (moleculeLoop.advanceLoopAndProcess() ) {
    Molecule_sp molecule = moleculeLoop.getMolecule();
    core::T_sp force_field_name = molecule->force_field_name();
    if (force_fields->gethash(force_field_name).nilp()) {
      core::T_sp combined_force_field = core::eval::funcall(_sym_find_force_field,force_field_name);
      force_fields->setf_gethash(force_field_name,combined_force_field);
    }
  }
  
  if (chem__verbose(0)) core::write_bf_stream(BF("Assigning atom types.\n"));
  moleculeLoop.loopTopGoal(matter,MOLECULES);
  while (moleculeLoop.advanceLoopAndProcess() ) {
    Molecule_sp molecule = moleculeLoop.getMolecule();
    core::T_sp force_field_name = molecule->force_field_name();
    core::T_sp combined_force_field = force_fields->gethash(force_field_name);
    if (chem__verbose(0)) core::write_bf_stream(BF("Assigning atom types for molecule %s using %s.\n") % _rep_(molecule->getName()) % _rep_(force_field_name));
    core::eval::funcall(_sym_assign_force_field_types,combined_force_field,molecule);
  }
  this->defineForMatterWithAtomTypes(matter,useExcludedAtoms,activeAtoms);
}


CL_LAMBDA((energy-function !) matter &key use-excluded-atoms active-atoms (assign-types t));
CL_DEFMETHOD void EnergyFunction_O::defineForMatterWithAtomTypes(Matter_sp matter, bool useExcludedAtoms, core::T_sp activeAtoms)
{_OF();

  if ( !(matter.isA<Aggregate_O>() || matter.isA<Molecule_O>() ) )
  {
    SIMPLE_ERROR(BF("You can only define energy functions for Aggregates or Molecules"));
  }

  //
  // Setup the extra information needed by the energy function to generate a topology file
  //
  // Get the name to pass that to the atomTable
  core::T_sp matterName = matter->getName();
  this->_AtomTable->setAggregateName(matterName);
  core::T_sp boundingBox = matter->getPropertyOrDefault(kw::_sym_bounding_box,_Unbound<core::T_O>());
  if (boundingBox.unboundp()) {
    this->_AtomTable->makUnboundBoundingBox();
  } else {
    this->_AtomTable->setBoundingBox(boundingBox);
  }
  if (chem__verbose(0)) core::write_bf_stream(BF("Assembling aggregate nonbond force-field.\n"));
  core::T_sp nonbondForceField = core::eval::funcall(chem::_sym_compute_merged_nonbond_force_field_for_aggregate,matter);
  this->_AtomTable->setNonbondForceFieldForAggregate(nonbondForceField);
  
	// 
	// Assign atom types
	//
  //
  // This should be user configurable from a dynamic variable that contains a bunch of
  // system parameters
  //
  this->_DielectricConstant = 80.0;

  // Separate the molecules for solute from the solvent and handle them solute first then solvent
  size_t final_solute_residue_iptres = 0;
  size_t number_of_molecules_nspm = 0;
  size_t first_solvent_molecule_nspsol = 0;
  bool solvent_exists = false;
  if (gc::IsA<Aggregate_sp>(matter)) {
    if (chem__verbose(0)) core::write_bf_stream(BF("Classifying solute/solvent.\n"));
    ql::list solute_molecules;
    ql::list solvent_molecules;
    {
      Loop moleculeLoop;
      moleculeLoop.loopTopGoal(matter,MOLECULES);
      while (moleculeLoop.advanceLoopAndProcess() ) {
        Molecule_sp molecule = moleculeLoop.getMolecule();
        if (molecule->molecule_type() == kw::_sym_solvent) {
          solvent_molecules << molecule;
          solvent_exists = true;
        } else {
          solute_molecules << molecule;
        }
      }
    }
    core::List_sp solute = solute_molecules.cons();
    for ( auto cur_solute : solute ) {
      Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solute));
      core::T_sp force_field_name = onemol->force_field_name();
      core::T_sp forceField = core::eval::funcall(chem::_sym_find_force_field,force_field_name);
      this->_AtomTable->constructFromMatter(onemol,nonbondForceField,activeAtoms);
      if (chem__verbose(0)) core::write_bf_stream(BF("Generating parameters for %s using %s force-field.\n") % _rep_(onemol->getName()) % _rep_(force_field_name) );
      core::eval::funcall(_sym_generate_molecule_energy_function_tables,this->asSmartPtr(),onemol,forceField,activeAtoms);
      final_solute_residue_iptres += onemol->contentSize();
      ++number_of_molecules_nspm;
    }    
    core::List_sp solvent = solvent_molecules.cons();
    first_solvent_molecule_nspsol = number_of_molecules_nspm+1;
    if (chem__verbose(0)) {
      if (core::cl__length(solvent)>0) {
        Molecule_sp onemol = gc::As<Molecule_sp>(oCar(oCar(solvent)));
        core::T_sp force_field_name = onemol->force_field_name();
        core::write_bf_stream(BF("Generating parameters for solvent %s using %s force-field.\n") % _rep_(onemol->getName()) % _rep_(force_field_name) );
      }
    }
    for ( auto cur_solvent : solvent ) {
      Molecule_sp onemol = gc::As_unsafe<Molecule_sp>(CONS_CAR(cur_solvent));
      core::T_sp force_field_name = onemol->force_field_name();
      core::T_sp forceField = core::eval::funcall(chem::_sym_find_force_field,force_field_name);
      this->_AtomTable->constructFromMatter(onemol,nonbondForceField,activeAtoms);
      core::eval::funcall(_sym_generate_molecule_energy_function_tables,this->asSmartPtr(),onemol,forceField,activeAtoms);
      ++number_of_molecules_nspm;
    }
  } else {
    Molecule_sp molecule = gc::As<Molecule_sp>(matter);
    core::T_sp force_field_name = molecule->force_field_name();
    core::T_sp forceField = core::eval::funcall(chem::_sym_find_force_field,force_field_name);
    this->_AtomTable->constructFromMatter(molecule,nonbondForceField,activeAtoms);
    if (chem__verbose(0)) core::write_bf_stream(BF("Generating parameters for %s using %s force-field.\n") % _rep_(molecule->getName()) % _rep_(force_field_name) );
    core::eval::funcall(_sym_generate_molecule_energy_function_tables,this->asSmartPtr(),molecule,forceField,activeAtoms);
    final_solute_residue_iptres = molecule->contentSize();
    number_of_molecules_nspm = 1;
    first_solvent_molecule_nspsol = 2;
  }
  if (solvent_exists) {
    printf("%s:%d  solvent_exists NSPSOL %lu  IPTRES %lu   NSM %lu \n", __FILE__, __LINE__,
           first_solvent_molecule_nspsol, final_solute_residue_iptres, number_of_molecules_nspm);
    this->_AtomTable->set_firstSolventMoleculeNSPSOL(first_solvent_molecule_nspsol);
    this->_AtomTable->set_finalSoluteResidueIPTRES(final_solute_residue_iptres);
    this->_AtomTable->set_totalNumberOfMoleculesNSPM(number_of_molecules_nspm);
  } else {
    printf("%s:%d  NO SOLVENT\n", __FILE__, __LINE__);
  }
  {
    core::T_sp nonbondForceField = this->_AtomTable->nonbondForceFieldForAggregate();
    this->generateNonbondEnergyFunctionTables(useExcludedAtoms,matter,nonbondForceField,activeAtoms);
    this->generateRestraintEnergyFunctionTables(matter,nonbondForceField,activeAtoms);
  }
  core::eval::funcall(_sym_report_parameter_warnings);
}







CL_LAMBDA((energy-function !) molecule force-field &key active-atoms);
CL_DOCSTRING("Generate the standard energy function tables. The atom types, and CIP priorities need to be precalculated.");
CL_DEFMETHOD void EnergyFunction_O::generateStandardEnergyFunctionTables(Matter_sp molecule, CombinedForceField_sp combinedForceField, core::T_sp activeAtoms )
{
  Loop loop;
  Atom_sp          a1, a2, a3, a4, aImproperCenter;
  core::Symbol_sp  t1, t2, t3, t4, t141, t144;
  EnergyAtom      *eaCenter, *ea1, *ea2, *ea3, *ea4;
  FFPtor_sp        ffPtor;
  FFItor_sp        ffItor;
  FFNonbond_sp	ffNonbond1, ffNonbond2;
  int             coordinateIndex;
  
  
  //
  // Define a Nonbond cross term table
  //
  // FIXME: I think this should only be done once after all of the molecules have been added
  // FIXME: Why are we using the forcefield->getNonbondDb()?  We should be using the aggregate nonbond force-field
#if 0  
  if (chem__verbose(0)) core::write_bf_stream(BF("Starting to build standard energy function tables\n"));
  GC_ALLOCATE(FFNonbondCrossTermTable_O, temp );
  this->_NonbondCrossTermTable = temp;
  this->_NonbondCrossTermTable->fillUsingFFNonbondDb(forceField->getNonbondDb());
#endif
  
    	//
	// Initialize the energy components
	//
  ALL_ENERGY_COMPONENTS(initialize());
//  this->_Matter= matter;
  this->_eraseMissingParameters();
  coordinateIndex = 0;
  ASSERTNOTNULL(forceField);
  // Merge the stretch terms
  FFStretchDb_sp ffstretches = FFStretchDb_O::create();
  core::List_sp parts = combinedForceField->forceFieldsAsList();
  for ( auto cur : parts ) {
    ForceField_sp other = gc::As<ForceField_sp>(CONS_CAR(cur));
    FFStretchDb_sp other_ffstretch = other->getStretchDb();
    ffstretches->forceFieldMerge(other_ffstretch);
  }
  // Search the stretch terms
  {_BLOCK_TRACE("Defining STRETCH");
    size_t terms = 0;
    size_t missing_terms = 0;
    loop.loopTopGoal(molecule,BONDS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getBondA1();
      a2 = loop.getBondA2();
      if ( activeAtoms.notnilp() &&
           (!inAtomSet(activeAtoms,a1) || !inAtomSet(activeAtoms,a2)) ) continue;
//      printf("%s:%d Looking at STRETCH term between %s - %s\n", __FILE__, __LINE__, _rep_(a1).c_str(), _rep_(a2).c_str());
//      t1 = a1->getType();
//      t2 = a2->getType();
      ea1 = this->getEnergyAtomPointer(a1);
      ea2 = this->getEnergyAtomPointer(a2);
      FFStretch_sp ffStretch = gc::As<FFStretch_sp>(ffstretches->findTerm(a1,a2));
      if ( ffStretch->level() != parameterized ) {
        this->_addMissingParameter(ffStretch);
        ++missing_terms;
      }
      if ( ffStretch->level() != unknown ) {
        EnergyStretch   energyStretch;
        energyStretch.defineFrom(ffStretch,ea1,ea2,this->_Stretch->getScale());
        this->_Stretch->addTerm(energyStretch);
        ++terms;
      } else {
        SIMPLE_WARN(BF("Could not find stretch parameter between %s and %s") % _rep_(a1) % _rep_(a2));
      }
    }
    if (chem__verbose(0)) core::write_bf_stream(BF("Built stretch table with %d terms added and %d missing terms\n") % terms % missing_terms);
  }
#ifdef	DEBUG_DEFINE_ENERGY
  lisp->print(BF("%s:%d There were %d stretch terms") % __FILE__ % __LINE__ % this->_Stretch.size() );
#endif
#if USE_ALL_ENERGY_COMPONENTS
  // Merge the angle terms
  FFAngleDb_sp ffangles = FFAngleDb_O::create();
  {
    core::List_sp parts = combinedForceField->forceFieldsAsList();
    for ( auto cur : parts ) {
      ForceField_sp other = gc::As<ForceField_sp>(CONS_CAR(cur));
      FFAngleDb_sp other_ffangle = other->getAngleDb();
      ffangles->forceFieldMerge(other_ffangle);
    }
  }
  // Search the angle terms
  {_BLOCK_TRACE("Defining ANGLES");
    size_t terms = 0;
    size_t missing_terms = 0;
    loop.loopTopGoal(molecule,ANGLES);
    while ( loop.advanceLoopAndProcess() )
    {
      a1 = loop.getAtom1();
      a2 = loop.getAtom2();
      a3 = loop.getAtom3();
      if ( activeAtoms.notnilp() &&
           ( !inAtomSet(activeAtoms,a1) || !inAtomSet(activeAtoms,a2) || !inAtomSet(activeAtoms,a3) )) continue;
      ea1 = this->getEnergyAtomPointer(a1);
      ea2 = this->getEnergyAtomPointer(a2);
      ea3 = this->getEnergyAtomPointer(a3);
      FFAngle_sp ffAngle = ffangles->findTerm(ffstretches,a1,a2,a3);
      if ( ffAngle->level() != parameterized ) {
        LOG(BF("Missing angle parameter between types: %s-%s-%s") % _rep_(a1->getType()) % _rep_(a2->getType()) % _rep_(a3->getType()) );
        this->_addMissingParameter(ffAngle);
        LOG(BF("Added to missing parameters") );
        ++missing_terms;
      }
      if ( ffAngle->level() != unknown ) {
        EnergyAngle energyAngle;
        energyAngle.defineFrom(ffAngle,ea1,ea2,ea3,this->_Angle->getScale());
        this->_Angle->addTerm(energyAngle);
        ++terms;
      }
    }
    if (chem__verbose(0)) core::write_bf_stream(BF("Built angle table with %d terms and %d missing terms\n") % terms % missing_terms);
  }
  // Merge the ptor terms
  FFPtorDb_sp ffptors = FFPtorDb_O::create();
  {
    core::List_sp parts = combinedForceField->forceFieldsAsList();
    for ( auto cur : parts ) {
      ForceField_sp other = gc::As<ForceField_sp>(CONS_CAR(cur));
      FFPtorDb_sp other_ffptor = other->getPtorDb();
      ffptors->forceFieldMerge(other_ffptor);
    }
  }
  // Search the ptor terms
  {_BLOCK_TRACE("Defining PROPERS");
    size_t terms = 0;
    size_t missing_terms = 0;
    loop.loopTopGoal(molecule,PROPERS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getAtom1();
      a2 = loop.getAtom2();
      a3 = loop.getAtom3();
      a4 = loop.getAtom4();
      if ( activeAtoms.notnilp() &&
           ( !inAtomSet(activeAtoms,a1)
             || !inAtomSet(activeAtoms,a2)
             || !inAtomSet(activeAtoms,a3)
             || !inAtomSet(activeAtoms,a4) )) continue;
      t1 = a1->getType();
      t2 = a2->getType();
      t3 = a3->getType();
      t4 = a4->getType();
//      core::write_bf_stream(BF("atoms types: %s-%s-%s-%s \n") % t1 % t2 % t3 % t4);
      ea1 = this->getEnergyAtomPointer(a1);
      ea2 = this->getEnergyAtomPointer(a2);
      ea3 = this->getEnergyAtomPointer(a3);
      ea4 = this->getEnergyAtomPointer(a4);
      core::T_sp tffPtor = ffptors->findBestTerm(t1,t2,t3,t4);
      if (tffPtor.notnilp()) {
        FFPtor_sp ffPtor = gc::As<FFPtor_sp>(tffPtor);
        int numPtors = 0;
        for ( int n=1;n<=FFPtor_O::MaxPeriodicity; n++ ) {
          if ( ffPtor->hasPeriodicity(n) ) {
            ++numPtors;
//            LOG(BF( "Adding proper term for atoms %s-%s-%s-%s types: %s-%s-%s-%s")%
//                ea1->getResidueAndName()
//                % ea2->getResidueAndName()
//                % ea3->getResidueAndName()
//                % ea4->getResidueAndName()
//                % t1 % t2 % t3 % t4
//                );
#if 0
            core::write_bf_stream(BF( "Adding proper term for atoms %s-%s-%s-%s types: %s-%s-%s-%s /n")%
                                  ea1->getResidueAndName()
                                  % ea2->getResidueAndName()
                                  % ea3->getResidueAndName()
                                  % ea4->getResidueAndName()
                                  % t1 % t2 % t3 % t4
                                  );
#endif
            EnergyDihedral energyDihedral;
            energyDihedral.defineFrom(n,ffPtor,ea1,ea2,ea3,ea4,this->_Dihedral->getScale());
            this->_Dihedral->addTerm(energyDihedral);
            ++terms;
          }
        }
        if ( numPtors == 0 ) {
          FFPtor_sp ptor = FFPtor_O::create_missing(t1,t2,t3,t4);
          this->_addMissingParameter(ptor);
          ++missing_terms;
        }
        if (ea1->inBondOrAngle(ea4->atom()) )
        {
#if 0
#ifdef	DEBUG_ON
          if ( _rep_(t1) < _rep_(t4) ) {
            t141 = t1;
            t144 = t4;
          }else{
            t141 = t4;
            t144 = t1;
          }
//          LOG(BF("Defining 1-4 interaction %-9s- %-9s   ") % t1 % t4 );
//          core::write_bf_stream(BF("Defining 1-4 interaction %-9s- %-9s   \n") % t1 % t4 );
#endif
#endif
        } else {
#if 0          
#ifdef	DEBUG_ON
          if ( _rep_(t1) < _rep_(t4) ) {
            t141 = t1;
            t144 = t4;
          }else{
            t141 = t4;
            t144 = t1;
          }
//          LOG(BF("Ignoring 1-4 interaction %-9s- %-9s    ") % t1 % t4 );
//          core::write_bf_stream(BF("Ignoring 1-4 interaction %-9s- %-9s    \n") % t1 % t4 );
#endif
#endif
        }
//		ea1->_CloserThan15.insert(ea4->_Atom);
//		ea4->_CloserThan15.insert(ea1->_Atom);
      }
    }
//    if (chem__verbose(0)) core::write_bf_stream(BF("Built dihedral table with %d terms and %d missing terms\n") % terms % missing_terms);
//   core::write_bf_stream(BF("Built dihedral table with %d terms and %d missing terms\n") % terms % missing_terms);
  }
  // Merge the itor terms
  FFItorDb_sp ffitors = FFItorDb_O::create();
  {
    core::List_sp parts = combinedForceField->forceFieldsAsList();
    for ( auto cur : parts ) {
      ForceField_sp other = gc::As<ForceField_sp>(CONS_CAR(cur));
      FFItorDb_sp other_ffitor = other->getItorDb();
      ffitors->forceFieldMerge(other_ffitor);
    }
  }
  // Search the itor terms
  {_BLOCK_TRACE("Defining IMPROPERS");
    EnergyDihedral energyDihedral;
    size_t terms = 0;
    loop.loopTopGoal(molecule,IMPROPERS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getAtom1();
      a2 = loop.getAtom2();
      a3 = loop.getAtom3();
      a4 = loop.getAtom4();
      FFItorDb_O::improperAtomSort(a1,a2,a4);
      if ( activeAtoms.notnilp() &&
           (!inAtomSet(activeAtoms,a1)
            || !inAtomSet(activeAtoms,a2)
            || !inAtomSet(activeAtoms,a3)
            || !inAtomSet(activeAtoms,a4)) ) continue;
      t1 = a1->getType();
      t2 = a2->getType();
      t3 = a3->getType();
      t4 = a4->getType();
      ea1 = this->getEnergyAtomPointer(a1);
      ea2 = this->getEnergyAtomPointer(a2);
      ea3 = this->getEnergyAtomPointer(a3);
      ea4 = this->getEnergyAtomPointer(a4);
      core::T_sp itor = ffitors->findBestTerm(t1,t2,t3,t4);
      if ( itor.notnilp() ) {
		    //
		    // Only one improper per central atom.
		    // We may not get exactly the same improper as AMBER does
		    //
        if ( a3 != aImproperCenter ) {
          ffItor = gc::As<FFItor_sp>(ffitors->findBestTerm(t1,t2,t3,t4));
          for ( int n=1;n<=6; n++ ) {
            if ( ffItor->hasPeriodicity(n) ) {
              energyDihedral.defineFrom(n,ffItor,ea1,ea2,ea3,ea4,this->_Dihedral->getScale());
              this->_Dihedral->addTerm(energyDihedral);
              aImproperCenter = a3;
              ++terms;
            }
          }
        }
      }
    }
    if (chem__verbose(0)) core::write_bf_stream(BF("Built improper table for %d terms\n") % terms);
  }
//  this->summarizeTerms();
}

SYMBOL_EXPORT_SC_(ChemPkg,prepare_amber_energy_nonbond);

CL_DOCSTRING("Generate the nonbond energy function tables. The atom types, and CIP priorities need to be precalculated.");
CL_DEFMETHOD void EnergyFunction_O::generateNonbondEnergyFunctionTables(bool useExcludedAtoms, Matter_sp matter, core::T_sp nonbondForceField, core::T_sp activeAtoms )
{_OF();
  if (chem__verbose(0))
    core::write_bf_stream(BF("Built atom table for %d atoms\n") % this->_AtomTable->getNumberOfAtoms());
#ifdef	DEBUG_DEFINE_ENERGY
  lisp->print(BF("%s:%d There were %d atoms") % __FILE__ % __LINE__ % this->_AtomTable.size() );
#endif
        // Nonbonds here!!!!!!!!!!!!!!
  if (useExcludedAtoms) {
// The nonbond parameters are calculated in Common Lisp
//    printf("atom number: %d\n", this->_AtomTable->getNumberOfAtoms());

//    printf("%s:%d in generateNonbondEnergyFunctionTables\n", __FILE__, __LINE__);
    core::List_sp parts = core::eval::funcall(_sym_prepare_amber_energy_nonbond,this->asSmartPtr(),nonbondForceField);
//      printf("%s:%d:%s    parts -> %s\n", __FILE__, __LINE__, __FUNCTION__, _rep_(parts).c_str());
//    printf("%s:%d in generateNonbondEnergyFunctionTables\n", __FILE__, __LINE__);
    this->_Nonbond->constructNonbondTermsFromAList(parts);
    // -----------
//    printf("%s:%d in generateNonbondEnergyFunctionTables\n", __FILE__, __LINE__);
      
    this->_Nonbond->constructExcludedAtomListFromAtomTable(this->_AtomTable, nonbondForceField);
    this->_Nonbond->construct14InteractionTerms(this->_AtomTable,matter,nonbondForceField,activeAtoms);
    LOG(BF("Done construct14InteractionTerms"));
//    printf("%s:%d:%s    nonbond -> %d\n", __FILE__, __LINE__, __FUNCTION__, _Nonbond->numberOfTerms());

  } else {
    this->_Nonbond->constructNonbondTermsFromAtomTable(false,this->_AtomTable, nonbondForceField);
  }
  if (chem__verbose(0)) core::write_bf_stream(BF("Built nonbond table for %d terms\n") % this->_Nonbond->numberOfTerms());
}


CL_LAMBDA((energy-function !) matter force-field &key active-atoms);
CL_DOCSTRING(R"doc(Generate the restraint energy function tables. The atom types, and CIP priorities need to be precalculated.
This should be called after generateStandardEnergyFunctionTables.)doc");
CL_DEFMETHOD void EnergyFunction_O::generateRestraintEnergyFunctionTables(Matter_sp matter, core::T_sp ffNonbond, core::T_sp activeAtoms )
{
  Loop loop;
  Atom_sp          a1, a2, a3, a4, aImproperCenter;
  core::Symbol_sp  t1, t2, t3, t4, t141, t144;
  EnergyAtom      *eaCenter, *ea1, *ea2, *ea3, *ea4;
  FFAngle_sp       ffAngle;
  FFPtor_sp        ffPtor;
  FFItor_sp        ffItor;
  FFNonbond_sp	ffNonbond1, ffNonbond2;
  int             coordinateIndex;
    	//
	// Setup the atom chiral restraints
	//
  {_BLOCK_TRACE("Defining chiral restraints");
    EnergyChiralRestraint	ichiral;
    Atom_sp	n1,n2,n3,n4;
    string	s1,s2,s3,s4;
    double	side ;
    loop.loopTopGoal(matter,ATOMS);
    while ( loop.advanceLoopAndProcess() ) {
      a1 = loop.getAtom();
      if ( activeAtoms.notnilp() && !inAtomSet(activeAtoms,a1) ) continue;
      if ( a1->getStereochemistryType() != undefinedCenter ) 
      {
        LOG(BF("Create a chiral restraint for %s") % a1->description()  );
			//
			// Figure out what the desired configuration should be
			// If it has been set then use that
			// Otherwise if its pro-chiral center then set it to "R"
			// If its a chiral center and it hasn't been set then throw an exception
			//
			// Configuration is set in variable side
			// "R" -> side = 1.0
			// "S" -> side = -1.0
			//
			// Calculate the angle for the chiral restraint based on
			// the stereochemistry that we want to impose
			// default R-stereochemistry (1-center)x(2-center).(3-center) is POSITIVE
			//
        if ( a1->getConfiguration() != undefinedConfiguration )
        {
          if ( a1->getConfiguration() == R_Configuration )
          {
            side = 1.0;
          } else
          {
            side = -1.0;
          }
        } else
        {
          if ( a1->getStereochemistryType() == prochiralCenter )
          {
            side = 1.0;
          } else
          {
            SIMPLE_ERROR(BF("Chiral center (%s) with configuration settings[%s] doesn't have its configuration set")
                         % a1->description()
                         % a1->getConfigurationAsString() );
          }
        }
        core::List_sp priority = a1->getNeighborsByRelativePriority();
        ASSERTP(core::cl__length(priority) == 4, "There must be 4 neighbors to assign stereochemistry");
        core::List_sp cur = priority;
        n1 = cur.asCons()->car<Atom_O>();
        cur = cur.asCons()->cdr();
        n2 = cur.asCons()->car<Atom_O>();
        cur = cur.asCons()->cdr();
        n3 = cur.asCons()->car<Atom_O>();
        cur = cur.asCons()->cdr();
        n4 = cur.asCons()->car<Atom_O>();
#if 0
        s1 = a1->getConfigurationPriorityHighest();
        s2 = a1->getConfigurationPriorityHigh();
        s3 = a1->getConfigurationPriorityLow();
        s4 = a1->getConfigurationPriorityLowest();
        n1 = a1->bondedNeighborWithName(s1);
        n2 = a1->bondedNeighborWithName(s2);
        n3 = a1->bondedNeighborWithName(s3);
        n4 = a1->bondedNeighborWithName(s4);
        ASSERTNOTNULLP(n1, "Atom("+a1->getName()+") does not have neighbor1("+s1+")");
        ASSERTNOTNULLP(n2, "Atom("+a2->getName()+") does not have neighbor2("+s2+")");
        ASSERTNOTNULLP(n3, "Atom("+a3->getName()+") does not have neighbor3("+s3+")");
        ASSERTNOTNULLP(n4, "Atom("+a4->getName()+") does not have neighbor4("+s4+")");
#endif
        eaCenter = this->getEnergyAtomPointer(a1);
        ea1 = this->getEnergyAtomPointer(n1);
        ea2 = this->getEnergyAtomPointer(n2);
        ea3 = this->getEnergyAtomPointer(n3);
        ea4 = this->getEnergyAtomPointer(n4);

			//
			// Setup chiral restraints for 1->2->center->3
			//			and 1->2->center->4
			//
        ichiral._Atom1 = ea1->atom();
        ichiral._Atom2 = ea2->atom();
        ichiral._Atom3 = eaCenter->atom();
        ichiral._Atom4 = ea3->atom();
        ichiral.term.I1 = ea1->coordinateIndexTimes3();
        ichiral.term.I2 = ea2->coordinateIndexTimes3();
        ichiral.term.I3 = eaCenter->coordinateIndexTimes3();
        ichiral.term.I4 = ea3->coordinateIndexTimes3();
        ichiral.term.K = this->_ChiralRestraintWeight * side;
        ichiral.term.CO = this->_ChiralRestraintOffset;
        this->_ChiralRestraint->addTerm(ichiral);
				// Now apply it to the other atom
				// on the chiral center, just flip the sign
				// of K
        ichiral._Atom4 = ea4->atom();
        ichiral.term.I4 = ea4->coordinateIndexTimes3();
					// flip the sign of the chiral restraint
        ichiral.term.K = this->_ChiralRestraintWeight * side * -1.0;
        ichiral.term.CO = this->_ChiralRestraintOffset;
        this->_ChiralRestraint->addTerm(ichiral);

			// To try and increase the number of molecules that
			// minimize into the correct configuration I'll add another
			// set of chiral restraints for this atom
			//
			// Setup chiral restraints for 2->4->center->3
			//			and 2->4->center->1
			//
        ichiral._Atom1 = ea2->atom();
        ichiral._Atom2 = ea4->atom();
        ichiral._Atom3 = eaCenter->atom();
        ichiral._Atom4 = ea3->atom();
        ichiral.term.I1 = ea2->coordinateIndexTimes3();
        ichiral.term.I2 = ea4->coordinateIndexTimes3();
        ichiral.term.I3 = eaCenter->coordinateIndexTimes3();
        ichiral.term.I4 = ea3->coordinateIndexTimes3();
        ichiral.term.K = this->_ChiralRestraintWeight * side;
        ichiral.term.CO = this->_ChiralRestraintOffset;
        this->_ChiralRestraint->addTerm(ichiral);
				// Now apply it to the other atom
				// on the chiral center, just flip the sign
				// of K
        ichiral._Atom4 = ea1->atom();
        ichiral.term.I4 = ea1->coordinateIndexTimes3();
					// flip the sign of the chiral restraint
        ichiral.term.K = this->_ChiralRestraintWeight * side * -1.0;
        ichiral.term.CO = this->_ChiralRestraintOffset;
        this->_ChiralRestraint->addTerm(ichiral);
      } else
      {
        LOG(BF("There is no chiral restraint for: %s") % a1->description()  );
      }
    }
    if (chem__verbose(0)) core::write_bf_stream(BF("Built chiral restraints table for %d terms\n") % this->_ChiralRestraint->numberOfTerms());
  }

	//
	// Set up the anchor restraints
	//
#if ATOMIC_ANCHOR
  {_BLOCK_TRACE("Defining anchor restraints");
    if ( this->_AtomTable->getNumberOfAtoms() > 0 )
    {
      EnergyAnchorRestraint	iterm(_lisp);
      Atom_sp	n1;
      string	s1;
      Vector3	anchorPos;
      gctools::Vec0<EnergyAtom>::iterator	iea1;
      for ( iea1 = this->_AtomTable->begin();
            iea1 != this->_AtomTable->end()-1; iea1++ )
      {
        a1 = iea1->atom();
        if ( activeAtoms.notnilp() && !inAtomSet(activeAtoms,a1) ) continue;
        if ( a1->isAnchorRestraintOn() )
        {
          LOG(BF("Create an anchor restraint for %s") % a1->description()  );
          anchorPos = a1->getAnchorPos();
          iterm.term.xa = anchorPos.getX();
          iterm.term.ya = anchorPos.getY();
          iterm.term.za = anchorPos.getZ();
          iterm.term.ka = this->_AnchorRestraintWeight;
          iterm.term.I1 = iea1->coordinateIndexTimes3();
          this->_AnchorRestraint->addTerm(iterm);
        }
      }
    } else
    {
      LOG(BF("There are no atoms"));
    }
    if (chem__verbose(0)) core::write_bf_stream(BF("Built anchor restraints table for %d terms\n") % this->_AnchorRestraint->numberOfTerms());
  }
#endif
	//
	// Define secondary amide restraints
	//
  if ( this->_RestrainSecondaryAmides )
  {_BLOCK_TRACE("Secondary amide restraints");
    gctools::Vec0<Atom_sp>	nitrogens;
    for ( AtomTable_O::iterator it=this->_AtomTable->begin(); it!=this->_AtomTable->end(); it++ )
    {
      Atom_sp a = it->atom();
      if ( activeAtoms.notnilp() && !inAtomSet(activeAtoms,a1) ) continue;
      if ( a->getElement() == element_N )
      {
        nitrogens.push_back(a);
      }
    }
    int startTerms = this->_ImproperRestraint->numberOfTerms();
    this->__createSecondaryAmideRestraints(nitrogens,activeAtoms);
    if (chem__verbose(0)) core::write_bf_stream(BF("Built secondary amide restraints including %d terms\n") % (this->_ImproperRestraint->numberOfTerms() - startTerms));
  } else
  {
    LOG(BF("Skipping Secondary amide restraints because _RestrainSecondaryAmides = %d") % this->_RestrainSecondaryAmides );
  }

#endif // USE_ALL_ENERGY_COMPONENTS

	//
	// Set up force-field restraints
	//
  {_BLOCK_TRACE("Defining force-field restraints");
    IterateRestraints_sp restraintIt = IterateRestraints_O::create(matter);
    int terms = this->_applyRestraints(ffNonbond,restraintIt,activeAtoms);
    if (chem__verbose(0)) core::write_bf_stream(BF("Built restraints including %d terms\n") % terms );
  }
  LOG(BF("Done terms") );
}



CL_LISPIFY_NAME("addTermsForListOfRestraints");
CL_DEFMETHOD void	EnergyFunction_O::addTermsForListOfRestraints(ForceField_sp forceField, core::List_sp restraintList, core::T_sp activeAtoms)
{
  adapt::IterateCons_sp	iterate;
  iterate = adapt::IterateCons_O::create(restraintList);
  this->_applyRestraints(forceField->getNonbondDb(),iterate,activeAtoms);
}

void	EnergyFunction_O::loadCoordinatesIntoVector(NVector_sp pos)
{
  int                             ci;
  gctools::Vec0<EnergyAtom>::iterator    ai;
  if ( pos->size() != this->getNVectorSize()) {
    SIMPLE_ERROR(BF("NVector is the incorrect length"));
  }
  for ( ai=this->_AtomTable->begin(); ai!=this->_AtomTable->end(); ai++ ) {
    ci = ai->coordinateIndexTimes3();
    pos->setElement(ci,ai->atom()->getPosition().getX());
    pos->setElement(ci+1, ai->atom()->getPosition().getY());
    pos->setElement(ci+2, ai->atom()->getPosition().getZ());
  }
}



void    EnergyFunction_O::saveCoordinatesFromVector(NVector_sp pos)
{
  int                             ci;
  double                          x,y,z;
  gctools::Vec0<EnergyAtom>::iterator    ai;
  Vector3                         v;
  for ( ai=this->_AtomTable->begin(); ai!=this->_AtomTable->end(); ai++ ) {
    ci = ai->coordinateIndexTimes3();
    x = pos->getElement(ci+0);
    y = pos->getElement(ci+1);
    z = pos->getElement(ci+2);
    LOG(BF("Set atom(%d) position (%lf,%lf,%lf)") % ci % x % y % z  );
    v.set(x,y,z);
    ai->atom()->setPosition(v);
  }
}

CL_LISPIFY_NAME("writeForceToAtoms");
CL_DEFMETHOD void    EnergyFunction_O::writeForceToAtoms(NVector_sp force)
{
  IMPLEMENT_ME();
  int                             ci;
  double                          x,y,z;
  gctools::Vec0<EnergyAtom>::iterator    ai;
  Vector3                         v;
  for ( ai=this->_AtomTable->begin(); ai!=this->_AtomTable->end(); ai++ ) {
    ci = ai->coordinateIndexTimes3();
    x = force->getElement(ci+0);
    y = force->getElement(ci+1);
    z = force->getElement(ci+2);
    v.set(x,y,z);
//    ai->atom()->setForce(v);
  }
}


void    EnergyFunction_O::saveCoordinatesAndForcesFromVectors(NVector_sp pos, NVector_sp force)
{
  this->saveCoordinatesFromVector(pos);
//  this->writeForceToAtoms(force);
}


EnergyAtom*	EnergyFunction_O::getEnergyAtomPointer(Atom_sp a)
{ 
  return this->_AtomTable->getEnergyAtomPointer(a);
};





void	EnergyFunction_O::dealWithProblem(core::Symbol_sp problem, core::T_sp error_args)
{_OF();
  core::List_sp atoms = _Nil<core::T_O>();
  if ( error_args.consp() && CONS_CAR(error_args) == kw::_sym_atoms ) {
    atoms = atoms;
  }
  for ( auto cur : atoms ) {
    Atom_sp a = gctools::As<Atom_sp>(oCar(cur));
    a->bumpPosition(0.1);
  }
}



CL_LISPIFY_NAME("summarizeBeyondThresholdInteractionsAsString");
CL_DEFMETHOD string EnergyFunction_O::summarizeBeyondThresholdInteractionsAsString()
{
  stringstream	ss;
  ss.str("");
  ss << "Beyond threshold ";
  ss << this->_Stretch->beyondThresholdInteractionsAsString();
#if USE_ALL_ENERGY_COMPONENTS
  ss << this->_Angle->beyondThresholdInteractionsAsString();
  ss << this->_Dihedral->beyondThresholdInteractionsAsString();
  ss << this->_Nonbond->beyondThresholdInteractionsAsString();
  ss << this->_ImproperRestraint->beyondThresholdInteractionsAsString();
  ss << this->_ChiralRestraint->beyondThresholdInteractionsAsString();
  ss << this->_AnchorRestraint->beyondThresholdInteractionsAsString();
  ss << this->_FixedNonbondRestraint->beyondThresholdInteractionsAsString();
#endif
  return ss.str();
}


CL_LISPIFY_NAME("summarizeEnergyAsString");
CL_DEFMETHOD string	EnergyFunction_O::summarizeEnergyAsString()
{
  stringstream	ss;
  ss.str("");
  ss << this->_Stretch->summarizeEnergyAsString();
#if USE_ALL_ENERGY_COMPONENTS
  ss << this->_Angle->summarizeEnergyAsString();
  ss << this->_Dihedral->summarizeEnergyAsString();
  ss << this->_Nonbond->summarizeEnergyAsString();
  ss << this->_ImproperRestraint->summarizeEnergyAsString();
  ss << this->_ChiralRestraint->summarizeEnergyAsString();
  ss << this->_AnchorRestraint->summarizeEnergyAsString();
  ss << this->_FixedNonbondRestraint->summarizeEnergyAsString();
#endif
  ss<< "             Total energy: " << this->_TotalEnergy << std::endl;
  return ss.str();
}



CL_LISPIFY_NAME("hasMissingParameters");
CL_DEFMETHOD bool	EnergyFunction_O::hasMissingParameters()
{
  return this->_MissingParameters.notnilp();
}


CL_LISPIFY_NAME("getMissingParameters");
CL_DEFMETHOD core::List_sp EnergyFunction_O::getMissingParameters()
{
  return this->_MissingParameters;
}

CL_LISPIFY_NAME("debugLogAsString");
CL_DEFMETHOD string	EnergyFunction_O::debugLogAsString()
{
  stringstream ss;
  ss << this->_Stretch->debugLogAsString() << std::endl;
#if USE_ALL_ENERGY_COMPONENTS
  ss << this->_Angle->debugLogAsString() << std::endl;
  ss << this->_Dihedral->debugLogAsString() << std::endl;
  ss << this->_Nonbond->debugLogAsString() << std::endl;
  ss << this->_ImproperRestraint->debugLogAsString() << std::endl;
  ss << this->_ChiralRestraint->debugLogAsString() << std::endl;
  ss << this->_AnchorRestraint->debugLogAsString() << std::endl;
  ss << this->_FixedNonbondRestraint->debugLogAsString() << std::endl;
#endif
  return ss.str();
}

SYMBOL_EXPORT_SC_(KeywordPkg,atom_table);
SYMBOL_EXPORT_SC_(KeywordPkg,stretch);
SYMBOL_EXPORT_SC_(KeywordPkg,angle);
SYMBOL_EXPORT_SC_(KeywordPkg,dihedral);
SYMBOL_EXPORT_SC_(KeywordPkg,nonbond);

CL_DEFUN void chem__fill_energy_function_from_alist(EnergyFunction_sp energy, core::List_sp alist)
{
  energy->_AtomTable = (safe_alist_lookup<AtomTable_sp>(alist,kw::_sym_atom_table));
  energy->_Stretch = (safe_alist_lookup<EnergyStretch_sp>(alist,kw::_sym_stretch));
  energy->_Angle = (safe_alist_lookup<EnergyAngle_sp>(alist,kw::_sym_angle));
  energy->_Dihedral = (safe_alist_lookup<EnergyDihedral_sp>(alist,kw::_sym_dihedral));
  energy->_Nonbond = (safe_alist_lookup<EnergyNonbond_sp>(alist,kw::_sym_nonbond));
}


void EnergyFunction_O::enableDebug()
{
  ALL_ENERGY_COMPONENTS(enableDebug());
}

void EnergyFunction_O::disableDebug()
{
  ALL_ENERGY_COMPONENTS(disableDebug());
}





};
