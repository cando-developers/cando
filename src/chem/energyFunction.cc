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

#include "core/common.h"
#include "energyFunction.h"
#include "loop.h"
#include "core/environment.h"
#include "minimizerLog.h"
#include "restraint.h"
#include "iterateRestraints.h"
#include "addon/iterateCons.h"
#include "ringFinder.h"
#include "cipPrioritizer.h"
#include "atom.h"
#include "virtualAtom.h"
#include "energyAtomTable.h"
#include "energyStretch.h"
#if USE_ALL_ENERGY_COMPONENTS
#include "energyAngle.h"
#include "energyDihedral.h"
#include "energyNonbond.h"
#include "energyImproperRestraint.h"
#include "energyChiralRestraint.h"
#include "energyAnchorRestraint.h"
#include "energyFixedNonbond.h"
#endif
#include "symbolTable.h"
#include "ffBaseDb.h"
#include "ffTypesDb.h"
#include "ffStretchDb.h"
#include "ffAngleDb.h"
#include "ffItorDb.h"
#include "ffPtorDb.h"
#include "ffNonbondDb.h"
#include "ffVdwDb.h"
#include "forceField.h"
#include "largeSquareMatrix.h"
#include "chemInfo.h"
#include "core/wrappers.h"


#if 0
#define	TURN_DEBUG_OFF
#define	MATHEMATICA_DEBUG_ON	1
#define	PARMDEBUG	1
#endif

#include "core/profiler.h"

#ifdef	darwin
#include "xmmintrin.h"
#endif





namespace chem
{

    SYMBOL_EXPORT_SC_(ChemPkg,nonbondTerm);
    SYMBOL_EXPORT_SC_(ChemPkg,restraintAnchor);



bool energyFunctionInitialized = false;
ChemInfo_sp secondaryAmide;


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

void energyFunction_initializeSmarts(core::Lisp_sp lisp)
{_G();
    if ( energyFunctionInitialized ) return;
    energyFunctionInitialized = true;
    secondaryAmide = ChemInfo_O::create();
    bool ok = secondaryAmide->compileSmarts("[$(N1(C2)(~[#1]3)~C4(=O5)C6)]");
    if ( !ok )
    {
	SIMPLE_ERROR(BF("Error compiling secondary amide MSMARTS"));
    }
}



#if INIT_TO_FACTORIES

#define ARGS_EnergyFunction_O_make "(matter force_field)"
#define DECL_EnergyFunction_O_make ""
#define DOCS_EnergyFunction_O_make "make EnergyFunction"
EnergyFunction_sp EnergyFunction_O::make(Matter_sp matter, ForceField_sp forceField)
  {_G();
      GC_ALLOCATE(EnergyFunction_O, me );
    if ( matter.notnilp() )
      {
	if ( forceField.notnilp() )
	  {
	    me->defineForMatter(matter,forceField);
	  } else
	  {
	      SIMPLE_ERROR(BF("You must provide a forceField if you provide a matter object"));
	  }
      }
    return me;
  };

#else

core::T_sp EnergyFunction_O::__init__(core::Function_sp exec, core::Cons_sp args,
				     core::Environment_sp env, core::Lisp_sp lisp)
{_OF();
    Matter_sp matter = env->lookup(_lisp->internWithPackageName(ChemPkg,"matter")).as<Matter_O>();
    ForceField_sp forceField = env->lookup(_lisp->internWithPackageName(ChemPkg,"forceField")).as<ForceField_O>();
    if ( matter.notnilp() )
    {
	if ( forceField.notnilp() )
	{
	    this->defineForMatter(matter,forceField);
	} else
	{
	    SIMPLE_ERROR(BF("You must provide a forceField if you provide a matter object"));
	}
    }
    return _Nil<core::T_O>();
}

#endif

    void EnergyFunction_O::exposeCando(core::Lisp_sp lisp)
    {
	Defun_maker(ChemPkg,EnergyFunction);
	core::class_<EnergyFunction_O>()
//	    .def_raw("core:__init__",&EnergyFunction_O::__init__,"(self &key matter forceField)")
    	.def("setOptions",&EnergyFunction_O::setOptions)
    	.def("setOption",&EnergyFunction_O::setOption)
	    .def("enableDebug",&EnergyFunction_O::enableDebug,"","","",false)
	    .def("disableDebug",&EnergyFunction_O::disableDebug,"","","",false)
    	.def("debugLogAsString",&EnergyFunction_O::debugLogAsString)
    	.def("energyComponentsAsString",&EnergyFunction_O::energyComponentsAsString)
	.def("defineForMatter",&EnergyFunction_O::defineForMatter)
	.def("useDefaultSettings",&EnergyFunction_O::useDefaultSettings)
	.def("checkForBeyondThresholdInteractions",&EnergyFunction_O::checkForBeyondThresholdInteractions)
	.def("addTermsForListOfRestraints",&EnergyFunction_O::addTermsForListOfRestraints)
	    .def("setName",&EnergyFunction_O::setName,"","","",false)
	    .def("getName",&EnergyFunction_O::getName,"","","",false)
	    .def("atomTable",&EnergyFunction_O::atomTable)
	.def("getStretchComponent",&EnergyFunction_O::getStretchComponent)
	.def("getAngleComponent",&EnergyFunction_O::getAngleComponent)
	.def("getDihedralComponent",&EnergyFunction_O::getDihedralComponent)
	.def("getNonbondComponent",&EnergyFunction_O::getNonbondComponent)
	.def("getImproperRestraintComponent",&EnergyFunction_O::getImproperRestraintComponent)
	.def("getChiralRestraintComponent",&EnergyFunction_O::getChiralRestraintComponent)
	.def("getAnchorRestraintComponent",&EnergyFunction_O::getAnchorRestraintComponent)
	.def("getFixedNonbondRestraintComponent",&EnergyFunction_O::getFixedNonbondRestraintComponent)
//	.def("getImproperRestraintComponentEnergy",&EnergyFunction_O::getImproperRestraintComponentEnergy)
//	.def("getDihedralComponentEnergy",&EnergyFunction_O::getDihedralComponentEnergy)
	.def("getMessage",&EnergyFunction_O::getMessage)
	.def("getMatter",&EnergyFunction_O::getMatter)
	.def("getNVectorSize",&EnergyFunction_O::getNVectorSize)
	.def("evaluateAll",&EnergyFunction_O::evaluateAll)
	.def("evaluateEnergy",&EnergyFunction_O::evaluateEnergy)
	.def("evaluateEnergyForce",&EnergyFunction_O::evaluateEnergyForce)
	.def("evaluateEnergyForceFullHessian",&EnergyFunction_O::evaluateEnergyForceFullHessian)
	.def("evaluateEnergyForceFullHessianForDebugging",&EnergyFunction_O::evaluateEnergyForceFullHessianForDebugging)
	.def("evaluateNumericalForce",&EnergyFunction_O::evaluateNumericalForce)
	.def("evaluateNumericalHessian",&EnergyFunction_O::evaluateNumericalHessian)
	.def("calculateEnergy",&EnergyFunction_O::calculateEnergy)
	.def("calculateEnergyAndForce",&EnergyFunction_O::calculateEnergyAndForce)
	.def("extractCoordinatesFromAtoms",&EnergyFunction_O::extractCoordinatesFromAtoms)
	.def("writeCoordinatesToAtoms",&EnergyFunction_O::writeCoordinatesToAtoms)
	.def("writeCoordinatesAndForceToAtoms",&EnergyFunction_O::writeCoordinatesAndForceToAtoms)
	.def("writeForceToAtoms",&EnergyFunction_O::writeForceToAtoms)
	.def("getTotalEnergy",&EnergyFunction_O::getTotalEnergy)
	.def("dumpTerms",&EnergyFunction_O::dumpTerms)
	.def("calculateNumericalDerivative",&EnergyFunction_O::calculateNumericalDerivative)
	.def("calculateNumericalSecondDerivative",
		&EnergyFunction_O::calculateNumericalSecondDerivative)
	.def("energyTermsEnabled",&EnergyFunction_O::energyTermsEnabled)
	.def("summarizeBeyondThresholdInteractionsAsString",
		&EnergyFunction_O::summarizeBeyondThresholdInteractionsAsString)
	.def("countTermsBeyondThreshold",&EnergyFunction_O::countTermsBeyondThreshold)
	.def("setupHessianPreconditioner",&EnergyFunction_O::setupHessianPreconditioner)
	.def("getMissingParameters",&EnergyFunction_O::getMissingParameters)
//	.def("countBadVdwInteractions",&EnergyFunction_O::countBadVdwInteractions)
	.def("hasMissingParameters",&EnergyFunction_O::hasMissingParameters)
	.def("summarizeEnergyAsString",&EnergyFunction_O::summarizeEnergyAsString)
	.def("compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition",
		&EnergyFunction_O::compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition)
	;
//	defInPackage(ChemPkg,"create_EnergyFunction", &EnergyFunction_O::create);
}

void EnergyFunction_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef USEBOOSTPYTHON
    PYTHON_CLASS(ChemPkg,EnergyFunction,"","",_lisp)
    	.def("setOptions",&EnergyFunction_O::setOptions)
    	.def("setOption",&EnergyFunction_O::setOption)
//    	.def("enableDebug",&EnergyFunction_O::enableDebug,"","","",false)
//    	.def("disableDebug",&EnergyFunction_O::disableDebug,"","","",false)
    	.def("debugLogAsString",&EnergyFunction_O::debugLogAsString)
    	.def("energyComponentsAsString",&EnergyFunction_O::energyComponentsAsString)
	.def("defineForMatter",&EnergyFunction_O::defineForMatter)
	.def("useDefaultSettings",&EnergyFunction_O::useDefaultSettings)
	.def("checkForBeyondThresholdInteractions",&EnergyFunction_O::checkForBeyondThresholdInteractions)
	.def("addTermsForListOfRestraints",&EnergyFunction_O::addTermsForListOfRestraints)
	.def("setName",&EnergyFunction_O::setName)
	.def("getName",&EnergyFunction_O::getName)
	.def("getStretchComponent",&EnergyFunction_O::getStretchComponent)
	.def("getAngleComponent",&EnergyFunction_O::getAngleComponent)
	.def("getDihedralComponent",&EnergyFunction_O::getDihedralComponent)
	.def("getNonbondComponent",&EnergyFunction_O::getNonbondComponent)
	.def("getImproperRestraintComponent",&EnergyFunction_O::getImproperRestraintComponent)
	.def("getChiralRestraintComponent",&EnergyFunction_O::getChiralRestraintComponent)
	.def("getAnchorRestraintComponent",&EnergyFunction_O::getAnchorRestraintComponent)
	.def("getFixedNonbondRestraintComponent",&EnergyFunction_O::getFixedNonbondRestraintComponent)
//	.def("getImproperRestraintComponentEnergy",&EnergyFunction_O::getImproperRestraintComponentEnergy)
//	.def("getDihedralComponentEnergy",&EnergyFunction_O::getDihedralComponentEnergy)
	.def("getMessage",&EnergyFunction_O::getMessage)
	.def("getMatter",&EnergyFunction_O::getMatter)
	.def("getNVectorSize",&EnergyFunction_O::getNVectorSize)
	.def("evaluateAll",&EnergyFunction_O::evaluateAll)
	.def("evaluateEnergy",&EnergyFunction_O::evaluateEnergy)
	.def("evaluateEnergyForce",&EnergyFunction_O::evaluateEnergyForce)
	.def("evaluateEnergyForceFullHessian",&EnergyFunction_O::evaluateEnergyForceFullHessian)
	.def("evaluateEnergyForceFullHessianForDebugging",&EnergyFunction_O::evaluateEnergyForceFullHessianForDebugging)
	.def("evaluateNumericalForce",&EnergyFunction_O::evaluateNumericalForce)
	.def("evaluateNumericalHessian",&EnergyFunction_O::evaluateNumericalHessian)
	.def("calculateEnergy",&EnergyFunction_O::calculateEnergy)
	.def("calculateEnergyAndForce",&EnergyFunction_O::calculateEnergyAndForce)
	.def("extractCoordinatesFromAtoms",&EnergyFunction_O::extractCoordinatesFromAtoms)
	.def("writeCoordinatesToAtoms",&EnergyFunction_O::writeCoordinatesToAtoms)
	.def("writeCoordinatesAndForceToAtoms",&EnergyFunction_O::writeCoordinatesAndForceToAtoms)
	.def("writeForceToAtoms",&EnergyFunction_O::writeForceToAtoms)
	.def("getTotalEnergy",&EnergyFunction_O::getTotalEnergy)
	.def("dumpTerms",&EnergyFunction_O::dumpTerms)
	.def("calculateNumericalDerivative",&EnergyFunction_O::calculateNumericalDerivative)
	.def("calculateNumericalSecondDerivative",
		&EnergyFunction_O::calculateNumericalSecondDerivative)
	.def("energyTermsEnabled",&EnergyFunction_O::energyTermsEnabled)
	.def("summarizeBeyondThresholdInteractionsAsString",
		&EnergyFunction_O::summarizeBeyondThresholdInteractionsAsString)
	.def("countTermsBeyondThreshold",&EnergyFunction_O::countTermsBeyondThreshold)
	.def("setupHessianPreconditioner",&EnergyFunction_O::setupHessianPreconditioner)
	.def("getMissingParameters",&EnergyFunction_O::getMissingParameters)
	.def("countBadVdwInteractions",&EnergyFunction_O::countBadVdwInteractions)
	.def("hasMissingParameters",&EnergyFunction_O::hasMissingParameters)
	.def("summarizeEnergyAsString",&EnergyFunction_O::summarizeEnergyAsString)
	.def("compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition",
		&EnergyFunction_O::compareAnalyticalAndNumericalForceAndHessianTermByTermAtCurrentPosition)
	;
#endif
//	defInPackage(ChemPkg,"create_EnergyFunction", &EnergyFunction_O::create);
}



    void EnergyFunction_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
	energyFunction_initializeSmarts(lisp);
    };


void	EnergyFunction_O::initialize()
{
    this->Base::initialize();
    this->_NonbondCrossTermTable = _Nil<FFNonbondCrossTermTable_O>();
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
    this->setName("");
    this->_Message = "";
    this->useDefaultSettings();
}

void	EnergyFunction_O::useDefaultSettings()
{_G();
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



uint	EnergyFunction_O::getNVectorSize() 
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
	SIMPLE_ERROR(BF("Unknown EnergyFunction setOption keyword[%s]") % option->__repr__() );
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
void	EnergyFunction_O::setOptions( core::Cons_sp options )
{_G();
    while ( options.notnilp() )
    {
	core::Symbol_sp option = core::oCar(options).as<core::Symbol_O>();
	options = options->cdr();
	core::T_sp val = core::oCar(options);
	options = options->cdr();
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

void	EnergyFunction_O::setupHessianPreconditioner(NVector_sp nvPosition,
					AbstractLargeSquareMatrix_sp m )
{_G();
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
		NVector_sp 	force,
       		bool		calcDiagonalHessian,
		bool		calcOffDiagonalHessian,
		AbstractLargeSquareMatrix_sp	hessian,
		NVector_sp	hdvec, NVector_sp dvec)
{_G()

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
bool	hasHessian = hessian.notnilp();
bool	hasHdAndD = (hdvec.notnilp())&&(dvec.notnilp());

#ifdef	DEBUG_ON //[
	// Summarize entry state for debugging
    LOG(BF("calcForce = %d") % calcForce  );
    LOG(BF("calcDiagonalHessian = %d") % calcDiagonalHessian  );
    LOG(BF("calcOffDiagonalHessian = %d") % calcOffDiagonalHessian  );
    LOG(BF("hasForce = %d") % hasForce  );
    LOG(BF("hasHessian = %d") % hasHessian  );
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
    try {
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

    } catch ( InteractionProblem ld ) {
////        _lisp->profiler().popTimerStates();
	SIMPLE_ERROR(BF("Interaction problem: %s")% ld.message() );
    }
////    _lisp->profiler().popTimerStates();

    	// More energy terms
    return this->_TotalEnergy;
}


string EnergyFunction_O::energyComponentsAsString()
{_G();
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
    NVector_sp pos = NVector_O::create(this->getNVectorSize(),_lisp);
    this->extractCoordinatesFromAtoms(pos);
    return this->compareAnalyticalAndNumericalForceAndHessianTermByTerm(pos);
}






/*!
 * Look for bad geometry
 *
 * Extract the coordinates of the atoms
 */
uint	EnergyFunction_O::checkForBeyondThresholdInteractions( )
{_G();
NVector_sp	pos;
stringstream	info;
int	fails = 0;

    info.str("");
    pos = NVector_O::create(this->getNVectorSize(),_lisp);
    this->extractCoordinatesFromAtoms(pos);
    fails = 0;
    fails += this->_Stretch->checkForBeyondThresholdInteractions(info,pos);
#if USE_ALL_ENERGY_COMPONENTS
    fails += this->_Angle->checkForBeyondThresholdInteractions(info,pos);
    fails += this->_Dihedral->checkForBeyondThresholdInteractions(info,pos);
    fails += this->_Nonbond->checkForBeyondThresholdInteractions(info,pos);
    fails += this->_ImproperRestraint->checkForBeyondThresholdInteractions(info,pos);
    fails += this->_ChiralRestraint->checkForBeyondThresholdInteractions(info,pos);
    fails += this->_AnchorRestraint->checkForBeyondThresholdInteractions(info,pos);
    fails += this->_FixedNonbondRestraint->checkForBeyondThresholdInteractions(info,pos);
#endif
    LOG(BF("Checking Done") );
    this->_Message = info.str();
    return fails;
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











double	EnergyFunction_O::evaluateEnergy( NVector_sp pos )
{
double		energy;
AbstractLargeSquareMatrix_sp	dummyMatrix;
NVector_sp	dummyForce, dummyHd, dummyD;
dummyMatrix = _Nil<AbstractLargeSquareMatrix_O>();
dummyForce = _Nil<NVector_O>();
dummyHd = _Nil<NVector_O>();
dummyD = _Nil<NVector_O>();
    energy = this->evaluateAll(pos,
	    			false, dummyForce,
				false, false, dummyMatrix,
				dummyHd, dummyD );
    return energy;
}



double	EnergyFunction_O::evaluateEnergyForce( NVector_sp pos, bool calcForce, NVector_sp force )
{_G();
double	energy;
NVector_sp	rawGrad;
NVector_sp	dummyHd;
NVector_sp	dummyD;
AbstractLargeSquareMatrix_sp		dummyMatrix;
dummyMatrix = _Nil<AbstractLargeSquareMatrix_O>();
dummyHd = _Nil<NVector_O>();
dummyD = _Nil<NVector_O>();
    if ( calcForce ) {
	rawGrad = force;
    } else {
	rawGrad = _Nil<NVector_O>();
    }
    energy = this->evaluateAll(pos,calcForce, rawGrad, false, false, dummyMatrix,
				dummyHd, dummyD );
    return energy;
}


double	EnergyFunction_O::evaluateEnergyForceFullHessian(
				NVector_sp pos,
				bool calcForce, NVector_sp force,
				bool calcDiagonalHessian,
				bool calcOffDiagonalHessian,
				AbstractLargeSquareMatrix_sp hessian )
{
double	energy;
NVector_sp	dummyHd, dummyD;
NVector_sp	rawPos, rawGrad;
dummyHd = _Nil<NVector_O>();
dummyD = _Nil<NVector_O>();
    if ( calcForce ) rawGrad = force;
    else rawGrad = _Nil<NVector_O>();
    energy = this->evaluateAll( pos,
	    			calcForce, rawGrad,
				calcDiagonalHessian,
				calcOffDiagonalHessian,
				hessian, dummyHd, dummyD );
    return energy;
}


double	EnergyFunction_O::evaluateEnergyForceFullHessianForDebugging()
{
NVector_sp	pos, force;
AbstractLargeSquareMatrix_sp	hessian;
double		energy;
NVector_sp	dummyHd, dummyD;
dummyHd = _Nil<NVector_O>();
dummyD = _Nil<NVector_O>();
    pos = NVector_O::create(this->getNVectorSize(),_lisp);
    force = NVector_O::create(this->getNVectorSize(),_lisp);
    hessian = FullLargeSquareMatrix_O::create(_lisp,this->getNVectorSize(),SymmetricDiagonalLower);
    this->extractCoordinatesFromAtoms(pos);
    energy = this->evaluateAll(pos, true, force,
				true, true, hessian,
				dummyHd, dummyD );
    return energy;
}





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
{_G();
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

    numForce = NVector_O::create(pos->size(),_lisp);
    this->evaluateNumericalForce(pos,numForce,DELTA);
    dot = numForce->dotProduct(analyticalForce);
    numericalMag = numForce->magnitude();
    analyticalMag = analyticalForce->magnitude();
    tempForce = NVector_O::create(pos->size(),_lisp);
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
	this->writeCoordinatesAndForceToAtoms(pos,analyticalForce);
	goto DONE;
    }
    if ( dot < 0.98 ) {
	result.str("");
	result << "The Lengths of Analytical and Numerical forces are very similar but they are not parallel%" << std::endl;
	result << "|analyticalForce| == " << analyticalMag << "  |numericalForce| == "<< numericalMag << std::endl;
	result << "(analyticalForce/|analyticalForce|).(numericalForce/|numericalForce|) = "<< dot << std::endl;
	report->_Message = result.str();
	this->writeCoordinatesAndForceToAtoms(pos,analyticalForce);
	goto DONE;
    }
    report->_Message = "Analytical and Numerical forces are virtually identical";
DONE:
    return report;
}







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




void EnergyFunction_O::_applyRestraints(ForceField_sp forceField, core::Iterator_sp restraintIterator )
{_G();
#if USE_ALL_ENERGY_COMPONENTS
    restraintIterator->first();
    while ( restraintIterator->notDone() )
    {
	Restraint_sp restraint = restraintIterator->current<Restraint_O>();
	if ( restraint.isA<RestraintDihedral_O>() )
	{
	    RestraintDihedral_sp dih = downcast<RestraintDihedral_O>(restraint);
	    EnergyImproperRestraint energyTerm;
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
	} else if ( restraint.isA<RestraintAnchor_O>() )
	{
	    RestraintAnchor_sp anchor = downcast<RestraintAnchor_O>(restraint);
	    EnergyAnchorRestraint	iterm;
	    Vector3		anchorPos;
	    Atom_sp a1 = anchor->getAtom();
	    EnergyAtom* ea1 = this->getEnergyAtomPointer(a1);
	    LOG(BF("Create an anchor restraint for %s") % a1->description()  );
	    anchorPos = anchor->getAnchorPos();
	    iterm.term.xa = anchorPos.getX();
	    iterm.term.ya = anchorPos.getY();
	    iterm.term.za = anchorPos.getZ();
	    iterm.term.ka = this->_AnchorRestraintWeight;
	    iterm.term.I1 = ea1->coordinateIndexTimes3();
	    this->_AnchorRestraint->addTerm(iterm);
	} else if ( restraint.isA<RestraintFixedNonbond_O>() )
	{
	    this->_FixedNonbondRestraint->setupForEvaluation(this->_AtomTable,this->_NonbondCrossTermTable);
	    RestraintFixedNonbond_sp fixedNonbond = restraint.as<RestraintFixedNonbond_O>();
	    Matter_sp matter = fixedNonbond->getMatter();
	    FFNonbondDb_sp nonbondDb = forceField->getNonbondDb();
//	    EnergyAtom	energyAtom(_lisp);
	    Loop loop;
	    {_BLOCK_TRACE("Defining FixedNonbond ATOMS");
		loop.loopTopGoal(matter,ATOMS);
		while ( loop.advanceLoopAndProcess() ) 
		{
		    Atom_sp a1 = loop.getAtom();
		    if ( a1->isAssignableTo<VirtualAtom_O>() ) continue; // skip virtuals
		    this->_FixedNonbondRestraint->addFixedAtom(nonbondDb,a1);
		}
	    }
	} else 
	{
	    IMPLEMENT_ME(); // Implement other types of restraints
	}
	restraintIterator->next();
    }
#endif
}

void EnergyFunction_O::_applyDihedralRestraint(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, double minDegrees, double maxDegrees, double weight)
{_OF();
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

void EnergyFunction_O::__createSecondaryAmideRestraints(gctools::Vec0<Atom_sp>& nitrogens)
{_G();
    gctools::Vec0<Atom_sp>::iterator ni;
    double transMin = -160.0;
    double transMax = 160.0;
    double cisMin = 20.0;
    double cisMax = -20.0;
    double weight = 1.0;
    for ( ni=nitrogens.begin(); ni!=nitrogens.end(); ni++ )
    {
	if ( secondaryAmide->matches(*ni) )
	{
	    ChemInfoMatch_sp match = secondaryAmide->getMatch();
	    Atom_sp ax = match->tag(chemkw::_sym_1);
	    Atom_sp ax1 = match->tag(chemkw::_sym_2);
	    Atom_sp ax2 = match->tag(chemkw::_sym_3);
	    Atom_sp ay = match->tag(chemkw::_sym_4);	// Carbonyl carbon
	    Atom_sp ay1 = match->tag(chemkw::_sym_5);
	    Atom_sp ay2 = match->tag(chemkw::_sym_6);
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
	    this->_applyDihedralRestraint(ax1,ax,ay,ay1,cisMin,cisMax,weight);
	    LOG(BF("Restrain cis %s - %s - %s -%s") % ax1->description() % ax->description() % ay->description() % ay1->description()  );
	    this->_applyDihedralRestraint(ax1,ax,ay,ay2,transMin,transMax,weight);
	    LOG(BF("Restrain trans %s - %s - %s -%s") % ax1->description() % ax->description() % ay->description() % ay2->description()  );
	    this->_applyDihedralRestraint(ax2,ax,ay,ay1,transMin,transMax,weight);
	    LOG(BF("Restrain trans %s - %s - %s -%s") % ax2->description() % ax->description() % ay->description() % ay1->description()  );
	    this->_applyDihedralRestraint(ax2,ax,ay,ay2,cisMin,cisMax,weight);
	    LOG(BF("Restrain cis %s - %s - %s -%s") % ax2->description() % ax->description() % ay->description() % ay2->description()  );
	}
    }
}





void EnergyFunction_O::defineForMatter(Matter_sp matter, ForceField_sp forceField)
{_G();
Loop            loop;
Atom_sp          a1, a2, a3, a4, aImproperCenter;
core::Symbol_sp  t1, t2, t3, t4, t141, t144;
EnergyAtom      *eaCenter, *ea1, *ea2, *ea3, *ea4;
EnergyStretch   energyStretch;
FFAngle_sp       ffAngle;
FFPtor_sp        ffPtor;
FFItor_sp        ffItor;
FFNonbond_sp	ffNonbond1, ffNonbond2;
int             coordinateIndex;



	//
	// Identify rings
	//
    {_BLOCK_TRACEF(BF("Identify rings"));
//	printf("%s:%d Skipping ring searching FIX ME\n", __FILE__, __LINE__);
	RingFinder_O::identifyRings(matter);
    }
	//
	// Assign relative Cahn-Ingold-Preylog priorities
	//

    {_BLOCK_TRACEF(BF("Assign CIP priorities"));
	CipPrioritizer_O::assignPriorities(matter); 
    }
	// 
	// Assign atom types
	//

    forceField->assignTypes(matter);
//    printf( "%s:%d dumping xml after assignTypes \n%s\n", __FILE__, __LINE__,matter->asXmlString() );

	//
	// Now create the energy function from all this info
	//

    if ( !(matter.isA<Aggregate_O>() || matter.isA<Molecule_O>() ) )
    {
        SIMPLE_ERROR(BF("You can only define energy functions for Aggregates or Molecules"));
    }
    	//
	// Define a Nonbond cross term table
	//
    GC_ALLOCATE(FFNonbondCrossTermTable_O, temp );
    this->_NonbondCrossTermTable = temp;
    this->_NonbondCrossTermTable->fillUsingFFNonbondDb(forceField->getNonbondDb());

    	//
	// Initialize the energy components
	//
    ALL_ENERGY_COMPONENTS(initialize());
    this->_Matter= matter;
    this->_DielectricConstant = 80.0;
    this->_eraseMissingParameters();
    coordinateIndex = 0;
    ASSERTNOTNULL(forceField);
#if 1
    this->_AtomTable->constructFromMatter(matter,forceField);
#else
    {_BLOCK_TRACE("Defining ATOMS");
	loop.loopTopGoal(matter,ATOMS);
	while ( loop.advanceLoopAndProcess() ) 
	{
	    a1 = loop.getAtom();
	    if ( a1->isAssignableTo<VirtualAtom_O>() ) continue; // skip virtuals
	    LOG(BF( "Atom = %s")% a1->description() );
	    a1->setAtomTableIndex(this->_AtomTable->getNumberOfAtoms());
	    EnergyAtom energyAtom(forceField,a1,coordinateIndex);
	    this->_AtomTable->add(energyAtom);
	    coordinateIndex += 3;
	}
    }
#endif
#ifdef	DEBUG_DEFINE_ENERGY
    lisp->print(BF("%s:%d There were %d atoms") % __FILE__ % __LINE__ % this->_AtomTable.size() );
#endif
    {_BLOCK_TRACE("Defining STRETCH");
	ASSERTNOTNULL(forceField->_Stretches);
	loop.loopTopGoal(matter,BONDS);
	while ( loop.advanceLoopAndProcess() ) {
	    a1 = loop.getBondA1();
	    a2 = loop.getBondA2();
	    t1 = a1->getType();
	    t2 = a2->getType();
	    ea1 = this->getEnergyAtomPointer(a1);
	    ea2 = this->getEnergyAtomPointer(a2);
	    FFStretch_sp ffStretch = forceField->_Stretches->findTerm(a1,a2);
	    ASSERTNOTNULL(ffStretch);
	    if ( ffStretch->level() != parameterized ) {
		this->_addMissingParameter(ffStretch->levelDescription());
	    }
	    energyStretch.defineFrom(ffStretch,ea1,ea2,this->_Stretch->getScale());
	    this->_Stretch->addTerm(energyStretch);
	}

    }
#ifdef	DEBUG_DEFINE_ENERGY
    lisp->print(BF("%s:%d There were %d stretch terms") % __FILE__ % __LINE__ % this->_Stretch.size() );
#endif
#if USE_ALL_ENERGY_COMPONENTS
    {_BLOCK_TRACE("Defining ANGLES");
	loop.loopTopGoal(matter,ANGLES);
	while ( loop.advanceLoopAndProcess() )
	{
	    a1 = loop.getAtom1();
	    a2 = loop.getAtom2();
	    a3 = loop.getAtom3();
	    ea1 = this->getEnergyAtomPointer(a1);
	    ea2 = this->getEnergyAtomPointer(a2);
	    ea3 = this->getEnergyAtomPointer(a3);
	    ffAngle = forceField->_Angles->findTerm(a1,a2,a3);
	    if ( ffAngle->level() != parameterized ) {
		LOG(BF("Missing angle parameter between types: %s-%s-%s") % a1->getType()% a2->getTypeString()% a3->getTypeString() );
		this->_addMissingParameter(ffAngle->levelDescription());
		LOG(BF("Added to missing parameters") );
	    }
	    EnergyAngle energyAngle;
	    energyAngle.defineFrom(ffAngle,ea1,ea2,ea3,this->_Angle->getScale());
	    this->_Angle->addTerm(energyAngle);
	}
    }
    {_BLOCK_TRACE("Defining PROPERS");
	EnergyNonbond energyNonbond;
	loop.loopTopGoal(matter,PROPERS);
	while ( loop.advanceLoopAndProcess() ) 
	{
	    a1 = loop.getAtom1();
	    a2 = loop.getAtom2();
	    a3 = loop.getAtom3();
	    a4 = loop.getAtom4();
	    t1 = a1->getType();
	    t2 = a2->getType();
	    t3 = a3->getType();
	    t4 = a4->getType();
	    ea1 = this->getEnergyAtomPointer(a1);
	    ea2 = this->getEnergyAtomPointer(a2);
	    ea3 = this->getEnergyAtomPointer(a3);
	    ea4 = this->getEnergyAtomPointer(a4);
	    if ( !forceField->_Ptors->hasBestTerm(t1,t2,t3,t4) ) 
	    {
		stringstream ss;
		ss << "Could not find ptor term(";
		ss << _rep_(t1) << "-";
		ss << _rep_(t2) << "-";
		ss << _rep_(t3) << "-";
		ss << _rep_(t4) << ")";
		ss << " for atoms(";
		ss << _rep_(a1->getName()) << "-";
		ss << _rep_(a2->getName()) << "-";
		ss << _rep_(a3->getName()) << "-";
		ss << _rep_(a4->getName()) << ")";
		_lisp->print(BF("%s") % ss.str() );
#if 0
		SIMPLE_ERROR(BF("%s")%ss.str());
//		forceField->_Ptors->cantFind(t1,t2,t3,t4);
//		energyDihedral.defineMissingProper(ea1,ea2,ea3,ea4);
//		this->_MissingDihedralTerms.push_back(energyDihedral);
#endif
	    } else 
	    {
		ffPtor = forceField->_Ptors->findBestTerm(t1,t2,t3,t4);
		for ( int n=1;n<=FFPtor_O::MaxPeriodicity; n++ ) 
		{
		    if ( ffPtor->hasPeriodicity(n) )
		    {
			LOG(BF( "Adding proper term for atoms %s-%s-%s-%s types: %s-%s-%s-%s")%
				    ea1->getResidueAndName()
				    % ea2->getResidueAndName()
				    % ea3->getResidueAndName()
				    % ea4->getResidueAndName()
				    % t1 % t2 % t3 % t4
			     );
			EnergyDihedral energyDihedral;
			energyDihedral.defineFrom(n,ffPtor,ea1,ea2,ea3,ea4,this->_Dihedral->getScale());
			this->_Dihedral->addTerm(energyDihedral);
		    }
		}
		if (ea1->inBondOrAngle(ea4->atom()) )
		{
#ifdef	DEBUG_ON
		    if ( t1 < t4 ) {
			t141 = t1;
			t144 = t4;
		    }else{
			t141 = t4;
			t144 = t1;
		    }
		    LOG(BF("Defining 1-4 interaction %-9s- %-9s   ") % t1 % t4 );
#endif
#if 0
		    if ( energyNonbond.defineFrom(forceField, true, ea1, ea4, this->_Nonbond ) ) 
		    {
		        this->_Nonbond->addTerm(energyNonbond);
		    }
#endif
		} else {
#ifdef	DEBUG_ON
		    if ( t1 < t4 ) {
			t141 = t1;
			t144 = t4;
		    }else{
			t141 = t4;
			t144 = t1;
		    }
		    LOG(BF("Ignoring 1-4 interaction %-9s- %-9s    ") % t1 % t4 );
#endif
		}
//		ea1->_CloserThan15.insert(ea4->_Atom);
//		ea4->_CloserThan15.insert(ea1->_Atom);
	    }
	}
    }
    {_BLOCK_TRACE("Defining IMPROPERS");
	aImproperCenter = _Nil<Atom_O>();
	EnergyDihedral energyDihedral;
	loop.loopTopGoal(matter,IMPROPERS);
	while ( loop.advanceLoopAndProcess() ) 
	{
	    a1 = loop.getAtom1();
	    a2 = loop.getAtom2();
	    a3 = loop.getAtom3();
	    a4 = loop.getAtom4();
	    t1 = a1->getType();
	    t2 = a2->getType();
	    t3 = a3->getType();
	    t4 = a4->getType();
	    ea1 = this->getEnergyAtomPointer(a1);
	    ea2 = this->getEnergyAtomPointer(a2);
	    ea3 = this->getEnergyAtomPointer(a3);
	    ea4 = this->getEnergyAtomPointer(a4);
	    if ( forceField->_Itors->hasBestTerm(t1,t2,t3,t4) ) 
	    {
		    //
		    // Only one improper per central atom.
		    // We may not get exactly the same improper as AMBER does
		    //
		if ( a3 != aImproperCenter ) {
		    ffItor = forceField->_Itors->findBestTerm(t1,t2,t3,t4);
		    for ( int n=1;n<=6; n++ ) {
			if ( ffItor->hasPeriodicity(n) )
			{
			    energyDihedral.defineFrom(n,ffItor,ea1,ea2,ea3,ea4,this->_Dihedral->getScale());
			    this->_Dihedral->addTerm(energyDihedral);
			    aImproperCenter = a3;
			}
		    }
		}
	    }
	}
    }


        // Nonbonds here!!!!!!!!!!!!!!
#if 1
    this->_Nonbond->constructFromAtomTable(this->_AtomTable, forceField);
#else
    {_BLOCK_TRACE("Defining NONBONDED");
	if ( this->_AtomTable->getNumberOfAtoms() > 2 )
	{
	    LOG(BF("Defining NONBONDS") );
            gctools::Vec0<EnergyAtom>::iterator	iea1;
            gctools::Vec0<EnergyAtom>::iterator	iea2;
	    for ( iea1 = this->_AtomTable->begin();
		  iea1 != this->_AtomTable->end()-1; iea1++ ) { //Was iea1 != this->_AtomTable->end()-1
		for ( iea2 = iea1+1; iea2 != this->_AtomTable->end(); iea2++ ) { // Was iea2 != this->_AtomTable->end()
		    if (!iea1->inBondOrAngle(iea2->_Atom)) 
		    {
			bool in14 = iea1->relatedBy14(iea2->_Atom);
			EnergyNonbond energyNonbond(_lisp);
			if ( energyNonbond.defineFrom(forceField, in14,
						      &(*iea1),&(*iea2),this->_Nonbond) )  {
			    this->_Nonbond->addTerm(energyNonbond);
			}
		    }
		}
	    }
	} else
	{
	    LOG(BF("There are no non-bonds"));
	}
    }
#endif


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
		core::Cons_sp priority = a1->getNeighborsByRelativePriority();
		ASSERTP(priority->length() == 4, "There must be 4 neighbors to assign stereochemistry");
		core::Cons_sp cur = priority;
		n1 = cur->car<Atom_O>();
		cur = cur->cdr();
		n2 = cur->car<Atom_O>();
		cur = cur->cdr();
		n3 = cur->car<Atom_O>();
		cur = cur->cdr();
		n4 = cur->car<Atom_O>();
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
	    if ( a->getElement() == element_N )
	    {
		nitrogens.push_back(a);
	    }
	}
        this->__createSecondaryAmideRestraints(nitrogens);
    } else
    {
	LOG(BF("Skipping Secondary amide restraints because _RestrainSecondaryAmides = %d") % this->_RestrainSecondaryAmides );
    }

#endif // USE_ALL_ENERGY_COMPONENTS

	//
	// Set up force-field restraints
	//
    {_BLOCK_TRACE("Defining force-field restraints");
	IterateRestraints_sp restraintIt = IterateRestraints_O::create(_lisp,matter);
	this->_applyRestraints(forceField,restraintIt);
    }
    LOG(BF("Done terms") );
}



void	EnergyFunction_O::addTermsForListOfRestraints(ForceField_sp forceField, core::Cons_sp restraintList)
{_G();
    addon::IterateCons_sp	iterate;
    iterate = addon::IterateCons_O::create(restraintList);
    this->_applyRestraints(forceField,iterate);
}

void	EnergyFunction_O::extractCoordinatesFromAtoms(NVector_sp pos)
{_G();
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



void    EnergyFunction_O::writeCoordinatesToAtoms(NVector_sp pos)
{_G();
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

void    EnergyFunction_O::writeForceToAtoms(NVector_sp force)
{
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
        ai->atom()->setForce(v);
    }
}


void    EnergyFunction_O::writeCoordinatesAndForceToAtoms(NVector_sp pos, NVector_sp force)
{
    this->writeCoordinatesToAtoms(pos);
    this->writeForceToAtoms(force);
}


EnergyAtom*	EnergyFunction_O::getEnergyAtomPointer(Atom_sp a)
{ 
    return this->_AtomTable->getEnergyAtomPointer(a);
};



double	EnergyFunction_O::calculateEnergy( )
{_G();
NVector_sp	pos;
    pos = NVector_O::create(this->getNVectorSize(),_lisp);
    this->extractCoordinatesFromAtoms(pos);
    return this->evaluateEnergy(pos);
}


double EnergyFunction_O::calculateEnergyAndForce( )
{
NVector_sp	pos;
NVector_sp	force;
double		energy;
    pos = NVector_O::create(this->getNVectorSize(),_lisp);
    force = NVector_O::create(this->getNVectorSize(),_lisp);
    this->extractCoordinatesFromAtoms(pos);
    energy = this->evaluateEnergyForce(pos,true,force);
    	// To calculate the force magnitude use force->magnitude();
    	// To calculate the force rmsMagnitude use force->rmsMagnitude();
    this->writeForceToAtoms(force);
    return energy;
}







void	EnergyFunction_O::dealWithProblem(const InteractionProblem& problem)
{_OF();
    switch ( problem._Type ) {
	case linearDihedral:
	    problem._Atom1->bumpPosition(0.1);
	    problem._Atom2->bumpPosition(0.1);
	    problem._Atom3->bumpPosition(0.1);
	    problem._Atom4->bumpPosition(0.1);
	    break;
	case linearImproperRestraint:
	    problem._Atom1->bumpPosition(0.1);
	    problem._Atom2->bumpPosition(0.1);
	    problem._Atom3->bumpPosition(0.1);
	    problem._Atom4->bumpPosition(0.1);
	    break;
	case linearAngle:
	    problem._Atom1->bumpPosition(0.1);
	    problem._Atom2->bumpPosition(0.1);
	    problem._Atom3->bumpPosition(0.1);
	    break;
	default:
	    SIMPLE_ERROR(BF("EnergyFunction_O::dealWithProblem>> I am not handling this problem yet"));
    }
}



string	EnergyFunction_O::summarizeBeyondThresholdInteractionsAsString()
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


string	EnergyFunction_O::summarizeEnergyAsString()
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



bool	EnergyFunction_O::hasMissingParameters()
{
    if ( this->_MissingParameters.size() > 0 ) {
	return true;
    }
    return false;
}


string	EnergyFunction_O::getMissingParameters()
{
set<string>::iterator	si;
stringstream		ss;
    for ( si=this->_MissingParameters.begin();
		si!=this->_MissingParameters.end(); si++ ) {
	ss << (*si) << std::endl;
    }
    return ss.str();
}


string	EnergyFunction_O::debugLogAsString()
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


void EnergyFunction_O::enableDebug()
{
    ALL_ENERGY_COMPONENTS(enableDebug());
}

void EnergyFunction_O::disableDebug()
{
    ALL_ENERGY_COMPONENTS(disableDebug());
}


#ifdef RENDER
uint	EnergyFunction_O::countBadVdwInteractions(double scaleSumOfVdwRadii, addon::DisplayList_sp displayIn)
{_G();
NVector_sp	pos;
EnergyNonbond_sp	nbComponent;
    uint i = 0;
#if USE_ALL_ENERGY_COMPONENTS
    pos = NVector_O::create(this->getNVectorSize(),_lisp);
    this->extractCoordinatesFromAtoms(pos);
    nbComponent = this->getNonbondComponent();
    i = nbComponent->countBadVdwOverlaps(scaleSumOfVdwRadii,pos,displayIn,displayIn->lisp());
#endif
    return i;
}
#endif




EXPOSE_CLASS_AND_GLOBALS(chem,EnergyFunction_O);
};





