// #undef USEBOOSTPYTHON

/* Currently using USE_POSIX_TIME=1 will fail to compile when USEBOOSTPYTHON is defined
   this is a bug in the boost libraries in including posix time header files.
   I have to fix this - for now, disable USE_POSIX_TIME */
#define USE_POSIX_TIME 0


#define	DEBUG_LEVEL_NONE

//#define	TURN_LINESEARCH_DETAILS_ON
#define	TURN_LINESEARCH_DETAILS_OFF


#define	VALIDATE_FORCE_OFF	// Turn this on to compare analytical force to numerical force at every step

#define	USE_NEW_CONJUGATE_GRADIENTS
#include <clasp/core/common.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/symbolTable.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/minimizer.h>
#include <iostream>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#if USE_POSIX_TIME
#include <clasp/core/posixTime.h>
#endif // USE_POSIX_TIME
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
//#include "core/environment.h"
#include <clasp/core/str.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/linearAlgebra.h>
//#include "core/xmlSaveArchive.h"
#include <cando/chem/matter.h>
//#include "core/lispCallback.h"
#include <cando/chem/forceField.h>
#include <cando/chem/energyStretch.h>
#if USE_ALL_ENERGY_COMPONENTS
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyImproperRestraint.h>
#include <cando/chem/energyFixedNonbond.h>
#endif
#include <clasp/core/wrappers.h>


SYMBOL_EXPORT_SC_(ChemKwPkg,message);

SYMBOL_SC_(ChemPkg,_PLUS_minimizerStatusConverter_PLUS_);
SYMBOL_SC_(ChemPkg,_PLUS_preconditionerTypeConverter_PLUS_);


SYMBOL_EXPORT_SC_(ChemPkg,MinimizerExceededMaxStepsError);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerStuckError);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerError);

namespace chem
{


#define MINIMIZER_ERROR(msg) ERROR(_sym_MinimizerError,core::lisp_createList(kw::_sym_message,core::Str_O::create(msg)))
#define MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(msg) ERROR(_sym_MinimizerExceededMaxStepsError,core::lisp_createList(kw::_sym_message,core::Str_O::create(msg)))
#define MINIMIZER_STUCK_ERROR(msg) ERROR(_sym_MinimizerStuckError,core::lisp_createList(kw::_sym_message,core::Str_O::create(msg)))



//#define		DEBUGINNER	1
//#define		DEBUGTORSION	1
//#define		DEBUGANGLE	1
//#define		DEBUGBOND	1
//#define		DEBUGTORSION	1




#define	ITMAX 100
#define	EPS	1.0e-10
#define	TOL	2.0e-4
#define DBRENTITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);
#define	LEFT_SHIFT_VALUES(a,b,c,d) {(a)=(b);(b)=(c);(c)=(d);}
#define	SWAP_VALUES(a,b,t) {(t)=(a);(a)=(b);(b)=(t);}
#undef SIGN
#define	SIGN(a,b)	((b)>0.0?fabs(a):-fabs(a))
//#define	MAX(a,b)	((a)>(b)?(a):(b))
#define	TINY	1.0e-20
#define	GOLD	1.618034
#define	GLIMIT	100.0

#ifdef	PARMDEBUG
#define	IFPARMDEBUG(x)	x
#else
#define	IFPARMDEBUG(x)	{}
#endif


#define NOTFOUND                        -1
#define MAXSTEEPESTDESCENTSTEPS         200           /* Should be 200 */
#define MAXCONJUGATEGRADIENTSTEPS       10000           /* Should be 200 */
#define MAXTRUNCATEDNEWTONSTEPS		500           /* Should be 200 */
#define MAXLINESEARCHSTEPS              10
#define STARTSTEPSIZE                   0.0001
#define MIN_GRADIENT_MEAN		0.000001





    string	stringForPreconditionerType(PreconditionerType t)
    {
	switch (t)
	{
	case noPreconditioner:
	    return "noPreconditioner";
	case hessianPreconditioner:
	    return "hessianPreconditioner";
	default:
	    return "-unknown preconditioner-";
	}
	return "-unknown preconditioner-";
    }

    string	shortStringForPreconditionerType(PreconditionerType t) {
	switch (t) {
	case noPreconditioner:
	    return "nP";
	case hessianPreconditioner:
	    return "hP";
	default:
	    return "??";
	}
	return "?P";
    }


    PreconditionerType	preconditionerTypeFromString(const string& t )
    {
	if ( t=="noPreconditioner" ) return noPreconditioner;
	if ( t=="hessianPreconditioner" ) return hessianPreconditioner;
	return unknownPreconditioner;
    }




    const char* minimizerOptions[] = {":showElapsedTime",""};
    void	Minimizer_O::changeOptions(core::List_sp optionCons)
    {_G();
    IMPLEMENT_ME();
#if 0
    core::HashTableEq_sp options = core::HashTableEq_O::createFromKeywordCons(optionCons,minimizerOptions,_lisp);
    options->setBoolValueIfAvailable(this->_ShowElapsedTime,":showElapsedTime");
#endif
}


#if INIT_TO_FACTORIES

#define ARGS_Minimizer_O_make "(matter force_field given_energy_function)"
#define DECL_Minimizer_O_make ""
#define DOCS_Minimizer_O_make "make Minimizer"
  Minimizer_sp Minimizer_O::make(Matter_sp matter, ForceField_sp forceField, EnergyFunction_sp givenEnergyFunction)
  {_G();
      GC_ALLOCATE(Minimizer_O, me );
    bool initialized = false;
    if ( givenEnergyFunction.nilp() )
      {
	if ( matter.notnilp() && forceField.notnilp() )
	  {
	    initialized = true;
	    EnergyFunction_sp ef = EnergyFunction_O::create();
	    ef->defineForMatter(matter,forceField);
	    me->setEnergyFunction(ef);
	  }
      } else
      {
	if ( matter.nilp() && forceField.nilp() )
	  {
	    initialized = true;
	    me->setEnergyFunction(givenEnergyFunction);
	  }
      }
    if ( !initialized )
      {
	  SIMPLE_ERROR(BF("When initializing a Minimizer object you must either provide an (energyFunction) or a (Matter/ForceField) pair"));
      }
    return me;
  }

#else

    core::T_sp  Minimizer_O::__init__(core::Function_sp exec, core::List_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	this->Base::__init__(exec,args,env,lisp);
	bool initialized = false;
	Matter_sp matter = env->lookup(lisp->internWithPackageName(ChemPkg,"matter")).as<Matter_O>();
	ForceField_sp forceField = env->lookup(lisp->internWithPackageName(ChemPkg,"forceField")).as<ForceField_O>();
	EnergyFunction_sp givenEnergyFunction = env->lookup(lisp->internWithPackageName(ChemPkg,"energyFunction")).as<EnergyFunction_O>();
	if ( givenEnergyFunction.nilp() )
	{
	    if ( matter.notnilp() && forceField.notnilp() )
	    {
		initialized = true;
		EnergyFunction_sp ef = EnergyFunction_O::create();
		ef->defineForMatter(matter,forceField);
		this->setEnergyFunction(ef);
	    }
	} else
	{
	    if ( matter.nilp() && forceField.nilp() )
	    {
		initialized = true;
		this->setEnergyFunction(givenEnergyFunction);
	    }
	}
	if ( !initialized )
	{
	    SIMPLE_ERROR(BF("When initializing a Minimizer object you must either provide an (energyFunction) or a (Matter/ForceField) pair"));
	}
	return _Nil<core::T_O>();
    }

#endif

    void Minimizer_O::exposeCando(core::Lisp_sp lisp)
    {
	Defun_maker(ChemPkg,Minimizer);
	core::class_<Minimizer_O>()
//	    .def_raw("core:__init__",&Minimizer_O::__init__,"(self &key matter forceField energyFunction)")
	    .def("throwMinimizerExceededMaxSteps",&Minimizer_O::throwMinimizerExceededMaxSteps)
	    .def("throwMinimizerStuck",&Minimizer_O::throwMinimizerStuck)
	    .def("throwMinimizerError",&Minimizer_O::throwMinimizerError)

//	    .def("setStepCallback",&Minimizer_O::setStepCallback)
	    .def("setEnergyFunction",&Minimizer_O::setEnergyFunction)
	    .def("enablePrintIntermediateResults",&Minimizer_O::enablePrintIntermediateResults)
	    .def("disablePrintIntermediateResults",&Minimizer_O::disablePrintIntermediateResults)
	    .def("minimizeSteepestDescent",&Minimizer_O::minimizeSteepestDescent)
	    .def("minimizeConjugateGradient",&Minimizer_O::minimizeConjugateGradient)
	    .def("resetAndMinimize",&Minimizer_O::resetAndMinimize)
	    .def("minimize",&Minimizer_O::minimize)
	    .def("evaluateEnergyAndForceManyTimes",&Minimizer_O::evaluateEnergyAndForceManyTimes)
//	.def("asXml",&Minimizer_O::asXml)
	    .def("getEnergyFunction",&Minimizer_O::getEnergyFunction)
	    .def("getIteration",&Minimizer_O::getIteration)
	    .def("status",&Minimizer_O::statusAsInt)
#ifdef XML_ARCHIVE
	    .def("debugStart",&Minimizer_O::debugStart)
	    .def("debugStop",&Minimizer_O::debugStop)
#endif
	    .def("statusAsString",&Minimizer_O::statusAsString)
	    .def("statusMessage",&Minimizer_O::statusMessage)
	    .def("configurationAsString",&Minimizer_O::configurationAsString)
	    .def("useDefaultSettings",&Minimizer_O::useDefaultSettings)
	    .def("getRMSForce",&Minimizer_O::getRMSForce)
	    .def("restart",&Minimizer_O::restart)
	    .def("setMaximumNumberOfSteepestDescentSteps", &Minimizer_O::setMaximumNumberOfSteepestDescentSteps)
	    .def("setMaximumNumberOfConjugateGradientSteps", &Minimizer_O::setMaximumNumberOfConjugateGradientSteps)
	    .def("setMaximumNumberOfTruncatedNewtonSteps", &Minimizer_O::setMaximumNumberOfTruncatedNewtonSteps)
	    .def("setSteepestDescentTolerance", &Minimizer_O::setSteepestDescentTolerance)
	    .def("setConjugateGradientTolerance", &Minimizer_O::setConjugateGradientTolerance)
	    .def("setTruncatedNewtonTolerance", &Minimizer_O::setTruncatedNewtonTolerance)
	    ;
    }

    void Minimizer_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Minimizer,"","",_lisp)
	    .def("throwMinimizerExceededMaxSteps",&Minimizer_O::throwMinimizerExceededMaxSteps)
	    .def("throwMinimizerStuck",&Minimizer_O::throwMinimizerStuck)
	    .def("throwMinimizerError",&Minimizer_O::throwMinimizerError)

//	    .def("setStepCallback",&Minimizer_O::setStepCallback)
	    .def("setEnergyFunction",&Minimizer_O::setEnergyFunction)
	    .def("enablePrintIntermediateResults",&Minimizer_O::enablePrintIntermediateResults)
	    .def("disablePrintIntermediateResults",&Minimizer_O::disablePrintIntermediateResults)
	    .def("minimizeSteepestDescent",&Minimizer_O::minimizeSteepestDescent)
	    .def("minimizeConjugateGradient",&Minimizer_O::minimizeConjugateGradient)
	    .def("resetAndMinimize",&Minimizer_O::resetAndMinimize)
	    .def("minimize",&Minimizer_O::minimize)
	    .def("evaluateEnergyAndForceManyTimes",&Minimizer_O::evaluateEnergyAndForceManyTimes)
//	.def("asXml",&Minimizer_O::asXml)
	    .def("getEnergyFunction",&Minimizer_O::getEnergyFunction)
	    .def("getIteration",&Minimizer_O::getIteration)
	    .def("status",&Minimizer_O::statusAsInt)
	    .def("debugStart",&Minimizer_O::debugStart)
	    .def("debugStop",&Minimizer_O::debugStop)
	    .def("statusAsString",&Minimizer_O::statusAsString)
	    .def("statusMessage",&Minimizer_O::statusMessage)
	    .def("configurationAsString",&Minimizer_O::configurationAsString)
	    .def("useDefaultSettings",&Minimizer_O::useDefaultSettings)
	    .def("getRMSForce",&Minimizer_O::getRMSForce)
	    .def("restart",&Minimizer_O::restart)
	    .def("setMaximumNumberOfSteepestDescentSteps", &Minimizer_O::setMaximumNumberOfSteepestDescentSteps)
	    .def("setMaximumNumberOfConjugateGradientSteps", &Minimizer_O::setMaximumNumberOfConjugateGradientSteps)
	    .def("setMaximumNumberOfTruncatedNewtonSteps", &Minimizer_O::setMaximumNumberOfTruncatedNewtonSteps)
	    .def("setSteepestDescentTolerance", &Minimizer_O::setSteepestDescentTolerance)
	    .def("setConjugateGradientTolerance", &Minimizer_O::setConjugateGradientTolerance)
	    .def("setTruncatedNewtonTolerance", &Minimizer_O::setTruncatedNewtonTolerance)
	    ;
#endif
    }



    string	Minimizer_O::statusAsString()
    {
	string	status, precon;
	switch (this->_Status) {
	case minimizerIdle:
	    return "minimizerIdle";
	case minimizerSucceeded:
	    return "minimizerSucceeded";
	case minimizerError:
	    return "minimizerError";
	case steepestDescentRunning:
	    status = "steepestDescentRunning";
	    break;
	case conjugateGradientRunning:
	    status = "conjugateGradientRunning";
	    break;
	case truncatedNewtonRunning:
	    status = "truncatedNewtonRunning";
	    break;
	default:
    	    return "unknownMinimizerStatus";
	}
	// Running, append preconditioner
	status = status + "-" + stringForPreconditionerType(this->_CurrentPreconditioner);
	return status;
    }


#ifdef XML_ARCHIVE
    void	Minimizer_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("_InitialLineSearchStep",this->_InitialLineSearchStep);
	node->attribute("_ShowElapsedTime",this->_ShowElapsedTime);
	node->attribute("_SteepestDescentTolerance",this->_SteepestDescentTolerance);
	node->attribute("_NumberOfSteepestDescentSteps",this->_NumberOfSteepestDescentSteps);
	node->attribute("_NumberOfConjugateGradientSteps",this->_NumberOfConjugateGradientSteps);
	node->attribute("_ConjugateGradientTolerance",this->_ConjugateGradientTolerance);
	node->attribute("_NumberOfTruncatedNewtonSteps",this->_NumberOfTruncatedNewtonSteps);
	node->attribute("_TruncatedNewtonTolerance",this->_TruncatedNewtonTolerance);
	node->attribute("energyFunction",this->_EnergyFunction );

//    node->attribute("_ZZ",this->_ZZ);
    }
#endif

    string	Minimizer_O::statusAsShortString()
    {
	string	status, precon;
	switch (this->_Status) {
	case minimizerIdle:
	    return "IDLE";
	case minimizerSucceeded:
	    return "DONE";
	case minimizerError:
	    return "ERR!";
	case steepestDescentRunning:
	    status = "SD";
	    break;
	case conjugateGradientRunning:
	    status = "CG";
	    break;
	case truncatedNewtonRunning:
	    status = "TN";
	    break;
	default:
	    return "?st?";
	}
	status = status + shortStringForPreconditionerType(this->_CurrentPreconditioner);
	return status;
    }



#if 0
    void	Minimizer_O::setStepCallback(core::LispCallback_sp cb)
    {_G();
    this->_StepCallback = cb;
    }
#endif


#ifdef XML_ARCHIVE
    void	Minimizer_O::debugStart()
    {_G();
    this->_DebugOn = true;
    LOG(BF("Status") );
    this->_Log = MinimizerLog_O::create();
    LOG(BF("Status") );
    this->_Log->_Minimizer = this->sharedThis<Minimizer_O>();
    LOG(BF("Status") );
    }


    void	Minimizer_O::debugStop(const string& fileName)
    {_G();
    LOG(BF("writing minimizer log to (%s)") % fileName.c_str()  );
    if ( !this->_DebugOn ) return;
    ASSERTNOTNULL(this->_Log);
    core::XmlSaveArchive_sp arc;
    LOG(BF("Creating archive") );
    arc = core::XmlSaveArchive_O::create();
    LOG(BF("putting minimizerLog") );
    arc->put("minimizerLog",this->_Log);
    LOG(BF("Calling saveAs") );
    arc->saveAs(fileName);
    LOG(BF("Done saveAs") );
    }
#endif


/*
 *	getPosition
 *
 *	Calculate the position from the nvOrigin+nvDirection*x
 */
    void	Minimizer_O::getPosition( NVector_sp 	nvResult,
					  NVector_sp	nvOrigin,
					  NVector_sp	nvDirection,
					  double		x )
    {
	nvResult->XPlusYTimesScalar(nvOrigin,nvDirection,x);
    }


/*
 *      dTotalEnergy
 *
 *      Evaluate the total amber.
 *      Evaluate the amber at nvPos and return the derivative
 *      in nvForce (which must be an initialized vector,
 *      but does not need to be filled with zeros).
 */
    double	Minimizer_O::dTotalEnergy( NVector_sp nvPos )
    {
	double          dEnergy;

	dEnergy = this->_EnergyFunction->evaluateEnergy( nvPos );
	return(dEnergy);
    }


/*
 *      dTotalEnergyForce
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Evaluate the total amber.
 *      Evaluate the amber at nvPos and return the derivative
 *      in nvForce (which must be an initialized vector,
 *      but does not need to be filled with zeros).
 */
    double	Minimizer_O::dTotalEnergyForce( NVector_sp nvPos, NVector_sp nvForce)
    {_G();
	return this->_EnergyFunction->evaluateEnergyForce(nvPos,true,nvForce);
    }


/*
 *	d1DTotalEnergy
 *
 *	Calculate the energy along a 1D coordinate
 */
    double	Minimizer_O::d1DTotalEnergy( double x )
    {
#ifdef	DEBUG_ON
//    this->nvP1DSearchOrigin->debugDump("origin");
//    LOG(BF("x = %lf") % (x ) );
//    this->nvP1DSearchDirection->debugDump("direction");
#endif
	this->getPosition(this->nvP1DSearchTemp1, this->nvP1DSearchOrigin, this->nvP1DSearchDirection,x);
	return this->_EnergyFunction->evaluateEnergy(this->nvP1DSearchTemp1);
    }


/*
 *	d1DTotalEnergyForce
 *
 *	Calculate the energy/derivative along a 1D coordinate
 */
    double	Minimizer_O::d1DTotalEnergyForce( double x, double* fx, double * dfx )
    {
#ifdef	DEBUG_ON
//    this->nvP1DSearchOrigin->debugDump("origin");
//    LOG(BF("x = %lf") % (x ) );
//    this->nvP1DSearchDirection->debugDump("direction");
#endif
	this->getPosition(this->nvP1DSearchTemp1,
			  this->nvP1DSearchOrigin,
			  this->nvP1DSearchDirection, x );
	*fx = this->_EnergyFunction->evaluateEnergyForce( this->nvP1DSearchTemp1,
							  true, this->nvP1DSearchTemp2 );
	*dfx = -(this->nvP1DSearchTemp2)->dotProduct(this->nvP1DSearchDirection);
	return *fx;
    }




/*
 *	minBracket
 *
 *	Bracket the minimum starting at nvOrigin along the
 *	positive direction nvDir.
 * 	Return the bracketing x values ax,bx,cx and the function
 *	evaluations fa,fb,fc at those points.
 *
 *	Start using the values in *ax, *bx;
 */
    void	Minimizer_O::minBracket(
	NVector_sp	nvOrigin,
	NVector_sp	nvDir,
	double		*dPxa,
	double		*dPxb,
	double		*dPxc,
	double		*dPfa,
	double		*dPfb,
	double		*dPfc )
    {_G();
    double		xa,xb,xc,fa,fb,fc, temp;
    double		r,q,u,ulim,fu;

    this->_MinBracketSteps = 0;
    xa = *dPxa;
    xb = *dPxb;
    fa = this->d1DTotalEnergy(xa);
    fb = this->d1DTotalEnergy(xb);
    // Make sure that we are going downhill a->b
    if ( fb > fa ) {
	SWAP_VALUES( xa, xb, temp );
	SWAP_VALUES( fa, fb, temp );
    }
    xc = xb+GOLD*(xb-xa);
    fc = this->d1DTotalEnergy(xc);
    LOG(BF("Start: xa(%lf) xb(%lf) xc(%lf) | fa(%lf) fb(%lf) fc(%lf)") % xa % xb % xc % fa % fb % fc );
    while ( fb > fc ) {
	this->_MinBracketSteps++;
        LOG(BF("Loop:  xa(%lf) xb(%lf) xc(%lf) | fa(%lf) fb(%lf) fc(%lf)")%xa%xb%xc%fa%fb%fc);
	r = (xb-xa)*(fb-fc);
	q = (xb-xc)*(fb-fa);
	u = xb-((xb-xc)*q-(xb-xa)*r)/
	    (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
	ulim = (xb)+GLIMIT*(xc-xb);
	if (( xb-u)*(u-xc)>0.0) {
	    fu = this->d1DTotalEnergy(u);
	    if ( fu<fc ) {
		xa = xb;
		xb = u;
		fa = fb;
		fb = fu;
		goto DONE;
	    } else if ( fu > fb ) {
		xc = u;
		fc = fu;
		goto DONE;
	    }
	    u = xc+GOLD*(xc-xb);
	    fu = this->d1DTotalEnergy(u);
	} else if ((xc-u)*(u-ulim) > 0.0) {
	    fu = this->d1DTotalEnergy(u);
	    if ( fu < fc ) {
		LEFT_SHIFT_VALUES( xb, xc, u, xc+GOLD*(xc-xb) );
		LEFT_SHIFT_VALUES( fb, fc, fu, this->d1DTotalEnergy(u));
	    }
	    // Limit parabolic u to maximum allowed value
	} else if ((u-ulim)*(ulim-xc)>=0.0 ){
	    u = ulim;
	    fu = this->d1DTotalEnergy(u);
	} else {
	    u = xc+GOLD*(xc-xb);
	    fu = this->d1DTotalEnergy(u);
	}
	LEFT_SHIFT_VALUES( xa, xb, xc, u );
	LEFT_SHIFT_VALUES( fa, fb, fc, fu );
    }
    DONE:
#ifdef	DEBUG_ON
    LOG(BF("minBracket  xa(%lf), xb(%lf), xc(%lf)") % xa % xb % xc );
    if ( fa > fb && fb < fc ) {
	LOG(BF("fa.GT.fb.LT.fc") );
    } else {
	LOG(BF("FAIL! minBracket FAILED!! It's not true that: fa.GT.fb.LT.fc") );
    }
    if ( xc < 0 ) {
	LOG(BF("ATTN!  xc is less than zero ") );
    }
#endif
    *dPxa = xa;
    *dPxb = xb;
    *dPxc = xc;
    *dPfa = fa;
    *dPfb = fb;
    *dPfc = fc;
    return;

    }


/*
 *	dbrent
 *
 *	Minimize the function within the interval ax,bx,cx
 *	Return the function value at the new point and the
 *	new minimum in (xmin)
 */
    double 	Minimizer_O::dbrent(	double ax, double bx, double cx,
					double tol,
					double& lineStep,
					int&	energyEvals,
					int&	forceEvals,
					int&	dbrentSteps
	)
    {_G();
    int	iter, ok1, ok2;
    double	_a,_b,_d,d1,d2,du,dv,dw,dx,_e=0.0;
    double	fu,fv,fw,fx,olde,tol1,tol2, u, u1, u2, v,w,x,xm, ft;
    double	retval;

    _d = 0.0;
    energyEvals = 0;
    forceEvals = 0;
    dbrentSteps = 0;
    LOG(BF("Incoming ax bx cx= %lf %lf %lf") % ax % bx % cx );
    _a = (ax<cx?ax:cx);
    _b = (ax>cx?ax:cx);
    LOG(BF("Initial _a _b= %lf %lf") % _a % _b );
    x = w = v = bx;
//
// Calculate the derivative along the search direction
//
    ft = d1DTotalEnergyForce( x, &fx, &dx );
    forceEvals++;
    LOG(BF("dbrent: derivative x,fx,dx = %lf %lf %lf") % x % fx % dx  );
    fw=fv=fx;
    dw=dv=dx;
    for (iter=1;iter<=DBRENTITMAX;iter++ ) {
	dbrentSteps++;
//	LOG(BF("dbrent: iter=%3d  _a, _b= %10.15lf %10.15lf") % iter % _a % _b  );
	xm=0.5*(_a+_b);
	tol1=tol*fabs(x)+ZEPS;
	tol2=2.0*tol1;
	if (fabs(x-xm) <=(tol2-0.5*(_b-_a))) {	// Stopping criterion
	    lineStep=x;
	    retval = fx;
	    LOG(BF("done due to (fabs(x-xm).LE.(tol2-0.5*(_b-_a)))") );
	    goto DONE;
	}
	if (fabs(_e)>tol1){
	    d1=2.0*(_b-_a);
	    d2=d1;
	    if (dw != dx ) d1=(w-x)*dx/(dx-dw); // Secant method, first on one, then on
	    if (dv != dx ) d2=(v-x)*dx/(dx-dv); // the other point
	    // Which of these two estimates of d shall we take? We will insist that
	    // they are within the bracket, and on the side pointed to by the
	    // derivative at x
	    u1=x+d1;
	    u2=x+d2;
	    ok1=(_a-u1)*(u1-_b)>0.0 && dx*d1 <= 0.0;
	    ok2=(_a-u2)*(u2-_b)>0.0 && dx*d2 <= 0.0;
	    olde=_e; // Movement on the lineStep before last
	    _e=_d;
	    if (ok1 || ok2) {	// Take only an acceptable d, and if both are acceptable
				// then take the smallest one.
	        if (ok1 && ok2)
		    _d = (fabs(d1) < fabs(d2) ? d1 : d2 );
		else if (ok1)
		    _d = d1;
		else
		    _d = d2;
		if (fabs(_d) <= fabs(0.5*olde)) {
		    u=x+_d;
		    if ( u-_a<tol2 || _b-u < tol2 )
		        _d = SIGN(tol1,xm-x);
		} else {
		    _e=(dx>=0.0?_a-x:_b-x); // Bisect, not golden section
		    			    // Decide which segment by the sign of the
					    // derivative
		    _d=0.5*_e;
		}
	    } else {
	        _e=(dx>=0.0?_a-x:_b-x);
	        _d=0.5*_e;
	    }
	} else {
	    _e=(dx>=0.0?_a-x:_b-x);
	    _d=0.5*_e;
	}
	if (fabs(_d) >= tol1) {
	    u=x+_d;
	    fu=d1DTotalEnergy(u);
	    energyEvals++;
	} else {
	    u = x+SIGN(tol1,_d);
	    fu=d1DTotalEnergy(u);
	    energyEvals++;
	    if ( fu>fx) { // If the minimum lineStep in the downhill direction takes us
			  // uphill, then we are done
		lineStep=x;
		retval = fx;
		LOG(BF("done due to fu>fx") );
		goto DONE;
	    }
	}
//
// Calculate the force along the search direction
//
	ft = d1DTotalEnergyForce( u, &ft, &du ); // Now housekeeping
	forceEvals++;
//
	if (fu<=fx) {
	    if (u>=x) _a=x; else _b=x;
	    MOV3( v,fv,dv, w,fw,dw );
	    MOV3( w,fw,dw, x,fx,dx );
	    MOV3( x,fx,dx, u,fu,du );
	} else {
	    if (u<x ) _a=u; else _b=u;
	    if (fu<=fw || w==x ) {
		MOV3( v,fv,dv, w,fw,dw );
		MOV3( w,fw,dw, u,fu,du );
	    } else if ( fu<=fv || v == x || v == w ) {
		MOV3( v,fv,dv, u,fu,du );
	    }
	}
    }
    LOG(BF("dbrent: ERROR Exceeded max iterations") );
    retval = -1.0;
    DONE:
    LOG(BF("dbrent evaluated energy(%d) and force(%d) times") % energyEvals % forceEvals  );
    return retval;
    }




    void Minimizer_O::lineSearchInitialReport( StepReport_sp report,
					       NVector_sp nvPos, NVector_sp nvDir, NVector_sp nvForce,
					       double xa, double xb, double xc,
					       double fa, double fb, double fc )
    {
	double lenForce, lenDir, angle, cosAngle;
	lenForce = nvForce->magnitude();
	lenDir = nvDir->magnitude();
	if ( lenForce == 0.0 || lenDir == 0.0 ) {
	    angle = 200.0;
	} else {
	    cosAngle = nvDir->dotProduct(nvForce)/(lenForce*lenDir);
	    if ( cosAngle > 1.0 ) cosAngle = 1.0;
	    if ( cosAngle < -1.0 ) cosAngle = -1.0;
	    angle = acos(cosAngle);
	}
	report->_AngleBetweenDirectionAndForceDeg = angle/0.0174533;
	report->_Xa = xa;
	report->_Xb = xb;
	report->_Xc = xc;
	report->_Fa = fa;
	report->_Fb = fb;
	report->_Fc = fc;
	report->_MinBracketSteps = this->_MinBracketSteps;
	report->_EnergyTermsEnabled = this->_EnergyFunction->energyTermsEnabled();
	report->_TotalEnergy = this->d1DTotalEnergy(0.0);
	report->_DirectionMagnitude = nvDir->magnitude();
	report->_ForceMagnitude = nvForce->magnitude();
	report->_MinimizerStatus = this->statusAsString();
	double dxa,dxc;
	if ( xa < xc ) {
	    dxa = xa;
	    dxc = xc;
	    report->_Direction = "searchForward";
	} else {
	    dxa = xc;
	    dxc = xa;
	    report->_Direction = "searchBackward";
	}
	double zx, zy;
	double xmin = dxa;
	double xinc = (dxc-dxa)/100.0;
	report->_LineSearchPosition = nvPos->deepCopy().as<NVector_O>();
	report->_LineSearchDirection = nvDir->deepCopy().as<NVector_O>();
	report->_TotalEnergyFn = NumericalFunction_O::create("Alpha","Total",xmin,xinc,_lisp);
	report->_StretchEnergyFn = NumericalFunction_O::create("Alpha","Stretch",xmin,xinc,_lisp);
	report->_AngleEnergyFn = NumericalFunction_O::create("Alpha","Angle",xmin,xinc,_lisp);
	report->_DihedralEnergyFn = NumericalFunction_O::create("Alpha","Dihedral",xmin,xinc,_lisp);
	report->_NonbondEnergyFn = NumericalFunction_O::create("Alpha","Nonbond",xmin,xinc,_lisp);
	report->_ImproperEnergyFn = NumericalFunction_O::create("Alpha","Improper",xmin,xinc,_lisp);
	report->_ChiralRestraintEnergyFn = NumericalFunction_O::create("Alpha","ChiralRestraint",xmin,xinc,_lisp);
	report->_AnchorRestraintEnergyFn = NumericalFunction_O::create("Alpha","AnchorRestraint",xmin,xinc,_lisp);
	report->_FixedNonbondRestraintEnergyFn = NumericalFunction_O::create("Alpha","FixedNonbondRestraint",xmin,xinc,_lisp);

	for ( zx=dxa;zx<=dxc;zx+=(dxc-dxa)/100.0 ) {
	    zy = this->d1DTotalEnergy(zx);
	    report->_TotalEnergyFn->appendValue(zy);
	    report->_StretchEnergyFn->appendValue(
		this->_EnergyFunction->getStretchComponent()->getEnergy());
#if USE_ALL_ENERGY_COMPONENTS
	    report->_AngleEnergyFn->appendValue(
		this->_EnergyFunction->getAngleComponent()->getEnergy());
	    report->_DihedralEnergyFn->appendValue(
		this->_EnergyFunction->getDihedralComponent()->getEnergy());
	    report->_NonbondEnergyFn->appendValue(
		this->_EnergyFunction->getNonbondComponent()->getEnergy());
	    report->_ImproperEnergyFn->appendValue(
		this->_EnergyFunction->getImproperRestraintComponent()->getEnergy());
	    report->_ChiralRestraintEnergyFn->appendValue(
		this->_EnergyFunction->getChiralRestraintComponent()->getEnergy());
	    report->_AnchorRestraintEnergyFn->appendValue(
		this->_EnergyFunction->getAnchorRestraintComponent()->getEnergy());
	    report->_FixedNonbondRestraintEnergyFn->appendValue(
		this->_EnergyFunction->getFixedNonbondRestraintComponent()->getEnergy());
#endif
	}
    };


    void Minimizer_O::throwMinimizerExceededMaxSteps()
    {_OF();
	MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(BF("test throw of MinimizerExceededMaxSteps"));
    };

    void Minimizer_O::throwMinimizerStuck()
    {_OF();
	MINIMIZER_STUCK_ERROR(BF("test throw of MinimizerStuck"));
    };


    void Minimizer_O::throwMinimizerError()
    {_OF();
	MINIMIZER_ERROR(BF("test throw of MinimizerError"));
    };



    void Minimizer_O::lineSearchFinalReport( StepReport_sp report, double step, double fMin,
					     int energyEvals, int forceEvals, int dbrentSteps )
    {_OF();
    ASSERT(report->_Iteration == this->_Iteration);
    report->_DbrentSteps = dbrentSteps;
    report->_EnergyEvals = energyEvals;
    report->_ForceEvals = forceEvals;
    report->_Step = step;
    report->_FMin = fMin;
    }


    void	Minimizer_O::stepReport( StepReport_sp report, double energy, NVector_sp force )
    {_OF();
    ASSERT(report->_Iteration == this->_Iteration);
    report->_ForceMagnitude = force->magnitude();
    report->_TotalEnergy = energy;
    report->_IterationMessages = this->_IterationMessages.str();
    }

/*!
 *      Perform a line search of the amber function
 *	From the point nvOrigin along the positive direction nvDir.
 *
 *	Only *dPstep and *dPfnew are updated.
 */
    void	Minimizer_O::lineSearch(	double	*dPstep,
						double	*dPfnew,
						NVector_sp nvOrigin,
						NVector_sp nvDirection,
						NVector_sp nvForce,
						NVector_sp nvTemp1,
						NVector_sp nvTemp2,
						int	iteration,
						StepReport_sp	report )
    {_G();
    double	xa, xb, xc;
    double	fa, fb, fc;
    double	step = 0.0;
    int	functionEvals;

    functionEvals = 0;
    LOG(BF("Starting") );

    //
    // Setup the 1D search origin and direction
    //

    this->define1DSearch(nvOrigin,nvDirection,nvTemp1,nvTemp2);
    double directionMag = nvDirection->magnitude();
    if ( directionMag < VERYSMALLSQUARED ) {
	xb = this->_InitialLineSearchStep;
    } else {
	xb = this->_InitialLineSearchStep/directionMag;
    }
    xa = 0.0;

    //
    // Bracket the minimum
    //

    this->minBracket( nvOrigin, nvDirection,
		      &xa, &xb, &xc, &fa, &fb, &fc );
    LOG(BF("Bracketed minimum") );
    LOG(BF("xa,xb,xc = %lf %lf %lf ") % xa % xb % xc  );
    LOG(BF("fa,fb,fc = %lf %lf %lf ") % fa % fb % fc  );
    if ( this->_DebugOn )
    {
	this->lineSearchInitialReport(report,nvOrigin,nvDirection,nvForce,
				      xa,xb,xc,fa,fb,fc);
    }

    //
    // Now use the bracketed minimum to find the line minimum
    //
    //
    // THIS IS WHERE YOU ADD DBRENT
    int energyEvals = 0;
    int forceEvals = 0;
    int dbrentSteps = 0;
    fb = dbrent( xa, xb, xc, TOL, step, energyEvals, forceEvals, dbrentSteps );
    LOG(BF("brent bracketed step = %lf") % step  );
    if ( this->_DebugOn )
    {
	this->lineSearchFinalReport( report, step, fb, energyEvals,
				     forceEvals, dbrentSteps );
    }

    *dPstep = step;
    *dPfnew = fb;


    LOG(BF("Done") );
    }



    void	Minimizer_O::_displayIntermediateMessage(
	double			step,
	double			fnew,
	double			forceMag,
	double			forceRMSMag,
	double			cosAngle,
	bool			steepestDescent )
    {
#define	MAX_DISPLAY	1024
	double		angle;
	char		buffer[MAX_DISPLAY];
	stringstream	sout;
	if ( this->_Iteration%10 == 1 || this->_DebugOn ) 
	{
	    sout << "---Stage--";
	    if ( this->_ShowElapsedTime )
	    {
		sout << "Seconds-";
	    }
	    sout << "Step-----Alpha---Dir-------------Energy-----------RMSforce";
	    if ( this->_EnergyFunction->getName() != "" ) 
	    {
		sout << "-------Name";
	    }
	    sout << std::endl;
	}
	sout << BF(" min%4s") % this->statusAsShortString();
	if ( this->_ShowElapsedTime )
	{
#if USE_POSIX_TIME
	    core::PosixTimeDuration_sp elapsed = core::PosixTimeDuration_O::createDurationSince(this->_StartTime);
	    sout << BF(" %8ld") % elapsed->totalSeconds();
#endif
	}
	sout << BF(" %5d") % this->_Iteration;
	sout << BF(" %9.6lf") % step;
	if ( steepestDescent ) 
	{
	    sout << "StDesc";
	} else 
	{
	    if ( cosAngle < -1.0 ) cosAngle = -1.0;
	    if ( cosAngle > 1.0 ) cosAngle = 1.0;
	    angle = acos(cosAngle)/0.0174533;
	    if ( angle < 0.1 ) 
	    {
		angle = 0.0;
	    }
	    sout << (BF(" %5.1lf") % angle);
	}
	sout << (BF(" %18.8lf") % fnew );
	sout << (BF(" %18.8lf") % forceRMSMag);
	if ( this->_EnergyFunction->getName() != "" ) 
	{
	    sout << BF(" %s") % this->_EnergyFunction->getName();
	}
	_lisp->print(BF("%s") % sout.str() );
	if ( this->_DebugOn ) 
	{
	    this->_Log->addMessage(buffer);
	}
    }




/*!
 * Steepest descent minimizer with or without preconditioning
 */
    void	Minimizer_O::_steepestDescent( int numSteps,
					       NVector_sp x,
					       double forceTolerance )
    {_G();
	StepReport_sp	stepReport = _Nil<StepReport_O>();
	int		iRestartSteps;
	NVector_sp	force, m;
	NVector_sp	s,dir,tv1,tv2;
	double		forceMag, forceRmsMag,prevStep;
	double		delta0, deltaNew;
	double		eSquaredDelta0;
	double		step, fnew, dirMag;
	double		cosAngle = 0.0;
	bool		steepestDescent;
	int		innerSteps;
	int		localSteps, k;

	LOG(BF("Checking status") );
	if ( this->_Status == minimizerError ) return;
	this->_Status = steepestDescentRunning;
	this->_CurrentPreconditioner = noPreconditioner;
	LOG(BF("step") );

	localSteps = 0;
	k = 0;
	step = 0.0;

	LOG(BF("step") );
	/* Calculate how many conjugate gradient steps can be */
	/* taken before a restart must be done */

	iRestartSteps = x->size();
	// Define NVectors
	force = NVector_O::create(iRestartSteps,_lisp);
	s = NVector_O::create(iRestartSteps,_lisp);
	LOG(BF("step") );
	dir = NVector_O::create(iRestartSteps,_lisp);
	tv1 = NVector_O::create(iRestartSteps,_lisp);
	tv2 = NVector_O::create(iRestartSteps,_lisp);
	LOG(BF("step") );
	// Done
	innerSteps = MIN(iRestartSteps,ITMAX);
	LOG(BF("step") );
	double fp = this->dTotalEnergyForce( x, force );
	LOG(BF("step") );
//    r->inPlaceTimesScalar(-1.0);
	//  no preconditioning
	s->copyIntoMe(force);
	LOG(BF("Done initialization") );
#if 0 //[
	// TODO calculate preconditioner here
	// s = M^(-1)force rather than just copying it from force
	switch ( preconditioner )
	{
	case noPreconditioner:
	    s->copy(force);
	    break;
	case identityPreconditioner:
	    m = NVector_O::create(iRestartSteps);
	    m->fill(1.0);
	    this->_EnergyFunction->backSubstituteDiagonalPreconditioner(m,s,force);
	    break;
	case diagonalPreconditioner:
	    m = NVector_O::create(iRestartSteps);
	    this->_EnergyFunction->setupDiagonalPreconditioner(x,m);
	    LOG(BF("Preconditioner max value: %lf") % m->maxValue() );
	    LOG(BF("Preconditioner min value: %lf") % m->minValue() );
	    minVal = m->minValue();
	    if ( minVal < 0.0 ) {
		m->addScalar(m,fabs(minVal)+1.0);
	    }
	    this->_EnergyFunction->backSubstituteDiagonalPreconditioner(m,s,force);
	    break;
	default:
	    SIMPLE_ERROR(BF("Unsupported preconditioner"));
	}
#endif //]
	dir->copyIntoMe(s);
	deltaNew = force->dotProduct(dir);
	delta0 = deltaNew;
	eSquaredDelta0 = forceTolerance*delta0;
	LOG(BF("eSquaredDelta0 = %lf") % (eSquaredDelta0 ) );
	LOG(BF("forceTolerance = %lf") % (forceTolerance ) );
	LOG(BF("delta0 = %lf") % (delta0 ) );
	localSteps = 0;
	if ( this->_PrintIntermediateResults ) {
	    _lisp->print(BF( "======= Starting Steepest Descent Minimizer" ));
	}
	try
	{
	    while (1) 
	    { _BLOCK_TRACEF(BF("Step %d") %localSteps);

		//
		// Absolute gradient test
		//
		forceMag = force->magnitude();
		forceRmsMag = force->rmsMagnitude();
		this->_RMSForce = forceRmsMag;
		if ( forceRmsMag < forceTolerance ) {
		    if ( this->_PrintIntermediateResults ) {
			_lisp->print(BF("DONE absolute force test:\nforceRmsMag(%lf).LT.forceTolerance(%lf)") % forceRmsMag % forceTolerance );
		    }
		    break;
		}

		this->_IterationMessages.str("");
		if ( localSteps>=numSteps ) {
		    MinimizerCondition_ExceededNumSteps fail;
		    fail._Message = "steepest descent";
		    fail._NumberOfSteps = numSteps;
		    throw(fail);
		}


//
// Here we need to implement the three Progress tests
//
//  See MOE manual "Optimization Functions"
// Stop when the following tests are all satisfied simultaneously:
// 1. F(x(k-1)) - F(xk) < T (1 + |F(xk)|)
// 2. |x(k-1) - xk| < Sqrt(T)(1+|xk|)
// 3. |grad F(xk)| <= T^(1/3)(1+|F(xk)|)
//
// But I have a problem with the CG optimizer ending while there is still a substantial gradient.
// That may be a problem with the potential function and I may need to solve this by
// validating the energy function externally
//
//




		// ***************************************************
		// ***************************************************
		// ***************************************************
		//
		// Once all the tests are complete then we start the iteration in earnest
		//
		{ // Dont break out of this block unless you THROW a SERIOUS error
		    if ( this->_DebugOn )
		    {
			stepReport = StepReport_O::create();
			stepReport->_Iteration = this->_Iteration;
		    }

		    //
		    // Descent test
		    // If the cos(angle) between the search direction and
		    // the force (steepest descent dir)
		    // is less than _xxxDescentTest then copy the force into the search dir
		    //
		    dirMag = dir->magnitude();
		    steepestDescent = false;
		    LOG(BF("Starting descent test") );
		    if ( forceMag != 0.0 && dirMag != 0.0 ) {
			cosAngle = force->dotProduct(dir)/(forceMag*dirMag);
		    } else {
			LOG(BF("something was zero length Using force") );
			dir->copyIntoMe(force);
			steepestDescent = true;
		    }

		    //
		    // Print intermediate status
		    //
		    prevStep = step;
		    if ( this->_PrintIntermediateResults ) {
			this->_displayIntermediateMessage(step,fnew,forceMag,forceRmsMag,cosAngle,steepestDescent);
		    }

		    this->lineSearch( &step, &fnew, x, dir, force, tv1, tv2, localSteps, stepReport );


		    if ( prevStep == 0.0 && step == 0.0 ) {
			MinimizerCondition_Stuck fail("Two successive CG steps went nowhere");
			throw(fail);
		    }

		    x->inPlaceAddTimesScalar( dir, step );

		    // r = -f'(x)   r == force!!!!
		    fp = dTotalEnergyForce( x, force );

#if 0 //[
		    if ( this->_EnergyFunction->getMessages() != "" )
		    {
			_lisp->print(BF("Minimizer_O::steepestDescentWithPreconditioner>> energy function messages: %s")%
					    this->_EnergyFunction->getMessages());
			this->_IterationMessages << this->_EnergyFunction->getMessages();
		    }
#endif //]

		    if ( this->_DebugOn ) {
			this->stepReport(stepReport,fp,force);
		    }


#ifdef	VALIDATE_FORCE_ON //[
		    this->validateForce(x,force);
#endif //]


		    // Don't use preconditioning
		    s->copyIntoMe(force);

#if 0 //[
		    switch ( preconditioner ) {
		    case noPreconditioner:
			s->copy(force);
			break;
		    case identityPreconditioner:
			m->fill(1.0);
			this->_EnergyFunction->backSubstituteDiagonalPreconditioner(m,s,force);
			break;
		    case diagonalPreconditioner:
			LOG(BF("Calculating diagonal preconditioner") );
			this->_EnergyFunction->setupDiagonalPreconditioner(x,m);
			LOG(BF("Preconditioner max value: %lf") % m->maxValue()  );
			LOG(BF("Preconditioner min value: %lf") % m->minValue()  );
			minVal = m->minValue();
			if ( minVal < 0.0 ) {
			    m->addScalar(m,fabs(minVal)+1.0);
			}
			this->_EnergyFunction->backSubstituteDiagonalPreconditioner(m,s,force);
			//		    if ( s->dotProduct(force) < 0 ) {
			//			s->copy(force);
			//		    }
			break;
		    default:
			SIMPLE_ERROR(BF("Unsupported preconditioner"));
			break;
		    }
#endif //]
		    dir->copyIntoMe(s);
		    if ( this->_DebugOn )
		    {
			ASSERTNOTNULL(this->_Log);
			this->_Log->addReport(stepReport);
			stepReport = _Nil<StepReport_O>();
		    }
		    localSteps++;
		    this->_Iteration++;
		}
#if 0
		if ( this->_StepCallback.notnilp() )
		{
		    this->_StepCallback->invoke(this->sharedThis<Minimizer_O>());
		}
#endif
	    }
	} catch (MinimizerCondition_ExceededNumSteps fail) {
	    if ( this->_DebugOn )
	    {
		ANN(stepReport);
		if ( stepReport.notnilp() )
		{
		    stepReport->prematureTermination("ExceededNumSteps");
		    this->_Log->addReport(stepReport);
		    stepReport = _Nil<StepReport_O>();
		}
	    }

	    //
	    // Lets save the current conformation
	    // before throwing this higher
	    //
	    fp = dTotalEnergyForce( x, force );
	    this->_EnergyFunction->writeCoordinatesAndForceToAtoms(x,force);
	    MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(fail.message());
	}
	fp = dTotalEnergyForce( x, force );
	this->_EnergyFunction->writeCoordinatesAndForceToAtoms(x,force);
	LOG(BF("Wrote coordinates and force to atoms") );
	if ( this->_DebugOn )
	{
	    LOG(BF("About to update stepReport") );
	    ASSERTNOTNULL(stepReport);
	    if ( stepReport.notnilp() )
	    {
		LOG(BF("About to add to stepReport") );
		stepReport->prematureTermination("Stuck");
		ASSERTNOTNULL(this->_Log);
		LOG(BF("About to write report to debugLog") );
		this->_Log->addReport(stepReport);
		LOG(BF("Done writing report to debugLog") );
		stepReport = _Nil<StepReport_O>();
	    }
	}
	LOG(BF("Leaving _steepestDescent") );
    }



    void	Minimizer_O::_conjugateGradient(
	int numSteps,
	NVector_sp x,
	double forceTolerance )
    {_G();
	StepReport_sp	stepReport = _Nil<StepReport_O>();
    int		iRestartSteps;
    NVector_sp	force, diag;
    SparseLargeSquareMatrix_sp	m,ldlt;
    NVector_sp	s,d,tv1,tv2;
    double		forceMag;
    double		prevStep = 0.0;
    double		forceRmsMag;
    double		delta0, deltaNew, deltaMid, deltaOld;
    double		eSquaredDelta0;
    double		step, fnew ;
    int		innerSteps;
    int		localSteps, k;
    double		beta, cosAngle, dirMag;
    bool		steepestDescent;
    int		refactor;

    if ( this->_Status == minimizerError ) return;
    this->_Status = conjugateGradientRunning;
    this->_CurrentPreconditioner = noPreconditioner;
    localSteps = 0;
    k = 0;
    refactor = 0;

    /* Calculate how many conjugate gradient steps can be */
    /* taken before a restart must be done */

    iRestartSteps = x->size();
    // Define NVectors
    force = NVector_O::create(iRestartSteps,_lisp);
    s = NVector_O::create(iRestartSteps,_lisp);
    d = NVector_O::create(iRestartSteps,_lisp);
    tv1 = NVector_O::create(iRestartSteps,_lisp);
    tv2 = NVector_O::create(iRestartSteps,_lisp);
    // Done
    innerSteps = MIN(iRestartSteps,ITMAX);
    double fp = dTotalEnergyForce( x, force );
//    r->inPlaceTimesScalar(-1.0);
    // TODO calculate preconditioner here
    // s = M^(-1)r rather than just copying it from r
    s->copyIntoMe(force);
#if 0 //[
    switch ( preconditioner ) {
    case noPreconditioner:
	s->copy(force);
	break;
    case diagonalPreconditioner:
	diag = NVector_O::create(iRestartSteps,_lisp);
	this->_EnergyFunction->setupDiagonalPreconditioner(x,diag);
	LOG(BF("Preconditioner max value: %lf") % diag->maxValue() );
	LOG(BF("Preconditioner min value: %lf") % diag->minValue() );
	this->_EnergyFunction->backSubstituteDiagonalPreconditioner(diag,s,force);
	break;
    case hessianPreconditioner:
	m = new_SparseLargeSquareMatrix_sp(iRestartSteps,SymmetricUpperDiagonal);
	ldlt=new_SparseLargeSquareMatrix_sp(iRestartSteps,SymmetricUpperDiagonal);
	m->fill(0.0);
	ldlt->fill(0.0);
	this->_EnergyFunction->setupHessianPreconditioner(x,m);
	this->_EnergyFunction->unconventionalModifiedCholeskyFactorization(m,ldlt);
	this->_EnergyFunction->backSubstituteLDLt(ldlt,s,force);
	break;
    default:
	SIMPLE_ERROR(BF("Unknown preconditioner option"));
    }
#endif //]
    d->copyIntoMe(s);
    deltaNew = force->dotProduct(d);
    delta0 = deltaNew;
    eSquaredDelta0 = forceTolerance*delta0;
    LOG(BF("eSquaredDelta0 = %lf") % (eSquaredDelta0 ) );
    LOG(BF("forceTolerance = %lf") % (forceTolerance ) );
    LOG(BF("delta0 = %lf") % (delta0 ) );
    localSteps = 0;
    step = 0.0;
    if ( this->_PrintIntermediateResults ) {
	_lisp->print(BF( "======= Starting Conjugate Gradient Minimizer" ));
    }
    try {
	while (1) {


	    //
	    // Absolute gradient test
	    //
	    forceMag = force->magnitude();
	    forceRmsMag = force->rmsMagnitude();
	    this->_RMSForce = forceRmsMag;
	    if ( forceRmsMag < forceTolerance ) {
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF("DONE absolute force test:\nforceRmsMag(%lf).LT.forceTolerance(%lf)")% forceRmsMag % forceTolerance );
		}
		break;
	    }

	    this->_IterationMessages.str("");

	    if ( localSteps > 0 ) {
		if ( localSteps>=numSteps ) {
		    if ( this->_PrintIntermediateResults )
		    {
			_lisp->print(BF( "Exceeded max number of steps(%d)") % numSteps );
		    }
		    MinimizerCondition_ExceededNumSteps fail;
		    fail._Message = "conjugate gradient";
		    fail._NumberOfSteps = numSteps;
		    throw(fail);
		}

		if ( prevStep == 0.0 && step == 0.0 ) {
		    MinimizerCondition_Stuck fail("CG-Two successive steps went nowhere!!!!");
		    throw(fail);
		}
	    }

//
// Here we need to implement the three Progress tests
//
//  See MOE manual "Optimization Functions"
// Stop when the following tests are all satisfied simultaneously:
// 1. F(x(k-1)) - F(xk) < T (1 + |F(xk)|)
// 2. |x(k-1) - xk| < Sqrt(T)(1+|xk|)
// 3. |grad F(xk)| <= T^(1/3)(1+|F(xk)|)
//
// But I have a problem with the CG optimizer ending while there is still a substantial gradient.
// That may be a problem with the potential function and I may need to solve this by
// validating the energy function externally
//
//



	    // ***************************************************
	    // ***************************************************
	    // ***************************************************
	    //
	    // Once all the tests are complete then we start the iteration in earnest
	    //
	    {  // Dont break out of this block unless its a SERIOUS ERROR
		if ( this->_DebugOn )
		{
		    stepReport = StepReport_O::create();
		    stepReport->_Iteration = this->_Iteration;
		}

		// j = 0
		// deltaD = d.d
		//	deltaD = d->squared();	// Used for secant method
		//
		//
		// Descent test
		// If the cos(angle) between the search direction and
		// the force (steepest descent dir)
		// is less than _xxxDescentTest then copy the force into the search dir
		//
		dirMag = d->magnitude();
		steepestDescent = false;
		LOG(BF("Starting descent test") );
		cosAngle = 0.0;
		if ( forceMag != 0.0 && dirMag != 0.0 ) {
		    LOG(BF("forceMag = %lf") % forceMag  );
		    LOG(BF("dirMag = %lf") % dirMag  );
		    cosAngle = force->dotProduct(d)/(forceMag*dirMag);
		} else {
		    LOG(BF("some magnitude was zero Using force") );
		    d->copyIntoMe(force);
		    steepestDescent = true;
		}

		prevStep = step;
		if ( this->_PrintIntermediateResults ) {
		    this->_displayIntermediateMessage(prevStep,fnew,forceMag,forceRmsMag,cosAngle,steepestDescent);
		}

		this->lineSearch( &step, &fnew, x, d, force,
				  tv1, tv2, localSteps,
				  stepReport);


		// x = x + (step)d
		// r = -f'(x)   r == force!!!!
		x->inPlaceAddTimesScalar( d, step );
		fp = dTotalEnergyForce( x, force );


		if ( this->_DebugOn )
		{
		    this->stepReport(stepReport,fp,force);
		}


#ifdef	VALIDATE_FORCE_ON
		this->validateForce(x,force);
#endif
		deltaOld = deltaNew;
		deltaMid = force->dotProduct(s);

		// No preconditioning
		s->copyIntoMe(force);

#if 0 //[
		// Calculate preconditioner M = f''(x)
		// s = M^(-1)r
		switch ( preconditioner ) {
		case noPreconditioner:
		    s->copy(force);
		    break;
		case diagonalPreconditioner:
		    this->_EnergyFunction->setupDiagonalPreconditioner(x,diag);
		    LOG(BF("Preconditioner max value: %lf") % diag->maxValue() );
		    LOG(BF("Preconditioner min value: %lf") % diag->minValue() );
		    this->_EnergyFunction->backSubstituteDiagonalPreconditioner(diag,s,force);
		    break;
		case hessianPreconditioner:
		    //		refactor++;
		    //		if ( refactor >= 5 ) {
		    this->_EnergyFunction->setupHessianPreconditioner(x,m);
		    this->_EnergyFunction->unconventionalModifiedCholeskyFactorization(m,ldlt);
		    refactor = 0;
		    //		}
		    this->_EnergyFunction->backSubstituteLDLt(ldlt,s,force);
		    break;
		default:
		    SIMPLE_ERROR(BF("Unknown preconditioner option"));
		}
#endif //]
		deltaNew = force->dotProduct(s);		// deltaNew = r.r
		beta = (deltaNew-deltaMid)/deltaOld;
		k = k + 1;
		if ( k == iRestartSteps || beta <= 0.0 ) {
		    d->copyIntoMe(s);
		    k = 0;
		    prevStep = 0.0;
		} else {
		    d->XPlusYTimesScalar(s,d,beta);
		}
		if ( this->_DebugOn )
		{
		    ASSERTNOTNULL(this->_Log);
		    this->_Log->addReport(stepReport);
		    stepReport = _Nil<StepReport_O>();
		}
		localSteps++;
		this->_Iteration++;
	    }
#if 0
	    if ( this->_StepCallback.notnilp() )
		this->_StepCallback->invoke(this->sharedThis<Minimizer_O>());
#endif
	}
    } catch (MinimizerCondition_ExceededNumSteps fail) {
	if ( this->_DebugOn )
	{
	    if ( stepReport.notnilp() )
	    {
		stepReport->prematureTermination("ExceededNumSteps");
		this->_Log->addReport(stepReport);
		stepReport = _Nil<StepReport_O>();
	    }
	}
	//
	// Lets save the current conformation
	// before throwing this higher
	//
	fp = dTotalEnergyForce( x, force );
	this->_EnergyFunction->writeCoordinatesAndForceToAtoms(x,force);
	MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(fail.message());
    } catch (MinimizerCondition_Stuck fail) {
	if ( this->_DebugOn )
	{
	    if ( stepReport.notnilp() )
	    {
		stepReport->prematureTermination("Stuck");
		this->_Log->addReport(stepReport);
		stepReport = _Nil<StepReport_O>();
	    }
	}
	//
	// Lets save the current conformation
	// before throwing this higher
	//
	fp = dTotalEnergyForce( x, force );
	this->_EnergyFunction->writeCoordinatesAndForceToAtoms(x,force);
	MINIMIZER_STUCK_ERROR(fail.message());
    }
    fp = dTotalEnergyForce( x, force );
    this->_EnergyFunction->writeCoordinatesAndForceToAtoms(x,force);
    if ( this->_DebugOn )
    {
	if ( stepReport.notnilp() )
	{
	    stepReport->prematureTermination("Stuck");
	    this->_Log->addReport(stepReport);
	    stepReport = _Nil<StepReport_O>();
	}
    }
    }





    void	Minimizer_O::_truncatedNewtonInnerLoop(
	int				kk,
	NVector_sp			xk,
	SparseLargeSquareMatrix_sp	mprecon,
	SparseLargeSquareMatrix_sp	ldlt,
	NVector_sp 			force,
	double				rmsForceMag,
	NVector_sp			pj,
	NVector_sp			pjNext,
	NVector_sp			rj,
	NVector_sp			dj,
	NVector_sp			zj,
	NVector_sp			qj )
    {_G();
    int	j, ITpcg;
    double	cr, delta, forceDotpj, crOverk, nk, alphaj;
    double	nkTimesRmsForceMag, rjDotzj;
    double				djDotqj, forceDotpjNext;
    double				rmsRjMag, rjDotzjNext, betaj;
    NVector_sp			nvDummy;
    SparseLargeSquareMatrix_sp	nmDummy;

    ASSERTNOTNULL(this->_EnergyFunction);
    LOG(BF("Resetting dummy vector and matrix") );
    nvDummy = _Nil<NVector_O>();
    nmDummy = _Nil<SparseLargeSquareMatrix_O>();

    if ( this->_DebugOn )
    {
	this->_Log->addMessage("_truncatedNewtonInnerLoop>>Starting\n");
    }
    LOG(BF("Setting up") );
    j = 1;
    pj->zero();	// NVector
    cr = 0.5;
    delta = 10.0e-10;
    forceDotpj = 0.0;	// The initial value of force.pj
    ITpcg = 40;
    rj->copyIntoMe(force);
    ASSERT(kk>0);
    crOverk = cr/((float)(kk));
    nk = MIN(crOverk,rmsForceMag);
    nkTimesRmsForceMag = nk*rmsForceMag;
    // 2a.
    // Perform the UMC of M so that the resulting effective
    // preconditioner is M~=LDL^t with a chosen parameter tao
    // (schlick uses tao=0). The factor L is stored in the
    // same sparse row format used for M.
    //
    LOG(BF("Carrying out UMC") );
//
//  Carry out UMC in the outer loop
//    this->_EnergyFunction->unconventionalModifiedCholeskyFactorization(mprecon,ldlt);
//
    // 2b.
    // Solve for zj in (M~)(zj)=(rj) by using the triangular
    // systems: (L)(x)=(rj) and (L^T)(zj)=(D^-1)(x)
    //
    LOG(BF("Back substitute") );
    backSubstituteLDLt(ldlt,zj,rj);

    // 2c.
    // set dj = zj
    //
    dj->copyIntoMe(zj);

    rjDotzj = rj->dotProduct(zj);
    while ( 1 ) 
    { _BLOCK_TRACEF(BF("_truncatedNewtonInnerLoop j=%d") % j );
    // 3. Singularity test
    // Compute the matrix-vector product qj=(H)(dj)
    // If either |((rj)^T).(zj)|<=delta
    // or |((dj)^T).(qj)|<=delta (e.g., delta=10^-10)
    // exit PCG loop with pk=pj ( for j=1, set pk=force)
    //

    this->_EnergyFunction->evaluateAll( xk, true, nvDummy,
					true, true, nmDummy,
					qj, dj );
    // MOVE rjDotzj calculation above this loop because
    // 	its calculated in step 6
    // rjDotzj = rj->dotProduct(zj);
    djDotqj = dj->dotProduct(qj);
    if ( fabs(rjDotzj) <= delta || fabs(djDotqj) <= delta ) {
	if ( j==1 ) {
	    pj->copyIntoMe(force);
	}
	if ( this->_DebugOn ) {
	    stringstream ss;
	    ss << "rjDotzj("<<rjDotzj<<").LT.delta("<<delta<<") or djDotqj("<<djDotqj<<").LT.delta("<<delta<<")"<<std::endl;
	    this->_Log->addMessage(ss.str().c_str());
	    this->_Log->addMessage("_truncatedNewtonInnerLoop>>Singularity test was true\n" );
	}
	LOG(BF("Singularity test was true") );
	goto DONE;
    }

    // 4. Implement the descent direction test and not
    // 	the negative curvature test.
    // [Descent direction test]
    //    Update the quantities
    //      alphaj = rjDotzj/djDotqj;
    //	{pj+1}->pj + alphaj*dj
    // If (gk.{pj+1}) >= (gk.pj) + delta then {
    //	exit inner loop with pk = pj (for j = 1, set pk = force )
    // }
    //   Note: since I use force rather than gradient(gk)
    //      I'm pretty sure that I need to invert the inequality test.
    //
    alphaj = rjDotzj/djDotqj;
    pjNext->XPlusYTimesScalar(pj,dj,alphaj);
    LOG(BF("pjNext angle with force=%lf(deg)") % pjNext->angleWithVector(force)/0.0174533 );
    forceDotpjNext = force->dotProduct(pjNext);
    if ( forceDotpjNext <= (forceDotpj + delta) ) {
	if ( j == 1 ) {
	    pj->copyIntoMe(force);
	} else {
	    // pk->copy(pj);  pk is pj
	}
	if ( this->_DebugOn ) {
	    this->_Log->addMessage("_truncatedNewtonInnerLoop>>Descent direction test was true\n" );
	}
	LOG(BF("Descent direction test was true") );
	goto DONE;
    }

    // 5. Truncation test
    // Compute rjNext = rj - (alphaj)(qj)
    // If "Truncation test" ||rjNext||<nk||g||
    // 	or "Step Limit" j+1>ITpcg then {
    // 	exit inner loop with search direction pk = pjNext
    // }
    //
    rj->inPlaceAddTimesScalar(qj,-alphaj);
    rmsRjMag = rj->rmsMagnitude();
    if ( rmsRjMag < nkTimesRmsForceMag || (j+1)>ITpcg ) {
	LOG(BF("rmsRjMag(%lf) < nkTimesRmsForceMag(%lf)") % rmsRjMag % nkTimesRmsForceMag );
	pj->copyIntoMe(pjNext);
	if ( this->_DebugOn ) {
	    this->_Log->addMessage("_truncatedNewtonInnerLoop>>Truncation test was true\n" );
	}
	LOG(BF("Truncation test was true") );
	goto DONE;
    }
    if ( (j+1)>ITpcg ) {
	LOG(BF("j+1(%d)>ITpcg(%d)") % j+1 % ITpcg );
	pj->copyIntoMe(pjNext);
	if ( this->_DebugOn ) {
	    this->_Log->addMessage("_truncatedNewtonInnerLoop>>Step limit test was true\n" );
	}
	LOG(BF("Step limit test was true") );
	goto DONE;
    }

    // 6. Continuation of PCG
    // Solve for zjNext as in step 2 in (M~)zjNext = rj
    // Update the quantities
    // 	betaj = (rjNext.zjNext)/(rj.zj)
    // 	djNext = zjNext + betaj*dj
    // j = j + 1 goto step 3
    //
    backSubstituteLDLt(ldlt,zj,rj);
    rjDotzjNext = rj->dotProduct(zj);
    betaj = rjDotzjNext/rjDotzj;
    rjDotzj = rjDotzjNext;
    dj->XPlusYTimesScalar(zj,dj,betaj);
    j = j + 1;
    pj->copyIntoMe(pjNext);
    }
    DONE:
    LOG(BF("Exiting inner loop with j = %d") % j );
    return;
    }



#define	EPSILONF	1.0e-10
#define	EPSILONG	1.0e-8
#define	SQRT_EPSILONF	1.0e-5
#define	CUBERT_EPSILONF	4.6416e-4

    void	Minimizer_O::_truncatedNewton(
	int numSteps,
	NVector_sp xK,
	double forceTolerance )
    {_G();
	StepReport_sp		stepReport = _Nil<StepReport_O>();
    int	iDimensions;
    double			fp;
    NVector_sp	forceK, pK, pjNext, rj, dj, zj, qj, xKNext, kSum;
    SparseLargeSquareMatrix_sp	mprecon, ldlt;
    double				energyXk, energyXkNext;
    double				rmsForceMag;
    int				kk;
    double				alphaK, delta, rmsMagXKNext;
    bool				b1aTest, b1bTest;
    double		prevAlphaK, dirMag, forceMag, cosAngle;
#define	TENEMINUS8	10.0e-8

    alphaK = 1.0;
    LOG(BF("Starting TN") );
    if ( this->_Status == minimizerError ) return;
    this->_Status = truncatedNewtonRunning;
    this->_CurrentPreconditioner = this->_TruncatedNewtonPreconditioner;
    kk = 1;
    // Define NVectors
    LOG(BF("Defining NVectors") );
    iDimensions = xK->size();
    forceK = NVector_O::create(iDimensions,_lisp);
    LOG(BF("Defining NVectors xKNext") );
    xKNext = NVector_O::create(iDimensions,_lisp);
    LOG(BF("status") );
    pK = NVector_O::create(iDimensions,_lisp);
    LOG(BF("status") );
    pK->zero();
    LOG(BF("Defining NVectors pjNext") );
    pjNext = NVector_O::create(iDimensions,_lisp);
    LOG(BF("Defining NVectors rj,dj,zj,qj") );
    rj = NVector_O::create(iDimensions,_lisp);
    dj = NVector_O::create(iDimensions,_lisp);
    zj = NVector_O::create(iDimensions,_lisp);
    qj = NVector_O::create(iDimensions,_lisp);
    kSum = NVector_O::create(iDimensions,_lisp);
    mprecon = SparseLargeSquareMatrix_O::create(_lisp,iDimensions,SymmetricUpperDiagonal);
    ldlt=SparseLargeSquareMatrix_O::create(_lisp,iDimensions,SymmetricUpperDiagonal);
    mprecon->fill(0.0);
    ldlt->fill(0.0);
    //
    // Evaluate initial energy and force
    //
    LOG(BF("Evaluating initial energy and force") );
    energyXkNext = dTotalEnergyForce( xK, forceK );
    rmsForceMag = forceK->rmsMagnitude();

    //
    // Setup the preconditioner and carry out UMC
    //
    LOG(BF("Setting up preconditioner") );
    this->_EnergyFunction->setupHessianPreconditioner(xK,mprecon);
    unconventionalModifiedCholeskySymbolicFactorization(mprecon,ldlt);
    unconventionalModifiedCholeskyFactorization(mprecon,ldlt,kSum);

    if ( this->_PrintIntermediateResults ) {
//	_lisp->print(BF("======= Starting Steepest Descent Minimizer\nn" );
//	_lisp->print(BF( "======= Starting Conjugate Gradient Minimizer" ));
	_lisp->print(BF( "======= Starting Truncated Newton Minimizer" ));
    }

    try {
	LOG(BF("Starting loop") );
	while ( 1 ) {


	    if ( this->_DebugOn )
	    {
		stepReport = StepReport_O::create();
		stepReport->_Iteration = this->_Iteration;
	    }


	    //
	    // Inner loop
	    //
	    _truncatedNewtonInnerLoop( kk, xK, mprecon, ldlt,
				       forceK, rmsForceMag, pK,
				       pjNext, rj, dj, zj, qj );


	    //
	    // Line Search

	    prevAlphaK = alphaK;
	    if ( this->_PrintIntermediateResults ) {
		dirMag = pK->magnitude();
		forceMag = forceK->magnitude();
		LOG(BF("Starting descent test") );
		if ( forceMag != 0.0 && dirMag != 0.0 ) {
		    LOG(BF("forceMag = %lf") % forceMag  );
		    LOG(BF("dirMag = %lf") % dirMag  );
		    cosAngle = forceK->dotProduct(pK)/(forceMag*dirMag);
		} else {
		    cosAngle = 0.0;
		}
		if ( this->_PrintIntermediateResults )
		{
		    this->_displayIntermediateMessage(prevAlphaK,energyXkNext,forceMag,rmsForceMag,cosAngle,false);
		}
	    }

	    energyXk = energyXkNext;
	    this->lineSearch( &alphaK, &energyXkNext, xK, pK, forceK,
			      zj, qj, kk, stepReport );
	    xKNext->XPlusYTimesScalar(xK,pK,alphaK);
	    //
	    // Evaluate the force at the new position
	    //
	    fp = dTotalEnergyForce( xKNext, forceK );
	    if ( this->_DebugOn )
	    {
		this->stepReport(stepReport,fp,forceK);
	    }

	    //
	    // Convergence tests
	    //

	    b1aTest=fabs(energyXkNext-energyXk)<EPSILONF*(1.0+fabs(energyXk));
	    LOG(BF("b1aTest = %d") % b1aTest  );
	    LOG(BF("energyXkNext(%lf)") % energyXkNext );
	    LOG(BF("energyXk(%lf)") % energyXk );
	    LOG(BF("fabs[energyXkNext-energyXk]=%le") % fabs(energyXkNext-energyXk) );
	    LOG(BF("[ EPSILONF*(1.0+fabs(energyXk))=%le") % EPSILONF*(1.0+fabs(energyXk)) );
	    if ( b1aTest ) {
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF( "search complete according to b1aTest" ));
		}
		break;
	    }
	    delta = xKNext->rmsDistanceFrom(xK);
	    rmsMagXKNext = xKNext->rmsMagnitude();
	    b1bTest=(delta<SQRT_EPSILONF*(1.0+rmsMagXKNext)/100.0);
	    if ( b1bTest ) {
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF( "search complete according to b1bTest" ));
		}
		break;
	    }

	    rmsForceMag = forceK->rmsMagnitude();
	    if ( rmsForceMag < forceTolerance ) {
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF( "search complete according to absolute force test" ));
		}
		break;
	    }

#if 0 //[
	    b1cTest=rmsForceMag<CUBERT_EPSILONF*(1.0+fabs(energyXkNext));
	    LOG(BF("b1cTest = %d") % b1cTest  );
	    LOG(BF("rmsForceMag [[%le]]<CUBERT_EPSILONF*(1.0+fabs(energyXkNext))[[%le]]") % rmsForceMag % CUBERT_EPSILONF*(1.0+fabs(energyXkNext))  );
	    if ( b1cTest ) {
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF( "search complete according to b1cTest" ));
		}
		break;
	    }
	    b1dTest=rmsForceMag<EPSILONG*(1.0+fabs(energyXkNext));
	    if ( b1dTest ) {
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF( "search complete according to b1dTest" ));
		}
		break;
	    }
#endif //]

	    //
	    // Preparation for next Newton step
	    //
	    // Compute the preconditioner M at X{k+1}
	    //
	    this->_EnergyFunction->setupHessianPreconditioner(xK,mprecon);
	    unconventionalModifiedCholeskyFactorization(mprecon,ldlt,kSum);
	    xK->copyIntoMe(xKNext);
	    kk++;
	    this->_Iteration++;
	    if ( this->_DebugOn )
	    {
		ASSERTNOTNULL(this->_Log);
		this->_Log->addReport(stepReport);
		stepReport = _Nil<StepReport_O>();
	    }
	    if ( kk > numSteps ) {
		MinimizerCondition_ExceededNumSteps fail;
		fail._Message = "truncated newton";
		fail._NumberOfSteps = numSteps;
		throw(fail);
	    }
#if 0
	    if ( this->_StepCallback.notnilp() )
		this->_StepCallback->invoke(this->sharedThis<Minimizer_O>());
#endif
	}
    } catch (MinimizerCondition_ExceededNumSteps fail) {
	if ( this->_DebugOn )
	{
	    if ( stepReport.notnilp() )
	    {
		stepReport->prematureTermination("ExceededNumSteps");
		this->_Log->addReport(stepReport);
		stepReport = _Nil<StepReport_O>();
	    }
	}
	//
	// Lets save the current conformation
	// before throwing this higher
	//
	dTotalEnergyForce( xK, forceK );
	this->_EnergyFunction->writeCoordinatesAndForceToAtoms(xK,forceK);
	MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(fail.message());
    }
    xK->copyIntoMe(xKNext);
    dTotalEnergyForce( xK, forceK );
    this->_EnergyFunction->writeCoordinatesAndForceToAtoms(xK,forceK);
    if ( this->_DebugOn )
    {
	if ( stepReport.notnilp() )
	{
	    stepReport->prematureTermination("Stuck");
	    this->_Log->addReport(stepReport);
	    stepReport = _Nil<StepReport_O>();
	}
    }
    }















/*
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------

*/

    void	Minimizer_O::_evaluateEnergyAndForceManyTimes( int numSteps,  NVector_sp nvPos )
    {_G();
#define MINSLOPE        0.000001
#define MINCHANGE       0.01

    NVector_sp         nvDir, nvNewPos, nvNewForce, nvTempPos, nvTempForce;
    double          dEnergy, dX, dRms;
    int             iCount;
    int		iSize;

    /* Create the Direction vector */
    dRms = 0.0;
    dX = 1.0;
    iSize = nvPos->size();
    nvDir       = NVector_O::create( iSize ,_lisp);
    nvNewPos    = NVector_O::create( iSize ,_lisp);
    nvNewForce = NVector_O::create( iSize ,_lisp);
    nvTempPos   = NVector_O::create( iSize ,_lisp);
    nvTempForce= NVector_O::create( iSize ,_lisp);

    iCount = 0;
    this->_Iteration = 1;
    do {
        dEnergy = this->dTotalEnergyForce( nvPos, nvNewForce);
	if ( iCount % 10000 == 0 ) {
	    _lisp->print(BF("Evaluating energy step#%d") % iCount );
	}
	iCount++;
    } while ( iCount < numSteps );

    LOG(BF("Exceeded maximum number of steps") );

    }


    void	Minimizer_O::validateForce(NVector_sp pos, NVector_sp force)
    {
	ForceMatchReport_sp	report;
	if ( this->_DebugOn ) {
	    report = this->_EnergyFunction->checkIfAnalyticalForceMatchesNumericalForce(pos,force);
	    this->_Log->addReport(report);
	}
    }






/*
  ==========================================================================

  Public routines

*/

/*
 *      Constructor
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Create a new amber amber.
 *      Create the containers for atoms, bonds, angles, and torsions.
 */

    void	Minimizer_O::initialize()
    {
	this->Base::initialize();
#ifdef	USE_CALLBACKS
	this->bFCallback   = NULL;
#endif
	this->_DebugOn = false;
	this->useDefaultSettings();
	this->restart();
//	this->_StepCallback = _Nil<core::LispCallback_O>();
    }



    void	Minimizer_O::useDefaultSettings()
    {
	this->_InitialLineSearchStep = 0.01;
	this->_NumberOfSteepestDescentSteps = MAXSTEEPESTDESCENTSTEPS;
	this->_SteepestDescentTolerance = 2000.0;
	this->_NumberOfConjugateGradientSteps = MAXCONJUGATEGRADIENTSTEPS;
	this->_ConjugateGradientTolerance = 10.0;		//	Use this for now, later add TN minimizer and switch to that when this is <10.0
	this->_NumberOfTruncatedNewtonSteps = MAXTRUNCATEDNEWTONSTEPS;
	this->_TruncatedNewtonTolerance = 0.00000001;
	this->_TruncatedNewtonPreconditioner = hessianPreconditioner;
	this->_PrintIntermediateResults = false;
	LOG(BF("_PrintIntermediateResults = %d") % this->_PrintIntermediateResults  );
	this->_ReportEverySteps = -1;
	this->_Status = minimizerIdle;
	this->_ShowElapsedTime = true;
#ifdef	USE_CALLBACKS
#ifdef	USEBOOSTPYTHON
	this->iPythonCallbackEverySteps = 1;
	this->POPythonCallback = NULL;
#endif
#endif
	this->_MinGradientMean = MIN_GRADIENT_MEAN;
    }





    void	Minimizer_O::enablePrintIntermediateResults()
    {_G();
    this->_PrintIntermediateResults = true;
    }




    void	Minimizer_O::disablePrintIntermediateResults()
    {_G();
    this->_PrintIntermediateResults = false;
    }




    void	Minimizer_O::setEnergyFunction(EnergyFunction_sp f)
    {
	this->_EnergyFunction = f;
	this->_Iteration = 1;
    }

    void	Minimizer_O::evaluateEnergyAndForceManyTimes(int numSteps)
    {_OF();
    NVector_sp	pos;
    if ( this->_EnergyFunction.nilp() )
    {
	SIMPLE_ERROR(BF("You must define an energy function to minimize"));
    }
    this->_Iteration = 1;
    pos = NVector_O::create(this->_EnergyFunction->getNVectorSize(),_lisp);
    this->_EnergyFunction->extractCoordinatesFromAtoms(pos);
    this->_evaluateEnergyAndForceManyTimes(numSteps,pos);
    }



    void	Minimizer_O::minimizeSteepestDescent()
    {_OF();
    NVector_sp	pos;
    int		retries;
    bool		sawProblem;
#if USE_POSIX_TIME
    this->_StartTime = core::PosixTime_O::createNow(_lisp);
#endif
    if ( this->_EnergyFunction.nilp() )
    {
	SIMPLE_ERROR(BF("You must define an energy function to minimize"));
    }
    pos = NVector_O::create(this->_EnergyFunction->getNVectorSize(),_lisp);
    retries = 100;
    try {
	do {
	    sawProblem = false;
	    try {
		this->_EnergyFunction->extractCoordinatesFromAtoms(pos);
		this->_steepestDescent(this->_NumberOfSteepestDescentSteps,pos,
				       this->_SteepestDescentTolerance );
	    } catch ( InteractionProblem ld ) {
		this->_EnergyFunction->dealWithProblem(ld);
		retries--;
		sawProblem = true;
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF("Dealt with problem(%s)") % ld._Message.c_str());
		}
	    }
	} while ( sawProblem && retries > 0 );
    } catch ( MinimizerCondition_ExceededNumSteps fail ) {
	this->_Message.str("");
	this->_Message<<"Steepest descent minimizer exceeded max("
		      << fail._NumberOfSteps<<") number of steps";
	this->_Status = minimizerError;
	MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(this->_Message.str());
    }
    if ( sawProblem ) {
	this->_Status = minimizerError;
	this->_Message << "There was a persistant interaction problem";
	MINIMIZER_ERROR(this->_Message.str());
    }
    }



    void	Minimizer_O::minimizeConjugateGradient()
    {_G();
	NVector_sp	pos;
	int		retries;
	bool		sawProblem;
#if USE_POSIX_TIME
	this->_StartTime = core::PosixTime_O::createNow(_lisp);
#endif
	if ( this->_EnergyFunction.nilp() )
	{
	    SIMPLE_ERROR(BF("You must define an energy function to minimize"));
	}
	pos = NVector_O::create(this->_EnergyFunction->getNVectorSize(),_lisp);
	retries = 100;
	try {
	    do {
		sawProblem = false;
		try {
		    this->_EnergyFunction->extractCoordinatesFromAtoms(pos);
		    this->_conjugateGradient(this->_NumberOfConjugateGradientSteps,pos,
					     this->_ConjugateGradientTolerance );
		} catch ( InteractionProblem ld ) {
		    this->_EnergyFunction->dealWithProblem(ld);
		    retries--;
		    sawProblem = true;
		    if ( this->_PrintIntermediateResults ) {
			_lisp->print(BF("Dealt with problem(%s)") % ld._Message.c_str());
		    }
		}
	    } while ( sawProblem && retries > 0 );
	} catch ( MinimizerCondition_ExceededNumSteps fail ) {
	    this->_Message.str("");
	    this->_Message<<"Conjugate gradient minimizer exceeded max("
			  << fail._NumberOfSteps<<") number of steps";
	    this->_Status = minimizerError;
	    MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(this->_Message.str());
	} catch ( MinimizerCondition_Stuck fail ) {
	    this->_Message.str("");
	    this->_Message<< fail.message() ;
	    this->_Status = minimizerError;
	    MINIMIZER_STUCK_ERROR(this->_Message.str());
	}
	if ( sawProblem ) {
	    this->_Status = minimizerError;
	    this->_Message << "There was a persistant interaction problem";
	    MINIMIZER_ERROR(this->_Message.str());
	}
    }


    void	Minimizer_O::resetAndMinimize()
    {_G();
    this->_Status = minimizerIdle;
    this->minimize();
    }


    void	Minimizer_O::minimize()
    {_G();
    NVector_sp	pos;
    int		retries;
    bool		sawProblem;
#if USE_POSIX_TIME
    this->_StartTime = core::PosixTime_O::createNow(_lisp);
#endif
    if ( this->_EnergyFunction.nilp() )
    {
	SIMPLE_ERROR(BF("You must define an energy function to minimize"));
    }
    pos = NVector_O::create(this->_EnergyFunction->getNVectorSize(),_lisp);
    retries = 100;
    try {
	do {
	    sawProblem = false;
	    try {
		this->_EnergyFunction->extractCoordinatesFromAtoms(pos);
		if ( this->_NumberOfSteepestDescentSteps > 0 ) {
		    this->_steepestDescent( this->_NumberOfSteepestDescentSteps,
					    pos, this->_SteepestDescentTolerance );
		}
		if ( this->_NumberOfConjugateGradientSteps > 0 ) {
		    this->_conjugateGradient( this->_NumberOfConjugateGradientSteps,
					      pos, this->_ConjugateGradientTolerance );
		}
		if ( this->_NumberOfTruncatedNewtonSteps > 0 ) {
		    this->_truncatedNewton( this->_NumberOfTruncatedNewtonSteps,
					    pos, this->_TruncatedNewtonTolerance );
		}
	    } catch ( InteractionProblem ld ) {
		this->_EnergyFunction->dealWithProblem(ld);
		retries--;
		sawProblem = true;
		if ( this->_PrintIntermediateResults ) {
		    _lisp->print(BF("Dealt with problem(%s)") % ld._Message.c_str());
		}
	    }
	} while ( sawProblem && retries > 0 );
    } catch ( MinimizerCondition_ExceededNumSteps fail ) {
	this->_Message.str("");
	this->_Message<<"Minimizer exceeded max("
		      << fail._NumberOfSteps<<") number of steps";
	this->_Status = minimizerError;
	MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(this->_Message.str());
    } catch ( MinimizerCondition_Stuck fail ) {
	this->_Message.str("");
	this->_Message<< fail.message();
	this->_Status = minimizerError;
	MINIMIZER_STUCK_ERROR(this->_Message.str());
    }
    if ( sawProblem ) {
	this->_Status = minimizerError;
	this->_Message << "There was a persistant interaction problem";
	MINIMIZER_ERROR(this->_Message.str());
    }
    }


    geom::QDomNode_sp	Minimizer_O::asXml()
    {
	geom::QDomNode_sp	xml;
	xml = geom::QDomNode_O::create(_lisp,"Minimizer");
	xml->addAttributeDoubleScientific("InitialLineSearchStep",this->_InitialLineSearchStep);
	xml->addAttributeInt("MaximumNumberOfSteepestDescentSteps",this->_NumberOfSteepestDescentSteps);
	xml->addAttributeInt("MaximumNumberOfConjugateGradientSteps",this->_NumberOfConjugateGradientSteps);
	xml->addAttributeInt("MaximumNumberOfTruncatedNewtonSteps",this->_NumberOfTruncatedNewtonSteps);
	xml->addAttributeDoubleScientific("SteepestDescentTolerance",this->_SteepestDescentTolerance);
	xml->addAttributeDoubleScientific("ConjugateGradientTolerance",this->_ConjugateGradientTolerance);
	xml->addAttributeDoubleScientific("TruncatedNewtonTolerance",this->_TruncatedNewtonTolerance);
	xml->addAttributeString("TruncatedNewtonPreconditioner",stringForPreconditionerType(this->_TruncatedNewtonPreconditioner));
	return xml;
    }


    string	Minimizer_O::configurationAsString()
    {
	stringstream	ss;
	ss.str("");
	ss << "Minimizer state:"<<std::endl;
	ss << "InitialLineSearchStep:                 " <<this->_InitialLineSearchStep << std::endl;
	ss << "MaximumNumberOfSteepestDescentSteps:   "<<this->_NumberOfSteepestDescentSteps << std::endl;
	ss << "SteepestDescentTolerance:              "<<this->_SteepestDescentTolerance << std::endl;
	ss << "MaximumNumberOfConjugateGradientSteps: "<<this->_NumberOfConjugateGradientSteps << std::endl;
	ss << "ConjugateGradientTolerance:            "<<this->_ConjugateGradientTolerance << std::endl;
	ss << "MaximumNumberOfTruncatedNewtonSteps: "<<this->_NumberOfTruncatedNewtonSteps << std::endl;
	ss << "TruncatedNewtonTolerance:            "<<this->_TruncatedNewtonTolerance << std::endl;
	ss << "TruncatedNewtonPreconditioner:       "<<stringForPreconditionerType(this->_TruncatedNewtonPreconditioner)<<std::endl;
	return ss.str();
    }

    void	Minimizer_O::restart()
    {
	this->_Status = minimizerIdle;
	this->_Message.str("");
	this->_Iteration = 1;
    }

    void Minimizer_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();

	SYMBOL_EXPORT_SC_(ChemPkg,noPreconditioner);
	SYMBOL_EXPORT_SC_(ChemPkg,hessianPreconditioner);
	core::enum_<PreconditionerType>(_sym__PLUS_preconditionerTypeConverter_PLUS_,"PreconditionerType")
	    .value(_sym_noPreconditioner,noPreconditioner)
	    .value(_sym_hessianPreconditioner,hessianPreconditioner)
	    ;
	SYMBOL_EXPORT_SC_(ChemPkg,minimizerError);
	SYMBOL_EXPORT_SC_(ChemPkg,minimizerSucceeded);
	SYMBOL_EXPORT_SC_(ChemPkg,truncatedNewtonRunning);
	SYMBOL_EXPORT_SC_(ChemPkg,conjugateGradientRunning);
	SYMBOL_EXPORT_SC_(ChemPkg,steepestDescentRunning);
	SYMBOL_EXPORT_SC_(ChemPkg,minimizerIdle);
	core::enum_<MinimizerStatus>(_sym__PLUS_minimizerStatusConverter_PLUS_,"MinimizerStatus")
	    .value(_sym_minimizerIdle,minimizerIdle)
	    .value(_sym_steepestDescentRunning,steepestDescentRunning)
	    .value(_sym_conjugateGradientRunning,conjugateGradientRunning)
	    .value(_sym_truncatedNewtonRunning,truncatedNewtonRunning)
	    .value(_sym_minimizerSucceeded,minimizerSucceeded)
	    .value(_sym_minimizerError,minimizerError)
	    ;
    }


    EXPOSE_CLASS_AND_GLOBALS(chem,Minimizer_O);

};

