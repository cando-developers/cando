/*
    File: minimizer.cc
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
// #undef USEBOOSTPYTHON

#define	DEBUG_LEVEL_NONE

//#define	TURN_LINESEARCH_DETAILS_ON
#define	TURN_LINESEARCH_DETAILS_OFF

// #define DEBUG_FACTORIZATION

#define	VALIDATE_FORCE_OFF	// Turn this on to compare analytical force to numerical force at every step

#define	USE_NEW_CONJUGATE_GRADIENTS
#include <clasp/core/common.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/minimizer.h>
#include <iostream>
#include <cando/geom/vector3.fwd.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/chem/ffBaseDb.h>
#include <clasp/core/array.h>
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
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyPointToLineRestraint.h>
#include <cando/chem/energyOutOfZPlane.h>
#include <cando/chem/energyDihedralRestraint.h>
#include <cando/chem/energyFixedNonbond.h>
#include <clasp/core/wrappers.h>


SYMBOL_EXPORT_SC_(ChemKwPkg,message);

SYMBOL_SC_(ChemPkg,_PLUS_minimizerStatusConverter_PLUS_);
SYMBOL_SC_(ChemPkg,_PLUS_preconditionerTypeConverter_PLUS_);

SYMBOL_EXPORT_SC_(KeywordPkg,number_of_steps);
SYMBOL_EXPORT_SC_(KeywordPkg,minimizer);
SYMBOL_EXPORT_SC_(KeywordPkg,localSteps);

SYMBOL_EXPORT_SC_(ChemPkg,steepest_descent);
SYMBOL_EXPORT_SC_(ChemPkg,conjugate_gradient);
SYMBOL_EXPORT_SC_(ChemPkg,truncated_newton);
SYMBOL_EXPORT_SC_(ChemPkg,truncated_newton_debug);

SYMBOL_EXPORT_SC_(ChemPkg,MinimizerExceededSD_MaxSteps);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerExceededCG_MaxSteps);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerExceededTN_MaxSteps);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerExceededMaxSteps);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerStuck);
SYMBOL_EXPORT_SC_(ChemPkg,MinimizerError);

SYMBOL_EXPORT_SC_(ChemPkg,STARsave_positionsSTAR);
namespace chem
{


#define MINIMIZER_ERROR(msg) ERROR(_sym_MinimizerError,core::lisp_createList(kw::_sym_message,core::Str_O::create(msg)))
#define MINIMIZER_EXCEEDED_MAX_STEPS_ERROR(msg) ERROR(_sym_MinimizerExceededMaxSteps,core::lisp_createList(kw::_sym_minimizer,msg._Minimizer, kw::_sym_number_of_steps, core::make_fixnum(msg._NumberOfSteps)));
#define MINIMIZER_STUCK_ERROR(msg,localSteps) ERROR(_sym_MinimizerStuck,core::lisp_createList(kw::_sym_message,core::SimpleBaseString_O::make(msg), kw::_sym_localSteps, core::clasp_make_fixnum(localSteps)))



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
{
  IMPLEMENT_ME();
#if 0
  core::HashTable_sp options = core::HashTable_O::createFromKeywordCons(optionCons,minimizerOptions);
  options->setBoolValueIfAvailable(this->_ShowElapsedTime,":showElapsedTime");
#endif
}

CL_LAMBDA(energy-function);
CL_LISPIFY_NAME(make-minimizer);
DOCGROUP(cando);
CL_DEFUN Minimizer_sp Minimizer_O::make(ScoringFunction_sp givenEnergyFunction)
{
  auto  me  = gctools::GC<Minimizer_O>::allocate_with_default_constructor();
  bool initialized = false;
  initialized = true;
  me->setEnergyFunction(givenEnergyFunction);
  return me;
}




CL_LISPIFY_NAME("statusAsString");
CL_DEFMETHOD     string	Minimizer_O::statusAsString()
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

void Minimizer_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,InitialLineSearchStep),this->_InitialLineSearchStep);
  node->field(INTERN_(kw,ShowElapsedTime),this->_ShowElapsedTime);
  node->field(INTERN_(kw,SteepestDescentTolerance),this->_SteepestDescentTolerance);
  node->field(INTERN_(kw,NumberOfSteepestDescentSteps),this->_NumberOfSteepestDescentSteps);
  node->field(INTERN_(kw,NumberOfConjugateGradientSteps),this->_NumberOfConjugateGradientSteps);
  node->field(INTERN_(kw,ConjugateGradientTolerance),this->_ConjugateGradientTolerance);
  node->field(INTERN_(kw,NumberOfTruncatedNewtonSteps),this->_NumberOfTruncatedNewtonSteps);
  node->field(INTERN_(kw,TruncatedNewtonTolerance),this->_TruncatedNewtonTolerance);
  node->field(INTERN_(kw,ScoringFunction),this->_ScoringFunction );
  node->field(INTERN_(kw,PrintIntermediateResults),this->_PrintIntermediateResults);
}

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

CL_LISPIFY_NAME(minimizer-set-initial-line-search-step);
CL_DEFMETHOD
void Minimizer_O::set_initial_line_search_step(double step)
{
  this->_InitialLineSearchStep = step;
}

/*
 *	getPosition
 *
 *	Calculate the position from the nvOrigin+nvDirection*x
 */
void	Minimizer_O::getPosition( NVector_sp 	nvResult,
                                  NVector_sp	nvOrigin,
                                  NVector_sp	nvDirection,
                                  double		x,
                                  core::T_sp activeAtomMask )
{
  XPlusYTimesScalarWithActiveAtomMask(nvResult,nvOrigin,nvDirection,x,activeAtomMask);
}


void maybeSavePosition(NVector_sp pos) {
  if (chem::_sym_STARsave_positionsSTAR->symbolValue().notnilp()) {
    core::T_sp ovec = chem::_sym_STARsave_positionsSTAR->symbolValue();
    size_t iDimensions = pos->size();
    auto tempPos = core::SimpleVector_float_O::create(iDimensions);
    if (gc::IsA<core::ComplexVector_sp>(ovec)) {
      auto vec = gc::As_unsafe<core::ComplexVector_sp>(ovec);
      for (size_t idx=0; idx<iDimensions; idx++ ) {
        (*tempPos)[idx] = (*pos)[idx];
      }
      vec->vectorPushExtend(tempPos);
      return;
    }
    SIMPLE_ERROR("Add support for pos with the value {}", _rep_(pos));
  }
}


/*
 *      dTotalEnergy
 *
 *      Evaluate the total amber.
 *      Evaluate the amber at nvPos and return the derivative
 *      in nvForce (which must be an initialized vector,
 *      but does not need to be filled with zeros).
 */
double	Minimizer_O::dTotalEnergy( NVector_sp pos,
                                   core::T_sp energyScale,
                                   core::T_sp activeAtomMask )
{
  double          dEnergy;

  dEnergy = this->_ScoringFunction->evaluateEnergy( pos,
                                                    energyScale,
                                                    nil<core::T_O>(),
                                                    activeAtomMask );
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
double	Minimizer_O::dTotalEnergyForce( NVector_sp pos,
                                        core::T_sp energyScale,
                                        NVector_sp nvForce, core::T_sp activeAtomMask )
{
  return this->_ScoringFunction->evaluateEnergyForce(pos,
                                                     energyScale,
                                                     true,nvForce,activeAtomMask);
}


/*
 *	d1DTotalEnergy
 *
 *	Calculate the energy along a 1D coordinate
 */
double	Minimizer_O::d1DTotalEnergy( double x, core::T_sp energyScale, core::T_sp activeAtomMask )
{
#ifdef	DEBUG_ON
//    this->nvP1DSearchOrigin->debugDump("origin");
//    LOG("x = {}" , (x ) );
//    this->nvP1DSearchDirection->debugDump("direction");
#endif
  this->getPosition(this->nvP1DSearchTemp1, this->nvP1DSearchOrigin, this->nvP1DSearchDirection,x, activeAtomMask );
  return this->_ScoringFunction->evaluateEnergy(this->nvP1DSearchTemp1,energyScale,nil<core::T_O>(),activeAtomMask);
}


/*
 *	d1DTotalEnergyForce
 *
 *	Calculate the energy/derivative along a 1D coordinate
 */
double	Minimizer_O::d1DTotalEnergyForce( double x, core::T_sp energyScale, double* fx, double * dfx, core::T_sp activeAtomMask )
{
#ifdef	DEBUG_ON
//    this->nvP1DSearchOrigin->debugDump("origin");
//    LOG("x = {}" , (x ) );
//    this->nvP1DSearchDirection->debugDump("direction");
#endif
  this->getPosition(this->nvP1DSearchTemp1,
                    this->nvP1DSearchOrigin,
                    this->nvP1DSearchDirection, x, activeAtomMask );
  *fx = this->_ScoringFunction->evaluateEnergyForce( this->nvP1DSearchTemp1,
                                                     energyScale,
                                                     true, this->nvP1DSearchTemp2,
                                                     activeAtomMask );
  *dfx = -dotProductWithActiveAtomMask(this->nvP1DSearchTemp2,this->nvP1DSearchDirection,activeAtomMask);
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
                                core::T_sp      energyScale,
                                NVector_sp	nvDir,
                                double		*dPxa,
                                double		*dPxb,
                                double		*dPxc,
                                double		*dPfa,
                                double		*dPfb,
                                double		*dPfc,
                                core::T_sp activeAtomMask )
{
  double		xa,xb,xc,fa,fb,fc, temp;
  double		r,q,u,ulim,fu;
  size_t numSteps = 0;
  this->_MinBracketSteps = 0;
  xa = *dPxa;
  xb = *dPxb;
  fa = this->d1DTotalEnergy(xa,energyScale,activeAtomMask);
  fb = this->d1DTotalEnergy(xb,energyScale,activeAtomMask);
    // Make sure that we are going downhill a->b
  if ( fb > fa ) {
    SWAP_VALUES( xa, xb, temp );
    SWAP_VALUES( fa, fb, temp );
  }
  xc = xb+GOLD*(xb-xa);
  fc = this->d1DTotalEnergy(xc,energyScale,activeAtomMask);
  LOG("Start: xa({}) xb({}) xc({}) | fa({}) fb({}) fc({})" , xa , xb , xc , fa , fb , fc );
  while ( fb > fc ) {
    numSteps++;
    if (numSteps%1000==0) {
      printf("%s:%d The minBracket function is failing with %lu steps xa=%f xb=%f  xc=%f fa=%f fb=%f fc=%f\n",
             __FILE__, __LINE__, numSteps, xa, xb, xc, fa, fb, fc );
    }
    this->_MinBracketSteps++;
    LOG("Loop:  xa({}) xb({}) xc({}) | fa({}) fb({}) fc({})"%xa%xb%xc%fa%fb%fc);
    r = (xb-xa)*(fb-fc);
    q = (xb-xc)*(fb-fa);
    u = xb-((xb-xc)*q-(xb-xa)*r)/
      (2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
    ulim = (xb)+GLIMIT*(xc-xb);
    if (( xb-u)*(u-xc)>0.0) {
      fu = this->d1DTotalEnergy(u,energyScale,activeAtomMask);
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
      fu = this->d1DTotalEnergy(u,energyScale,activeAtomMask);
    } else if ((xc-u)*(u-ulim) > 0.0) {
      fu = this->d1DTotalEnergy(u,energyScale,activeAtomMask);
      if ( fu < fc ) {
        LEFT_SHIFT_VALUES( xb, xc, u, xc+GOLD*(xc-xb) );
        LEFT_SHIFT_VALUES( fb, fc, fu, this->d1DTotalEnergy(u,energyScale,activeAtomMask));
      }
	    // Limit parabolic u to maximum allowed value
    } else if ((u-ulim)*(ulim-xc)>=0.0 ){
      u = ulim;
      fu = this->d1DTotalEnergy(u,energyScale,activeAtomMask);
    } else {
      u = xc+GOLD*(xc-xb);
      fu = this->d1DTotalEnergy(u,energyScale,activeAtomMask);
    }
    LEFT_SHIFT_VALUES( xa, xb, xc, u );
    LEFT_SHIFT_VALUES( fa, fb, fc, fu );
  }
 DONE:
#ifdef	DEBUG_ON
  LOG("minBracket  xa({}), xb({}), xc({})" , xa , xb , xc );
  if ( fa > fb && fb < fc ) {
    LOG("fa.GT.fb.LT.fc" );
  } else {
    LOG("FAIL! minBracket FAILED!! It's not true that: fa.GT.fb.LT.fc" );
  }
  if ( xc < 0 ) {
    LOG("ATTN!  xc is less than zero " );
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


//__attribute__((optnone))
extern "C" {
__attribute__((optnone))
static void safe_division(double& result, double x, double y)
{
  if (y!=0) result = x/y;
}
}

/*
 *	dbrent
 *
 *	Minimize the function within the interval ax,bx,cx
 *	Return the function value at the new point and the
 *	new minimum in (xmin)
 */
double 	Minimizer_O::dbrent(	double ax, double bx, double cx,
                                core::T_sp energyScale,
                                double tol,
                                double& lineStep,
                                int&	energyEvals,
                                int&	forceEvals,
                                int&	dbrentSteps,
                                core::T_sp activeAtomMask
                                )
{
  int	iter, ok1, ok2;
  double	_a,_b,_d,d1,d2,du,dv,dw,dx,_e=0.0;
  double	fu,fv,fw,fx,olde,tol1,tol2, u, u1, u2, v,w,x,xm, ft;
  double	retval;

  _d = 0.0;
  energyEvals = 0;
  forceEvals = 0;
  dbrentSteps = 0;
  LOG("Incoming ax bx cx= {} {} {}" , ax , bx , cx );
  if (chem__verbose(4) ) {
    core::clasp_write_string(fmt::format("Incoming ax bx cx -> {} {} {}\n" , ax , bx , cx));
  }
  _a = (ax<cx?ax:cx);
  _b = (ax>cx?ax:cx);
  LOG("Initial _a _b= {} {}" , _a , _b );
  x = w = v = bx;
//
// Calculate the derivative along the search direction
//
  ft = d1DTotalEnergyForce( x, energyScale, &fx, &dx, activeAtomMask );
  forceEvals++;
  LOG("dbrent: derivative x,fx,dx = {} {} {}" , x , fx , dx  );
  fw=fv=fx;
  dw=dv=dx;
  for (iter=1;iter<=DBRENTITMAX;iter++ ) {
    dbrentSteps++;
//	LOG("dbrent: iter=%3d  _a, _b= %10.15lf %10.15lf" , iter , _a , _b  );
    if (chem__verbose(4) ) {
      core::clasp_write_string(fmt::format("dbrent: iter = {:3}  _a, _b -> {} {}\n" , iter , _a , _b ));
    }
    xm=0.5*(_a+_b);
    tol1=tol*fabs(x)+ZEPS;
    tol2=2.0*tol1;
    if (fabs(x-xm) <=(tol2-0.5*(_b-_a))) {	// Stopping criterion
      lineStep=x;
      retval = fx;
      LOG("done due to (fabs(x-xm)<=(tol2-0.5*(_b-_a)))" );
      if (chem__verbose(4) ) {
        core::clasp_write_string(fmt::format("dbrent: done due to (fabs(x-xm) <= (tol2-0.5*(_b-_a))) -> {} <= {}\n" , fabs(x-xm) , (tol2-0.5*(_b-_a))));
      }
      goto DONE;
    }
    if (fabs(_e)>tol1){
      d1=2.0*(_b-_a);
      d2=d1;

      double numerator;
      double denomenator;
      double limit;
      //if (dw != dx ) d1=(w-x)*dx/(dx-dw); // Original code Secant method, first on one, then on
#if 0
# define WARN_BAD_VAL(xx) if (std::isnan(xx) || std::isinf(xx)) printf("%s:%d:%s problem with " #xx " - it is %e\n", __FILE__, __LINE__, __FUNCTION__, xx );
#else
# define WARN_BAD_VAL(xx)
#endif
      WARN_BAD_VAL(dw);
      WARN_BAD_VAL(dx);
      WARN_BAD_VAL(dv);
      // if ((dw-dx) != 0.0) d1=(w-x)*dx/(dx-dw); // Secant method, first on one, then on
      safe_division( d1, (w-x)*dx, dx-dw );

      // if (dv!=dx) d2=(v-x)*dx/(dx-dv); // the other point
      // if ((dv-dx) != 0.0) d2=(v-x)*dx/(dx-dv); // the other point
      safe_division( d2, (v-x)*dx, dx-dv );

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
      fu=d1DTotalEnergy(u,energyScale,activeAtomMask);
      energyEvals++;
    } else {
      u = x+SIGN(tol1,_d);
      fu=d1DTotalEnergy(u,energyScale,activeAtomMask);
      energyEvals++;
      if ( fu>fx) { // If the minimum lineStep in the downhill direction takes us
			  // uphill, then we are done
        lineStep=x;
        retval = fx;
        LOG("done due to fu>fx" );
        if (chem__verbose(4) ) {
          core::clasp_write_string(fmt::format("dbrent: done due to fu>fx -> {} > {}\n" , fu , fx ));
        }
        goto DONE;
      }
    }
//
// Calculate the force along the search direction
//
    ft = d1DTotalEnergyForce( u, energyScale, &ft, &du, activeAtomMask ); // Now housekeeping
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
  LOG("dbrent: ERROR Exceeded max iterations" );
  if (chem__verbose(4) ) {
    core::clasp_write_string(fmt::format("dbrent: exceeded max number of iterations {}\n" , DBRENTITMAX));
  }
  retval = -1.0;
 DONE:
  LOG("dbrent evaluated energy({}) and force({}) times" , energyEvals , forceEvals  );
  return retval;
}




void Minimizer_O::lineSearchInitialReport( StepReport_sp report,
                                           NVector_sp nvPos,
                                           core::T_sp energyScale,
                                           NVector_sp nvDir, NVector_sp nvForce,
                                           double xa, double xb, double xc,
                                           double fa, double fb, double fc,
                                           core::T_sp activeAtomMask )
{
  double lenForce, lenDir, angle, cosAngle;
  lenForce = magnitudeWithActiveAtomMask(nvForce,activeAtomMask);
  lenDir = magnitudeWithActiveAtomMask(nvDir,activeAtomMask);
  if ( lenForce == 0.0 || lenDir == 0.0 ) {
    angle = 200.0;
  } else {
    cosAngle = dotProductWithActiveAtomMask(nvDir,nvForce,activeAtomMask)/(lenForce*lenDir);
    if ( cosAngle > 1.0 ) cosAngle = 1.0;
    if ( cosAngle < -1.0 ) cosAngle = -1.0;
    angle = ::safe_acos(cosAngle);
  }
  report->_AngleBetweenDirectionAndForceDeg = angle/0.0174533;
  report->_Xa = xa;
  report->_Xb = xb;
  report->_Xc = xc;
  report->_Fa = fa;
  report->_Fb = fb;
  report->_Fc = fc;
  report->_MinBracketSteps = this->_MinBracketSteps;
  report->_EnergyTermsEnabled = this->_ScoringFunction->energyTermsEnabled();
  report->_TotalEnergy = this->d1DTotalEnergy(0.0,energyScale,activeAtomMask);
  report->_DirectionMagnitude = magnitudeWithActiveAtomMask(nvDir,activeAtomMask);
  report->_ForceMagnitude = magnitudeWithActiveAtomMask(nvForce,activeAtomMask);
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
  report->_LineSearchPosition = copy_nvector(nvPos);
  report->_LineSearchDirection = copy_nvector(nvDir);
  report->_TotalEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("Total"),xmin,xinc);
  report->_StretchEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("Stretch"),xmin,xinc);
  report->_AngleEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("Angle"),xmin,xinc);
  report->_DihedralEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("Dihedral"),xmin,xinc);
  report->_NonbondEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("Nonbond"),xmin,xinc);
  report->_ImproperEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("Improper"),xmin,xinc);
  report->_ChiralRestraintEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("ChiralRestraint"),xmin,xinc);
  report->_AnchorRestraintEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("AnchorRestraint"),xmin,xinc);
  report->_PointToLineRestraintEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("PointToLineRestraint"),xmin,xinc);
  report->_OutOfZPlaneEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("OutOfZPlane"),xmin,xinc);
  report->_FixedNonbondRestraintEnergyFn = NumericalFunction_O::create(core::SimpleBaseString_O::make("Alpha"),core::SimpleBaseString_O::make("FixedNonbondRestraint"),xmin,xinc);

  for ( zx=dxa;zx<=dxc;zx+=(dxc-dxa)/100.0 ) {
    zy = this->d1DTotalEnergy(zx,energyScale,activeAtomMask);
    report->_TotalEnergyFn->appendValue(zy);
#if 0
            // skipping components - it's not general
    report->_StretchEnergyFn->appendValue(
                                          this->_EnergyFunction->getStretchComponent()->getEnergy());
    report->_AngleEnergyFn->appendValue(this->_EnergyFunction->getAngleComponent()->getEnergy());
    report->_DihedralEnergyFn->appendValue(this->_EnergyFunction->getDihedralComponent()->getEnergy());
    report->_NonbondEnergyFn->appendValue(this->_EnergyFunction->getNonbondComponent()->getEnergy());
    report->_ImproperEnergyFn->appendValue(this->_EnergyFunction->getImproperRestraintComponent()->getEnergy());
    report->_ChiralRestraintEnergyFn->appendValue(this->_EnergyFunction->getChiralRestraintComponent()->getEnergy());
    report->_AnchorRestraintEnergyFn->appendValue(this->_EnergyFunction->getAnchorRestraintComponent()->getEnergy());
    report->_PointToLineRestraintEnergyFn->appendValue(this->_EnergyFunction->getPointToLineRestraintComponent()->getEnergy());
    report->_OutOfZPlaneEnergyFn->appendValue(this->_EnergyFunction->getOutOfZPlaneComponent()->getEnergy());
    report->_FixedNonbondRestraintEnergyFn->appendValue(this->_EnergyFunction->getFixedNonbondRestraintComponent()->getEnergy());
#endif
  }
};


CL_LISPIFY_NAME("throwMinimizerError");
CL_DEFMETHOD     void Minimizer_O::throwMinimizerError()
{
  MINIMIZER_ERROR("test throw of MinimizerError");
};



void Minimizer_O::lineSearchFinalReport( StepReport_sp report, double step, double fMin,
                                         int energyEvals, int forceEvals, int dbrentSteps )
{
  ASSERT(report->_Iteration == this->_Iteration);
  report->_DbrentSteps = dbrentSteps;
  report->_EnergyEvals = energyEvals;
  report->_ForceEvals = forceEvals;
  report->_Step = step;
  report->_FMin = fMin;
}


void	Minimizer_O::stepReport( StepReport_sp report, double energy, NVector_sp force, core::T_sp activeAtomMask )
{
  ASSERT(report->_Iteration == this->_Iteration);
  report->_ForceMagnitude = magnitudeWithActiveAtomMask(force,activeAtomMask);
  report->_TotalEnergy = energy;
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
                                        core::T_sp energyScale,
                                        NVector_sp nvDirection,
                                        NVector_sp nvForce,
                                        NVector_sp nvTemp1,
                                        NVector_sp nvTemp2,
                                        int	iteration,
                                        StepReport_sp	report,
                                        core::T_sp activeAtomMask )
{
  double	xa, xb, xc;
  double	fa, fb, fc;
  double	step = 0.0;
  int	functionEvals;

  functionEvals = 0;
  LOG("Starting" );

    //
    // Setup the 1D search origin and direction
    //

  this->define1DSearch(nvOrigin,nvDirection,nvTemp1,nvTemp2);
  double directionMag = magnitudeWithActiveAtomMask(nvDirection,activeAtomMask);
  if ( directionMag < VERYSMALLSQUARED ) {
    xb = this->_InitialLineSearchStep;
  } else {
    xb = this->_InitialLineSearchStep/directionMag;
  }
  xa = 0.0;

    //
    // Bracket the minimum
    //
  if ( chem__verbose(4) ) {
    core::clasp_write_string(fmt::format("Trying to bracket min this->_InitialLineSearchStep,directionMag,xb -> {},{},{}\n" , this->_InitialLineSearchStep , directionMag , xb));
  }

  this->minBracket( nvOrigin, energyScale, nvDirection,
                    &xa, &xb, &xc, &fa, &fb, &fc, activeAtomMask );
  LOG("Bracketed minimum" );
  LOG("xa,xb,xc = {} {} {} " , xa , xb , xc  );
  LOG("fa,fb,fc = {} {} {} " , fa , fb , fc  );
  if ( chem__verbose(4) ) {
    core::clasp_write_string(fmt::format("Bracketed min xa,xb,xc -> {},{},{}\n" , xa , xb , xc ));
    core::clasp_write_string(fmt::format("              fa,fb,fc -> {},{},{}\n" , fa , fb , fc ));
  }

  if ( this->_DebugOn )
  {
    this->lineSearchInitialReport(report,
                                  nvOrigin,
                                  energyScale,
                                  nvDirection,
                                  nvForce,
                                  xa,xb,xc,fa,fb,fc,activeAtomMask);
  }

    //
    // Now use the bracketed minimum to find the line minimum
    //
    //
    // THIS IS WHERE YOU ADD DBRENT
  int energyEvals = 0;
  int forceEvals = 0;
  int dbrentSteps = 0;
  fb = dbrent( xa, xb, xc, energyScale, TOL, step, energyEvals, forceEvals, dbrentSteps, activeAtomMask );
  LOG("brent bracketed step = {}" , step  );
  if ( this->_DebugOn )
  {
    this->lineSearchFinalReport( report, step, fb, energyEvals,
                                 forceEvals, dbrentSteps );
  }

  *dPstep = step;
  *dPfnew = fb;


  LOG("Done" );
}



bool	Minimizer_O::_displayIntermediateMessage(NVector_sp       pos,
                                                 double			step,
                                                 double			fenergy,
                                                 double			forceRMSMag,
                                                 double			cosAngle,
                                                 bool			steepestDescent,
                                                 bool                   forcePrint,
                                                 core::T_sp             activeAtomMask )
{
#define	MAX_DISPLAY	1024
  double		angle;
  char		buffer[MAX_DISPLAY];
  stringstream	sout;
  if ( forcePrint || this->_Iteration % this->_ReportEveryNSteps == 0 ) {
    if ( this->_Iteration%(10*this->_ReportEveryNSteps) == 1 || this->_DebugOn )
    {
      if (this->_DebugOn) {
        sout << "DBG";
      }
      sout << "---Stage-";
      if ( this->_ShowElapsedTime )
      {
        sout << "Seconds--";
      }
      sout << "Step-log(Alpha)--Prev.dir--------Energy-----------RMSforce";
      if ( this->_ScoringFunction->scoringFunctionName().notnilp() ) 
      {
        sout << "-------Name";
      }
      sout << std::endl;
    }
    if (this->_DebugOn) {
      sout << "   ";
    }
    sout << fmt::format(" min{:4}" , this->statusAsShortString());
    if ( this->_ShowElapsedTime )
    {
      auto EndTime = std::chrono::steady_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(EndTime - this->_StartTime._value);
      sout << fmt::format(" {:7}" , elapsed.count());
    }
    sout << fmt::format(" {:5}" , this->_Iteration);
    if (step==0.0) {
      sout << fmt::format(" {:9s}" , "nan" );
    } else {
      sout << fmt::format(" {:9.4f}" , log(step));
    }
    if ( steepestDescent ) 
    {
      sout << "StDesc";
    } else 
    {
      if ( cosAngle < -1.0 ) cosAngle = -1.0;
      if ( cosAngle > 1.0 ) cosAngle = 1.0;
      angle = ::safe_acos(cosAngle)/0.0174533;
      if ( angle < 0.1 ) 
      {
        angle = 0.0;
      }
      if (angle < 0.001) {
        core::clasp_writeln_string(fmt::format("  prev.dir has gone to zero\n"));
      }
      sout << (fmt::format(" {:5.1f}" , angle));
    }
    sout << (fmt::format(" {:18.5f}" , fenergy ));
    sout << (fmt::format(" {:18.5f}" , forceRMSMag));
    if ( this->_ScoringFunction->scoringFunctionName().notnilp() ) 
    {
      sout << fmt::format(" {}" , this->_ScoringFunction->scoringFunctionName());
    }
    core::clasp_writeln_string(sout.str());
    if (this->_Iteration%10==0) {
      std::string sstr = gc::As<core::String_sp>(core::eval::funcall(INTERN_(chem,summarize_energy_function),this->_ScoringFunction,pos,activeAtomMask))->get_std_string();
      core::clasp_writeln_string(sstr);
      core::cl__finish_output();
    }
    if ( this->_DebugOn ) {
      this->_Log->addMessage(buffer);
    }
    return true;
  }
  return false;
}




/*!
 * Steepest descent minimizer with or without preconditioning
 */
void	Minimizer_O::_steepestDescent( int numSteps,
                                       NVector_sp pos,
                                       core::T_sp energyScale,
                                       double forceTolerance,
                                       core::T_sp activeAtomMask,
                                       core::T_sp callback )
{
  StepReport_sp	stepReport = StepReport_O::create();
  int		iRestartSteps;
  NVector_sp	forceVec, m;
  NVector_sp	precon_dir,dirVec,dirPrev,tv1,tv2;
  double	forceMag, forceRmsMag,prevStep;
  double	delta0, deltaNew;
  double	eSquaredDelta0;
  double	step, fnew, dirMag;
  core::DoubleFloat_sp dstep;
  double	cosAngle = 0.0;
  bool		steepestDescent;
  int		innerSteps;
  int		localSteps, k;
  bool          printedLatestMessage;

  if ( this->_PrintIntermediateResults ) {
    core::clasp_write_string(fmt::format("++++ Steepest descent ++++ max_steps = {}     tolerance = {}\n", numSteps, forceTolerance ));
  }
  dstep = core::DoubleFloat_O::create(0.0);
  LOG("Checking status" );
  if ( this->_Status == minimizerError ) return;
  this->_Status = steepestDescentRunning;
  this->_CurrentPreconditioner = noPreconditioner;
  LOG("step" );

  localSteps = 0;
  k = 0;
  step = 0.0;

  LOG("step" );
	/* Calculate how many conjugate gradient steps can be */
	/* taken before a restart must be done */

  iRestartSteps = pos->size();
	// Define NVectors
  forceVec = NVector_O::create(iRestartSteps);
  this->_Force = forceVec;
  precon_dir = NVector_O::create(iRestartSteps);
  LOG("step" );
  dirVec = NVector_O::create(iRestartSteps);
  dirPrev = NVector_O::create(iRestartSteps);
  tv1 = NVector_O::create(iRestartSteps);
  tv2 = NVector_O::create(iRestartSteps);
  LOG("step" );
	// Done
  innerSteps = MIN(iRestartSteps,ITMAX);
  LOG("step" );
  double fp = this->dTotalEnergyForce( pos, energyScale, forceVec, activeAtomMask );
  LOG("step" );
//    r->inPlaceTimesScalar(-1.0);
	//  no preconditioning
  copyVector(dirVec,forceVec);
  deltaNew = dotProductWithActiveAtomMask(forceVec,dirVec,activeAtomMask);
  delta0 = deltaNew;
  eSquaredDelta0 = forceTolerance*delta0;
  LOG("eSquaredDelta0 = {}" , (eSquaredDelta0 ) );
  LOG("forceTolerance = {}" , (forceTolerance ) );
  LOG("delta0 = {}" , (delta0 ) );
  localSteps = 0;
  {
    while (1) { 
		//
		// Absolute gradient test
		//
      forceMag = magnitudeWithActiveAtomMask(forceVec,activeAtomMask);
      forceRmsMag = rmsMagnitudeWithActiveAtomMask(forceVec,activeAtomMask);
      this->_RMSForce = forceRmsMag;
      		    //
		    // Print intermediate status
		    //
      if ( forceRmsMag < forceTolerance ) {
        if ( this->_PrintIntermediateResults ) {
          if (!printedLatestMessage) {
            printedLatestMessage = this->_displayIntermediateMessage(pos,step,fp,forceRmsMag,cosAngle,steepestDescent,true,activeAtomMask);
          }
          core::clasp_writeln_string(fmt::format(" ! DONE absolute force test:\n ! forceRmsMag({}) < forceTolerance({})" , forceRmsMag , forceTolerance) );
        }
        break;
      }
      if ( localSteps>=numSteps ) {
        if ( this->_DebugOn ) {
          ANN(stepReport);
          if ( stepReport.notnilp() ) {
            stepReport->prematureTermination("ExceededNumSteps");
            this->_Log->addReport(stepReport);
          }
        }
	    //
	    // Lets save the current conformation
	    // before throwing this higher
	    //
        fp = dTotalEnergyForce( pos, energyScale, forceVec, activeAtomMask );
        this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(pos,forceVec);
        ERROR(_sym_MinimizerExceededSD_MaxSteps, (ql::list() 
                                               << kw::_sym_minimizer << this->asSmartPtr()
                                               << kw::_sym_number_of_steps << core::make_fixnum(localSteps)
                                               << kw::_sym_coordinates << pos).result());
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
        if ( this->_DebugOn ) {
          stepReport = StepReport_O::create();
          stepReport->_Iteration = this->_Iteration;
        }
		    // Descent test
		    // If the cos(angle) between the search direction and
		    // the force (steepest descent dir)
		    // is less than _xxxDescentTest then copy the force into the search dir
        dirMag = magnitudeWithActiveAtomMask(dirVec,activeAtomMask);
        steepestDescent = false;
        LOG("Starting descent test" );
        if ( forceMag != 0.0 && dirMag != 0.0 ) {
          // Calculate angle between dirVec and dirPrev
          //cosAngle = dotProductWithActiveAtomMask(forceVec,dirVec,activeAtomMask)/(forceMag*dirMag);
          cosAngle = dotProductWithActiveAtomMask(dirPrev,dirVec,activeAtomMask)/(forceMag*dirMag);
        } else {
          LOG("something was zero length Using force" );
          copyVector(dirVec,forceVec);
          steepestDescent = true;
          cosAngle = 0.0;
        }
        this->lineSearch( &step, &fnew, pos, energyScale, dirVec, forceVec, tv1, tv2, localSteps, stepReport, activeAtomMask );
        prevStep = step;
        printedLatestMessage = false;
        if ( this->_PrintIntermediateResults ) {
          double tforceRmsMag = rmsMagnitudeWithActiveAtomMask(forceVec,activeAtomMask);
          printedLatestMessage = this->_displayIntermediateMessage(pos,step,fp,tforceRmsMag,cosAngle,steepestDescent,false,activeAtomMask);
        }
        if (this->_StepCallback.notnilp()) {
          dstep = core::DoubleFloat_O::create(step);
          core::eval::funcall(this->_StepCallback, _sym_steepest_descent, pos, forceVec, dstep, dirVec );
        }
        if (chem__verbose(5)) {
          core::clasp_write_string(fmt::format(" {}  lineSearch step-> {}\n" , __FUNCTION__ , step ));
        }
        if ( prevStep == 0.0 && step == 0.0 ) {
          ERROR(_sym_MinimizerStuck, (ql::list()
                                      << kw::_sym_minimizer << this->asSmartPtr()
                                      << kw::_sym_localSteps << core::clasp_make_fixnum(localSteps)
                                      << kw::_sym_coordinates << pos).result());
        }
        // Advance to new pos
        inPlaceAddTimesScalarWithActiveAtomMask( pos, dirVec, step, activeAtomMask );
        maybeSavePosition(pos);

		    // r = -f'(x)   r == force!!!!
        fp = dTotalEnergyForce( pos, energyScale, forceVec, activeAtomMask );
        if (callback.notnilp()) core::eval::funcall( callback, _sym_steepest_descent, this->_ScoringFunction, pos, activeAtomMask );
        if ( this->_DebugOn ) {
          this->stepReport(stepReport,fp,forceVec,activeAtomMask);
        }
#ifdef	VALIDATE_FORCE_ON //[
        this->validateForce(x,force);
#endif //]
		    // Don't use preconditioning
        copyVector(dirPrev,dirVec);
        copyVector(dirVec,forceVec);
        if ( this->_DebugOn ) {
          ASSERTNOTNULL(this->_Log);
          this->_Log->addReport(stepReport);
        }
        localSteps++;
        this->_Iteration++;
      }
      // Handle queued interrupts
      gctools::handle_all_queued_interrupts();
    }
  }

  if ( this->_PrintIntermediateResults && !printedLatestMessage ) {
    this->_displayIntermediateMessage(pos,step,fnew,forceRmsMag,cosAngle,steepestDescent,false,activeAtomMask);
  }
  fp = dTotalEnergyForce( pos, energyScale, forceVec, activeAtomMask );
  this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(pos,forceVec);
  LOG("Wrote coordinates and force to atoms" );
  if ( this->_DebugOn )
  {
    LOG("About to update stepReport" );
    ASSERTNOTNULL(stepReport);
    if ( stepReport.notnilp() )
    {
      LOG("About to add to stepReport" );
      stepReport->prematureTermination("Stuck");
      ASSERTNOTNULL(this->_Log);
      LOG("About to write report to debugLog" );
      this->_Log->addReport(stepReport);
      LOG("Done writing report to debugLog" );
    }
  }
  LOG("Leaving _steepestDescent" );
}



void	Minimizer_O::_conjugateGradient(int numSteps,
                                        NVector_sp pos,
                                        core::T_sp energyScale,
                                        double forceTolerance,
                                        core::T_sp activeAtomMask,
                                        core::T_sp callback )
{
  StepReport_sp	stepReport = StepReport_O::create();
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
  size_t        printedLatestMessage;

  if ( this->_PrintIntermediateResults ) {
    core::clasp_write_string(fmt::format("++++ Conjugate gradient ++++ max_steps = {}     tolerance = {}\n", numSteps, forceTolerance ));
  }
  if ( this->_Status == minimizerError ) return;
  this->_Status = conjugateGradientRunning;
  this->_CurrentPreconditioner = noPreconditioner;
  localSteps = 0;
  k = 0;
  refactor = 0;

    /* Calculate how many conjugate gradient steps can be */
    /* taken before a restart must be done */

  iRestartSteps = pos->size();
    // Define NVectors
  force = NVector_O::create(iRestartSteps);
  this->_Force = force;
  s = NVector_O::create(iRestartSteps);
  d = NVector_O::create(iRestartSteps);
  tv1 = NVector_O::create(iRestartSteps);
  tv2 = NVector_O::create(iRestartSteps);
    // Done
  innerSteps = MIN(iRestartSteps,ITMAX);
  double fp = dTotalEnergyForce( pos, energyScale, force, activeAtomMask );
//    r->inPlaceTimesScalar(-1.0);
    // TODO calculate preconditioner here
    // s = M^(-1)r rather than just copying it from r
  copyVector(s,force);
#if 0 //[
  switch ( preconditioner ) {
  case noPreconditioner:
      s->copy(force);
      break;
  case diagonalPreconditioner:
      diag = NVector_O::create(iRestartSteps);
      this->_EnergyFunction->setupDiagonalPreconditioner(x,diag);
      LOG("Preconditioner max value: {}" , diag->maxValue() );
      LOG("Preconditioner min value: {}" , diag->minValue() );
      this->_EnergyFunction->backSubstituteDiagonalPreconditioner(diag,s,force);
      break;
  case hessianPreconditioner:
      m = new_SparseLargeSquareMatrix_sp(iRestartSteps,SymmetricDiagonalLower);
      ldlt=new_SparseLargeSquareMatrix_sp(iRestartSteps,SymmetricDiagonalLower);
      m->fill(0.0);
      ldlt->fill(0.0);
      this->_EnergyFunction->setupHessianPreconditioner(x,m);
      this->_EnergyFunction->unconventionalModifiedCholeskyFactorization(m,ldlt);
      this->_EnergyFunction->backSubstituteLDLt(ldlt,s,force);
      break;
  default:
      SIMPLE_ERROR("Unknown preconditioner option");
  }
#endif //]
  copyVector(d,s);
  deltaNew = dotProductWithActiveAtomMask(force,d,activeAtomMask);
  delta0 = deltaNew;
  eSquaredDelta0 = forceTolerance*delta0;
  LOG("eSquaredDelta0 = {}" , (eSquaredDelta0 ) );
  LOG("forceTolerance = {}" , (forceTolerance ) );
  LOG("delta0 = {}" , (delta0 ) );
  localSteps = 0;
  step = 0.0;
  {
    while (1) {
	    //
	    // Absolute gradient test
	    //
      forceMag = magnitudeWithActiveAtomMask(force,activeAtomMask);
      forceRmsMag = rmsMagnitudeWithActiveAtomMask(force,activeAtomMask);
      this->_RMSForce = forceRmsMag;
      prevStep = step;
      printedLatestMessage = false;
      if ( this->_PrintIntermediateResults ) {
        printedLatestMessage = this->_displayIntermediateMessage(pos,prevStep,fp,forceRmsMag,cosAngle,steepestDescent,false,activeAtomMask);
      }
      if ( forceRmsMag < forceTolerance ) {
        if ( this->_PrintIntermediateResults ) {
          if (!printedLatestMessage) {
            printedLatestMessage = this->_displayIntermediateMessage(pos,step,fp,forceRmsMag,cosAngle,steepestDescent,true,activeAtomMask);
          }
          core::clasp_writeln_string(fmt::format(" ! DONE absolute force test:\n ! forceRmsMag({})<forceTolerance({})" , forceRmsMag , forceTolerance ));
        }
        break;
      }
      if ( localSteps > 0 ) {
        if ( localSteps>=numSteps ) {
          if ( this->_DebugOn )
          {
            if ( stepReport.notnilp() )
            {
              stepReport->prematureTermination("ExceededNumSteps");
              this->_Log->addReport(stepReport);
            }
          }
	//
	// Lets save the current conformation
	// before throwing this higher
	//
          fp = dTotalEnergyForce( pos, energyScale, force, activeAtomMask );
          this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(pos,force);
          ERROR(_sym_MinimizerExceededCG_MaxSteps, (ql::list() 
                                                 << kw::_sym_minimizer << this->asSmartPtr()
                                                 << kw::_sym_number_of_steps << core::make_fixnum(localSteps)
                                                 << kw::_sym_coordinates << pos).result());
        }
        if ( prevStep == 0.0 && step == 0.0 ) {
          if ( this->_DebugOn )
          {
            if ( stepReport.notnilp() )
            {
              stepReport->prematureTermination("Stuck");
              this->_Log->addReport(stepReport);
            }
          }
	//
	// Lets save the current conformation
	// before throwing this higher
	//
          fp = dTotalEnergyForce( pos, energyScale, force, activeAtomMask );
          this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(pos,force);
          MINIMIZER_STUCK_ERROR("Stuck in conjugate gradients",localSteps);
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
        dirMag = magnitudeWithActiveAtomMask(d,activeAtomMask);
        steepestDescent = false;
        LOG("Starting descent test" );
        cosAngle = 0.0;
        if ( forceMag != 0.0 && dirMag != 0.0 ) {
          LOG("forceMag = {}" , forceMag  );
          LOG("dirMag = {}" , dirMag  );
          cosAngle = dotProductWithActiveAtomMask(force,d,activeAtomMask)/(forceMag*dirMag);
        } else {
          LOG("some magnitude was zero Using force" );
          copyVector(d,force);
          steepestDescent = true;
        }

        this->lineSearch( &step, &fnew, pos, energyScale, d, force,
                          tv1, tv2, localSteps,
                          stepReport, activeAtomMask );

        if (this->_StepCallback.notnilp()) {
          core::DoubleFloat_sp dstep = core::DoubleFloat_O::create(step);
          core::eval::funcall(this->_StepCallback, _sym_conjugate_gradient, pos, force, dstep, d );
        }

		// x = x + (step)d
		// r = -f'(x)   r == force!!!!
        inPlaceAddTimesScalarWithActiveAtomMask( pos, d, step, activeAtomMask );
        maybeSavePosition(pos);

        fp = dTotalEnergyForce( pos, energyScale, force, activeAtomMask );
        if (callback.notnilp()) core::eval::funcall( callback, _sym_conjugate_gradient, this->_ScoringFunction, pos, activeAtomMask );

        if ( this->_DebugOn )
        {
          this->stepReport(stepReport,fp,force,activeAtomMask);
        }


#ifdef	VALIDATE_FORCE_ON
        this->validateForce(x,force);
#endif
        deltaOld = deltaNew;
        deltaMid = dotProductWithActiveAtomMask(force,s,activeAtomMask);

		// No preconditioning
        copyVector(s,force);

#if 0 //[
		// Calculate preconditioner M = f''(x)
		// s = M^(-1)r
        switch ( preconditioner ) {
        case noPreconditioner:
            s->copy(force);
            break;
        case diagonalPreconditioner:
            this->_EnergyFunction->setupDiagonalPreconditioner(x,diag);
            LOG("Preconditioner max value: {}" , diag->maxValue() );
            LOG("Preconditioner min value: {}" , diag->minValue() );
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
            SIMPLE_ERROR("Unknown preconditioner option");
        }
#endif //]
        deltaNew = dotProductWithActiveAtomMask(force,s,activeAtomMask);		// deltaNew = r.r
        beta = (deltaNew-deltaMid)/deltaOld;
        k = k + 1;
        if ( k == iRestartSteps || beta <= 0.0 ) {
          copyVector(d,s);
          k = 0;
          prevStep = 0.0;
        } else {
          XPlusYTimesScalarWithActiveAtomMask(d,s,d,beta,activeAtomMask);
        }
        if ( this->_DebugOn )
        {
          ASSERTNOTNULL(this->_Log);
          this->_Log->addReport(stepReport);
        }
        localSteps++;
        this->_Iteration++;
      }
      // Handle queued interrupts
      gctools::handle_all_queued_interrupts();
    }
  }
  if ( this->_PrintIntermediateResults && !printedLatestMessage ) {
    this->_displayIntermediateMessage( pos,step,fnew,forceRmsMag,cosAngle,steepestDescent,false,activeAtomMask);
  }
  fp = dTotalEnergyForce( pos, energyScale, force, activeAtomMask );
  this->_ScoringFunction->saveCoordinatesAndForcesFromVectors( pos,force );
  if ( this->_DebugOn )
  {
    if ( stepReport.notnilp() )
    {
      stepReport->prematureTermination("Stuck");
      this->_Log->addReport(stepReport);
    }
  }
}

void	Minimizer_O::_truncatedNewtonInnerLoop(int				kk,
                                               NVector_sp			pos,
                                               core::T_sp                       energyScale,
                                               SparseLargeSquareMatrix_sp	mprecon,
                                               SparseLargeSquareMatrix_sp	ldlt,
                                               NVector_sp 			force,
                                               double				rmsForceMag,
                                               NVector_sp			dirVec,
                                               NVector_sp			dirVecNext,
                                               NVector_sp			rj,
                                               NVector_sp			dj,
                                               NVector_sp			zj,
                                               NVector_sp			qj,
                                               bool&                            innerLoopDeltaJ1,
                                               core::T_sp activeAtomMask )
{
  int	jindex, ITpcg;
  double	cr, delta, forceDotpj, crOverk, nk, alphaj;
  double	nkTimesRmsForceMag, rjDotzj;
  double				djDotqj, forceDotpjNext;
  double				rmsRjMag, rjDotzjNext, betaj;
  gc::Nilable<NVector_sp>			nvDummy;
  gc::Nilable<SparseLargeSquareMatrix_sp>	nmDummy;

  ASSERTNOTNULL(this->_EnergyFunction);
  LOG("Resetting dummy vector and matrix" );
  nvDummy = nil<T_O>();
  nmDummy = nil<T_O>();

  if ( this->_DebugOn )
  {
    this->_Log->addMessage("_truncatedNewtonInnerLoop>>Starting\n");
  }
  LOG("Setting up" );
  jindex = 1;
  dirVec->zero();	// NVector
  cr = 0.5;
  delta = 1.0e-10;
  forceDotpj = 0.0;	// The initial value of force.pj
  ITpcg = 40;
  copyVector(rj,force);
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
  LOG("Carrying out UMC" );
//
//  Carry out UMC in the outer loop
//    this->_EnergyFunction->unconventionalModifiedCholeskyFactorization(mprecon,ldlt);
//
    // 2b.
    // Solve for zj in (M~)(zj)=(rj) by using the triangular
    // systems: (L)(x)=(rj) and (L^T)(zj)=(D^-1)(x)
    //
  LOG("Back substitute" );
  backSubstituteLDLt(ldlt,zj,rj);

    // 2c.
    // set dj = zj
    //
  copyVector(dj,zj);

  rjDotzj = dotProductWithActiveAtomMask(rj,zj,activeAtomMask);
  while ( 1 ) { 
    // 3. Singularity test
    // Compute the matrix-vector product qj=(H)(dj)
    // If either |((rj)^T).(zj)|<=delta
    // or |((dj)^T).(qj)|<=delta (e.g., delta=10^-10)
    // exit PCG loop with pk=pj ( for j=1, set pk=force)
    //

    this->_ScoringFunction->evaluateAll( pos,
                                         energyScale,
                                         nil<core::T_O>(),
                                         true, nvDummy,
                                         true, true, nmDummy,
                                         qj, dj, activeAtomMask);
    // MOVE rjDotzj calculation above this loop because
    // 	its calculated in step 6
    // rjDotzj = rj->dotProduct(zj);
    djDotqj = dotProductWithActiveAtomMask(dj,qj,activeAtomMask);
//    core::clasp_writeln_string(fmt::format( "{}:{} j=={} fabs(rjDotzj)={}  fabs(djDotqj)={} delta={}",__FILE__,__LINE__, jindex, fabs(rjDotzj),fabs(djDotqj),delta));
    if ( fabs(rjDotzj) <= delta || fabs(djDotqj) <= delta ) {
      if ( jindex==1 ) {
        copyVector(dirVec,force);
        innerLoopDeltaJ1 = true;
//        core::clasp_writeln_string(fmt::format( "{}:{} dirVec set to force",__FILE__,__LINE__));
      }
      if ( this->_DebugOn ) {
        stringstream ss;
        ss << "rjDotzj("<<rjDotzj<<").LT.delta("<<delta<<") or djDotqj("<<djDotqj<<").LT.delta("<<delta<<")"<<std::endl;
        this->_Log->addMessage(ss.str().c_str());
        this->_Log->addMessage("_truncatedNewtonInnerLoop>>Singularity test was true\n" );
      }
      LOG("Singularity test was true" );
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
    XPlusYTimesScalarWithActiveAtomMask(dirVecNext,dirVec,dj,alphaj,activeAtomMask);
    LOG("dirVecNext angle with force={}(deg)" , dirVecNext->angleWithVector(force)/0.0174533 );
    forceDotpjNext = dotProductWithActiveAtomMask(force,dirVecNext,activeAtomMask);
    if ( forceDotpjNext <= (forceDotpj + delta) ) {
      if ( jindex == 1 ) {
        copyVector(dirVec,force);
      } else {
	    // pk->copy(pj);  pk is pj
      }
      if ( this->_DebugOn ) {
        this->_Log->addMessage("_truncatedNewtonInnerLoop>>Descent direction test was true\n" );
      }
      LOG("Descent direction test was true" );
      goto DONE;
    }

    // 5. Truncation test
    // Compute rjNext = rj - (alphaj)(qj)
    // If "Truncation test" ||rjNext||<nk||g||
    // 	or "Step Limit" j+1>ITpcg then {
    // 	exit inner loop with search direction pk = dirVecNext
    // }
    //
    inPlaceAddTimesScalarWithActiveAtomMask(rj,qj,-alphaj,activeAtomMask);
    rmsRjMag = rmsMagnitudeWithActiveAtomMask(rj,activeAtomMask);
    if ( rmsRjMag < nkTimesRmsForceMag || (jindex+1)>ITpcg ) {
      LOG("rmsRjMag({}) < nkTimesRmsForceMag({})" , rmsRjMag , nkTimesRmsForceMag );
      copyVector(dirVec,dirVecNext);
      if ( this->_DebugOn ) {
        this->_Log->addMessage("_truncatedNewtonInnerLoop>>Truncation test was true\n" );
      }
      LOG("Truncation test was true" );
      goto DONE;
    }
    if ( (jindex+1)>ITpcg ) {
      LOG("j+1({})>ITpcg({})" , j+1 , ITpcg );
      copyVector(dirVec,dirVecNext);
      if ( this->_DebugOn ) {
        this->_Log->addMessage("_truncatedNewtonInnerLoop>>Step limit test was true\n" );
      }
      LOG("Step limit test was true" );
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
    rjDotzjNext = dotProductWithActiveAtomMask(rj,zj,activeAtomMask);
    betaj = rjDotzjNext/rjDotzj;
    rjDotzj = rjDotzjNext;
    XPlusYTimesScalarWithActiveAtomMask(dj, zj,dj,betaj,activeAtomMask);
    jindex = jindex + 1;
    copyVector(dirVec,dirVecNext);
  }
 DONE:
  LOG("Exiting inner loop with j = {}" , j );
  return;
}



#define	EPSILONF	1.0e-10
#define	EPSILONG	1.0e-8
#define	SQRT_EPSILONF	1.0e-5
#define	CUBERT_EPSILONF	4.6416e-4

void	Minimizer_O::_truncatedNewton(int numSteps,
                                      NVector_sp pos,
                                      core::T_sp energyScale,
                                      double forceTolerance,
                                      core::T_sp activeAtomMask,
                                      core::T_sp callback )
{
  StepReport_sp	stepReport = StepReport_O::create();
  int	iDimensions;
  double			fp;
  NVector_sp	forceK, dirVec, dirVecNext, rj, dj, zj, qj, posNext, kSum;
  SparseLargeSquareMatrix_sp	mprecon;
  SparseLargeSquareMatrix_sp    opt_mprecon;
  SparseLargeSquareMatrix_sp    ldlt, opt_ldlt;
  double			energyXk, energyXkNext;
  double			rmsForceMag;
  int				kk;
  double				alphaK, delta, rmsMagXKNext;
  bool				b1aTest, b1bTest;
  double		prevAlphaK = 0.0;
  double                dirMag, forceMag;
  double                cosAngle = 0.0;
  bool                  innerLoopDeltaJ1 = false;
#define	TENEMINUS8	10.0e-8

  if ( this->_PrintIntermediateResults ) {
    core::clasp_write_string(fmt::format("++++ Truncated Newton ++++ max_steps = {}     tolerance = {}\n", numSteps, forceTolerance ));
  }
  alphaK = 1.0;
  LOG("Starting TN" );
  if ( this->_Status == minimizerError ) return;
  this->_Status = truncatedNewtonRunning;
  this->_CurrentPreconditioner = this->_TruncatedNewtonPreconditioner;
  kk = 1;
    // Define NVectors
  LOG("Defining NVectors" );
  iDimensions = pos->size();
  forceK = NVector_O::create(iDimensions);
  this->_Force = forceK;
  LOG("Defining NVectors xKNext" );
  posNext = NVector_O::create(iDimensions);
  LOG("status" );
  dirVec = NVector_O::create(iDimensions);
  LOG("status" );
  dirVec->zero();
  LOG("Defining NVectors pjNext" );
  dirVecNext = NVector_O::create(iDimensions);
  LOG("Defining NVectors rj,dj,zj,qj" );
  rj = NVector_O::create(iDimensions);
  dj = NVector_O::create(iDimensions);
  zj = NVector_O::create(iDimensions);
  qj = NVector_O::create(iDimensions);
  kSum = NVector_O::create(iDimensions);
  mprecon = SparseLargeSquareMatrix_O::create(iDimensions,SymmetricDiagonalLower);
  ldlt = SparseLargeSquareMatrix_O::create(iDimensions,SymmetricDiagonalLower);
    //
    // Evaluate initial energy and force
    //
  LOG("Evaluating initial energy and force" );
  energyXkNext = dTotalEnergyForce( pos, energyScale, forceK, activeAtomMask );
  rmsForceMag = rmsMagnitudeWithActiveAtomMask(forceK,activeAtomMask);
  if ( this->_PrintIntermediateResults )
  {
    this->_displayIntermediateMessage(pos,prevAlphaK,energyXkNext,rmsForceMag,cosAngle,false,false,activeAtomMask);
  }

    //
    // Setup the preconditioner and carry out UMC
    //
  LOG("Setting up preconditioner" );
  this->_ScoringFunction->setupHessianPreconditioner(pos,mprecon,activeAtomMask);
  opt_mprecon = mprecon->optimized();
  // Insert entries once into ldlt based on mprecon so it has entries that can accept the factorization 
  unconventionalModifiedCholeskySymbolicFactorization(opt_mprecon,ldlt);
  // Factorize mprecon into ldlt
  unconventionalModifiedCholeskyFactorization(opt_mprecon,ldlt,kSum);
  opt_ldlt = ldlt->optimized();

#ifdef DEBUG_FACTORIZATION
  core::clasp_writeln_string(fmt::format( "======= Debug ldlt at start"));
  auto ldlt_debug = SparseLargeSquareMatrix_O::create(iDimensions,SymmetricDiagonalLower);
  chem_old::unconventionalModifiedCholeskySymbolicFactorization(opt_mprecon,ldlt_debug);
  auto kSum_debug = NVector_O::create(iDimensions);
  unconventionalModifiedCholeskyFactorization(opt_mprecon,ldlt_debug,kSum_debug);
  ldlt->ensureIdentical(ldlt_debug);
  chem__nvector_ensure_identical(kSum,kSum_debug, 0.01);
  if (this->_StepCallback.notnilp()) {
    core::eval::funcall(this->_StepCallback,_sym_truncated_newton_debug,opt_mprecon,ldlt,opt_ldlt);
  }
  if (callback.notnilp()) core::eval::funcall( callback, _sym_truncated_newton_debug, this->_ScoringFunction, x, activeAtomMask )
#endif


  {
    LOG("Starting loop" );
    while ( 1 ) {
      if ( this->_DebugOn ) {
        stepReport = StepReport_O::create();
        stepReport->_Iteration = this->_Iteration;
      }


	    //
	    // Inner loop
	    //
      _truncatedNewtonInnerLoop( kk, pos, energyScale, opt_mprecon, opt_ldlt,
                                 forceK, rmsForceMag, dirVec,
                                 dirVecNext, rj, dj, zj, qj,
                                 innerLoopDeltaJ1,
                                 activeAtomMask );

      // Added Feb 5, 2024 - if the inner loop fails to find a direction
      // and sets it to the gradient then truncated Newton will
      // stall and run out the max iterations
      // So lets use this as a convergence test
      if (innerLoopDeltaJ1) {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string(fmt::format( "search complete step {} according to innerLoopDeltaJ1" , this->_Iteration ));
        }
        break;
      }

	    //
	    // Line Search

      prevAlphaK = alphaK;
      if ( this->_PrintIntermediateResults ) {
        dirMag = magnitudeWithActiveAtomMask(dirVec,activeAtomMask);
        forceMag = magnitudeWithActiveAtomMask(forceK,activeAtomMask);
//        core::clasp_writeln_string(fmt::format( "  dirMag = {}   forceMag = {} \n", dirMag, forceMag ));
        LOG("Starting descent test" );
        if ( forceMag != 0.0 && dirMag != 0.0 ) {
          LOG("forceMag = {}" , forceMag  );
          LOG("dirMag = {}" , dirMag  );
          cosAngle = dotProductWithActiveAtomMask(forceK,dirVec,activeAtomMask)/(forceMag*dirMag);
        } else {
          cosAngle = 0.0;
        }
      }

      energyXk = energyXkNext;
      this->lineSearch( &alphaK, &energyXkNext, pos, energyScale, dirVec, forceK, zj, qj, kk, stepReport, activeAtomMask );
      if (this->_StepCallback.notnilp()) {
        core::DoubleFloat_sp dstep = core::DoubleFloat_O::create(alphaK);
        core::eval::funcall(this->_StepCallback, _sym_truncated_newton, pos, forceK, dstep, dirVec, opt_mprecon, opt_ldlt  );
      }
      if (callback.notnilp()) core::eval::funcall( callback, _sym_truncated_newton, this->_ScoringFunction, pos, activeAtomMask );
      XPlusYTimesScalarWithActiveAtomMask(posNext, pos, dirVec, alphaK,activeAtomMask);
      maybeSavePosition(posNext);
	    //
	    // Evaluate the force at the new position
	    //
      fp = dTotalEnergyForce( posNext, energyScale, forceK, activeAtomMask );
      if ( this->_DebugOn )
      {
        this->stepReport(stepReport,fp,forceK,activeAtomMask);
      }

	    //
	    // Convergence tests
	    //

      b1aTest=fabs(energyXkNext-energyXk)<EPSILONF*(1.0+fabs(energyXk));
      LOG("b1aTest = {}" , b1aTest  );
      LOG("energyXkNext({})" , energyXkNext );
      LOG("energyXk({})" , energyXk );
      LOG("fabs[energyXkNext-energyXk]={}" , fabs(energyXkNext-energyXk) );
      LOG("[ EPSILONF*(1.0+fabs(energyXk))={}" , EPSILONF*(1.0+fabs(energyXk)) );
      if ( b1aTest ) {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string(fmt::format( "search complete step {} according to b1aTest" , this->_Iteration ));
          core::clasp_writeln_string(fmt::format("last alphaK = {} energyXkNext = {} energyXk = {}" , alphaK , energyXkNext , energyXk ));
          core::clasp_writeln_string(fmt::format("b1aTest=fabs(energyXkNext-energyXk)<EPSILONF*(1.0+fabs(energyXk)) {} < {};" , fabs(energyXkNext-energyXk) , (EPSILONF*(1.0+fabs(energyXk)))));
        }
        break;
      }
      delta = rmsDistanceFromWithActiveAtomMask(posNext,pos,activeAtomMask);
      rmsMagXKNext = rmsMagnitudeWithActiveAtomMask(posNext,activeAtomMask);
      b1bTest=(delta<SQRT_EPSILONF*(1.0+rmsMagXKNext)/100.0);
      if ( b1bTest ) {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string(fmt::format( "search complete step {} according to b1bTest delta[{}] < SQRT_EPSILONF*(1.0+rmsMagXKNext)/100.0 [{}]"  , this->_Iteration , delta , (SQRT_EPSILONF*(1.0+rmsMagXKNext)/100.0)));
        }
        break;
      }

      rmsForceMag = rmsMagnitudeWithActiveAtomMask(forceK,activeAtomMask);
      if ( rmsForceMag < forceTolerance ) {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string(fmt::format( "search complete [{}] according to absolute force test rmsForceMag({}) < forceTolerance({})"  , this->_Iteration , rmsForceMag , forceTolerance ));
        }
        break;
      }
	    //
	    // Preparation for next Newton step
	    //
	    // Compute the preconditioner M at X{k+1}
	    //
      this->_ScoringFunction->setupHessianPreconditioner(pos,mprecon,activeAtomMask);
      opt_mprecon = mprecon->optimized();
      unconventionalModifiedCholeskyFactorization(opt_mprecon,ldlt,kSum);
      opt_ldlt = ldlt->optimized();
#ifdef DEBUG_FACTORIZATION
  core::clasp_writeln_string(fmt::format( "======= Debug ldlt in iteration"));
  auto ldlt_debug = SparseLargeSquareMatrix_O::create(iDimensions,SymmetricDiagonalLower);
  chem_old::unconventionalModifiedCholeskySymbolicFactorization(opt_mprecon,ldlt_debug);
  auto kSum_debug = NVector_O::create(iDimensions);
  unconventionalModifiedCholeskyFactorization(opt_mprecon,ldlt_debug,kSum_debug);
  if (this->_StepCallback.notnilp()) {
    core::eval::funcall(this->_StepCallback,_sym_truncated_newton_debug,opt_mprecon,ldlt,opt_ldlt);
  }
  ldlt->ensureIdentical(ldlt_debug);
  chem__nvector_ensure_identical(kSum,kSum_debug, 0.01);
#endif

      copyVector(pos,posNext);
      kk++;
      this->_Iteration++;
      if ( this->_DebugOn )
      {
        ASSERTNOTNULL(this->_Log);
        this->_Log->addReport(stepReport);
      }
      if ( kk > numSteps ) {
        if ( this->_DebugOn )
        {
          if ( stepReport.notnilp() )
          {
            stepReport->prematureTermination("ExceededNumSteps");
            this->_Log->addReport(stepReport);
          }
        }
	//
	// Lets save the current conformation
	// before throwing this higher
	//
        dTotalEnergyForce( pos, energyScale, forceK, activeAtomMask );
        this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(pos,forceK);
        ERROR(_sym_MinimizerExceededTN_MaxSteps, (ql::list() 
                                               << kw::_sym_minimizer << this->asSmartPtr()
                                               << kw::_sym_number_of_steps << core::make_fixnum(kk)
                                               << kw::_sym_coordinates << pos).result());
      }
      if ( this->_PrintIntermediateResults ) {
        this->_displayIntermediateMessage(pos,prevAlphaK,fp,rmsForceMag,cosAngle,false,false,activeAtomMask);
      }

// Handle queued interrupts
      gctools::handle_all_queued_interrupts();
    }
  }
  copyVector(pos,posNext);
  dTotalEnergyForce( pos, energyScale, forceK, activeAtomMask );
  this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(pos,forceK);
  if ( this->_DebugOn )
  {
    if ( stepReport.notnilp() )
    {
      stepReport->prematureTermination("Stuck");
      this->_Log->addReport(stepReport);
    }
  }
}


/*
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------

*/

void	Minimizer_O::_evaluateEnergyAndForceManyTimes( int numSteps,  NVector_sp nvPos, core::T_sp energyScale, core::T_sp activeAtomMask )
{
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
  nvDir       = NVector_O::create( iSize );
  nvNewPos    = NVector_O::create( iSize );
  nvNewForce = NVector_O::create( iSize );
  nvTempPos   = NVector_O::create( iSize );
  nvTempForce= NVector_O::create( iSize );

  iCount = 0;
  this->_Iteration = 1;
  do {
    dEnergy = this->dTotalEnergyForce( nvPos, energyScale, nvNewForce, activeAtomMask );
    if ( iCount % 10000 == 0 ) {
      core::clasp_writeln_string(fmt::format("Evaluating energy step#{}" , iCount ));
    }
    iCount++;
  } while ( iCount < numSteps );

  LOG("Exceeded maximum number of steps" );

}


void	Minimizer_O::validateForce(NVector_sp pos, core::T_sp energyScale, NVector_sp force, core::T_sp activeAtomMask )
{
  ForceMatchReport_sp	report;
  if ( this->_DebugOn ) {
    report = this->_ScoringFunction->checkIfAnalyticalForceMatchesNumericalForce(pos,energyScale,force,activeAtomMask);
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
  this->_Position = nil<core::T_O>();
  this->_Force = nil<core::T_O>();
//	this->_StepCallback = _Nil<core::LispCallback_O>();
}



CL_LISPIFY_NAME("useDefaultSettings");
CL_DEFMETHOD     void	Minimizer_O::useDefaultSettings()
{
  this->_InitialLineSearchStep = 0.01;
  this->_NumberOfSteepestDescentSteps = MAXSTEEPESTDESCENTSTEPS;
  this->_SteepestDescentTolerance = 2000.0;
  this->_NumberOfConjugateGradientSteps = MAXCONJUGATEGRADIENTSTEPS;
  this->_ConjugateGradientTolerance = 10.0;		//	Use this for now, later add TN minimizer and switch to that when this is <10.0
  this->_NumberOfTruncatedNewtonSteps = MAXTRUNCATEDNEWTONSTEPS;
  this->_TruncatedNewtonTolerance = 0.00000001;
  this->_TruncatedNewtonPreconditioner = hessianPreconditioner;
  this->_PrintIntermediateResults = 0;
  LOG("_PrintIntermediateResults = {}" , this->_PrintIntermediateResults  );
  this->_ReportEveryNSteps = 100;
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



CL_LISPIFY_NAME(minimizer-set-debug-on);
CL_DEFMETHOD void Minimizer_O::setDebugOn(bool debugOn) {
  this->_DebugOn = debugOn;
  if (debugOn) {
    this->_Log = MinimizerLog_O::create();
  }
}


CL_LISPIFY_NAME("enablePrintIntermediateResults");
CL_LAMBDA((chem:minimizer chem:minimizer) cl:&optional (steps 1) (level 1));
CL_DEFMETHOD     void	Minimizer_O::enablePrintIntermediateResults(size_t steps, size_t level)
{
  if (steps < 1 ) {
    this->_ReportEveryNSteps = 100;
  } else {
    this->_ReportEveryNSteps = steps;
  }
  this->_PrintIntermediateResults = level;
}




CL_LISPIFY_NAME("disablePrintIntermediateResults");
CL_DEFMETHOD     void	Minimizer_O::disablePrintIntermediateResults()
{
  this->_PrintIntermediateResults = 0;
}




CL_LISPIFY_NAME("setEnergyFunction");
CL_DEFMETHOD     void	Minimizer_O::setEnergyFunction(ScoringFunction_sp f)
{
  this->_ScoringFunction = f;
  this->_Iteration = 1;
}

CL_LISPIFY_NAME("evaluateEnergyAndForceManyTimes");
CL_DEFMETHOD     void	Minimizer_O::evaluateEnergyAndForceManyTimes(core::T_sp energyScale, int numSteps, core::T_sp activeAtomMask)
{
  NVector_sp	pos;
  ASSERT(this->_ScoringFunction);
  this->_Iteration = 1;
  pos = NVector_O::create(this->_ScoringFunction->getNVectorSize());
  this->_ScoringFunction->loadCoordinatesIntoVector(pos);
  this->_evaluateEnergyAndForceManyTimes(numSteps,pos,energyScale,activeAtomMask);
}



CL_LISPIFY_NAME("resetAndMinimize");
CL_LAMBDA((minimizer chem:minimizer) &key energy-scale active-atom-mask callback );
CL_DEFMETHOD     void	Minimizer_O::resetAndMinimize(core::T_sp energyScale, core::T_sp activeAtomMask, core::T_sp callback )
{
  this->_Status = minimizerIdle;
  this->minimize(energyScale, activeAtomMask, callback );
}


CL_LISPIFY_NAME("minimize");
CL_LAMBDA((minimizer chem:minimizer) &key energy-scale active-atom-mask callback coords);
CL_DEFMETHOD core::T_mv Minimizer_O::minimize(core::T_sp energyScale, core::T_sp activeAtomMask, core::T_sp callback, core::T_sp coords)
{
  NVector_sp	pos;
  int		retries;
  this->_StartTime._value = std::chrono::steady_clock::now();
  ASSERT(this->_ScoringFunction);
  if (coords.nilp()) {
    pos = NVector_O::create(this->_ScoringFunction->getNVectorSize());
  } else if (gc::IsA<NVector_sp>(coords)) {
    // Use a copy of the coords so we don't overwrite the input
    pos = copy_nvector(gc::As_unsafe<NVector_sp>(coords),0,nil<core::T_O>());
    if (this->_ScoringFunction->getNVectorSize()!=pos->size()) {
      SIMPLE_ERROR("You provided coordinates but they have the wrong length {} elements - expected {} elements", pos->size(), this->_ScoringFunction->getNVectorSize());
    }
  } else {
    SIMPLE_ERROR("Coordinates must either be NIL or an nvector of length {}", this->_ScoringFunction->getNVectorSize());
  }
  this->_Position = pos;
  int maxRetries = 100;
  retries = maxRetries;
  do {
    try {
      // If coords were provided then don't read them out of the molecule.
      if (coords.nilp()) this->_ScoringFunction->loadCoordinatesIntoVector(pos);
      if (chem__verbose(4)) {
        for (size_t idx=0; idx<pos->length(); idx++ ) {
          core::clasp_write_string(fmt::format("Starting pos[{}] -> {}\n" , idx , (*pos)[idx]));
        }
      }
      if ( this->_NumberOfSteepestDescentSteps > 0 ) {
        this->_steepestDescent( this->_NumberOfSteepestDescentSteps,
                                pos, energyScale, this->_SteepestDescentTolerance, activeAtomMask, callback );
      } else {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string("======= Skipping Steepest Descent #steps = 0");
        }
      }
      if ( this->_NumberOfConjugateGradientSteps > 0 ) {
        this->_conjugateGradient( this->_NumberOfConjugateGradientSteps,
                                  pos, energyScale, this->_ConjugateGradientTolerance, activeAtomMask, callback );
      } else {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string("======= Skipping Conjugate Gradients #steps = 0");
        }
      }
      if ( this->_NumberOfTruncatedNewtonSteps > 0 ) {
        this->_truncatedNewton( this->_NumberOfTruncatedNewtonSteps,
                                pos, energyScale, this->_TruncatedNewtonTolerance, activeAtomMask, callback );
      } else {
        if ( this->_PrintIntermediateResults ) {
          core::clasp_writeln_string("======= Skipping Truncated Newton #steps = 0");
        }
      }
      if ( this->_PrintIntermediateResults ) {
        core::clasp_writeln_string("======= All three minimizers have completed or passed on minimization.");
      }
      goto DONE;
    } catch ( RestartMinimizer ld ) {
      retries--;
#if 0
      ERROR(_sym_MinimizerError, (ql::list()
                                  << kw::_sym_message << core::SimpleBaseString_O::make("RestartMinimizer")
                                  << kw::_sym_minimizer << this->asSmartPtr()
                                  << kw::_sym_coordinates << pos).result());
#endif
    }
  } while ( retries > 0 );
  ERROR(_sym_MinimizerError, (ql::list()
                              << kw::_sym_message << core::SimpleBaseString_O::make("too-many-RestartMinimizer")
                              << kw::_sym_minimizer << this->asSmartPtr()
                              << kw::_sym_coordinates << pos).result());
 DONE:
  return Values(pos,mk_double_float(this->dTotalEnergy(pos,energyScale,activeAtomMask)));
}


CL_LISPIFY_NAME("writeIntermediateResultsToEnergyFunction");
CL_DEFMETHOD void Minimizer_O::writeIntermediateResultsToEnergyFunction()
{
  if ( this->_Position.nilp() ) {
    SIMPLE_ERROR("There are no intermediate results");
  }
  this->_ScoringFunction->saveCoordinatesAndForcesFromVectors(this->_Position,this->_Force);
}

#if 0
adapt::QDomNode_sp	Minimizer_O::asXml()
{
  adapt::QDomNode_sp	xml;
  xml = adapt::QDomNode_O::create("Minimizer");
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
#endif


CL_LISPIFY_NAME("configurationAsString");
CL_DEFMETHOD     string	Minimizer_O::configurationAsString()
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

CL_LISPIFY_NAME("minimizer-restart");
CL_DEFMETHOD     void	Minimizer_O::restart()
{
  this->_Status = minimizerIdle;
  this->_Iteration = 1;
}


DOCGROUP(cando);
CL_DEFUN void chem__restart_minimizer()
{
  throw RestartMinimizer();
}


SYMBOL_EXPORT_SC_(ChemPkg,noPreconditioner);
SYMBOL_EXPORT_SC_(ChemPkg,hessianPreconditioner);
CL_BEGIN_ENUM(chem::PreconditionerType,_sym__PLUS_preconditionerTypeConverter_PLUS_,"PreconditionerType");
CL_VALUE_ENUM(_sym_noPreconditioner,noPreconditioner);
CL_VALUE_ENUM(_sym_hessianPreconditioner,hessianPreconditioner);
CL_END_ENUM(_sym__PLUS_preconditionerTypeConverter_PLUS_);
;
SYMBOL_EXPORT_SC_(ChemPkg,minimizerError);
SYMBOL_EXPORT_SC_(ChemPkg,minimizerSucceeded);
SYMBOL_EXPORT_SC_(ChemPkg,truncatedNewtonRunning);
SYMBOL_EXPORT_SC_(ChemPkg,conjugateGradientRunning);
SYMBOL_EXPORT_SC_(ChemPkg,steepestDescentRunning);
SYMBOL_EXPORT_SC_(ChemPkg,minimizerIdle);
CL_BEGIN_ENUM(chem::MinimizerStatus,_sym__PLUS_minimizerStatusConverter_PLUS_,"MinimizerStatus");
CL_VALUE_ENUM(_sym_minimizerIdle,chem::minimizerIdle);
CL_VALUE_ENUM(_sym_steepestDescentRunning,chem::steepestDescentRunning);
CL_VALUE_ENUM(_sym_conjugateGradientRunning,chem::conjugateGradientRunning);
CL_VALUE_ENUM(_sym_truncatedNewtonRunning,chem::truncatedNewtonRunning);
CL_VALUE_ENUM(_sym_minimizerSucceeded,chem::minimizerSucceeded);
CL_VALUE_ENUM(_sym_minimizerError,chem::minimizerError);
CL_END_ENUM(_sym__PLUS_minimizerStatusConverter_PLUS_);


};
