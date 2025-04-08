/*
    File: minimizer.h
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
 *	minimizer.h
 *
 *	Maintains a database of stretch types
 */

#ifndef FFMinimizer_H
#define	FFMinimizer_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <chrono>
#include <clasp/core/common.h>
//#include "bond.h"
#include <cando/geom/vector3.h>
#include <cando/chem/matter.fwd.h>


#include <cando/adapt/quickDom.fwd.h>// minimizer.h wants QDomNode needs quickDom.fwd.h
//#include "core/lispCallback.fwd.h"// minimizer.h wants LispCallback needs lispCallback.fwd.h

namespace       chem
{


//SMART(QDomNode);
  SMART(Atom);
  SMART(Residue);
  SMART(ForceField);
  SMART(MinimizerLog);
  SMART(ScoringFunction);
  SMART(SparseLargeSquareMatrix);
  SMART(StepReport);

#define	CALLBACK_ERROR	0
#define	CALLBACK_CONT	1
#define	CALLBACK_HALT	2


  typedef	enum { minimizerIdle,
		       steepestDescentRunning,
		       conjugateGradientRunning,
		       truncatedNewtonRunning,
		       minimizerSucceeded, minimizerError } MinimizerStatus;


  typedef	enum	{ unknownPreconditioner, noPreconditioner, hessianPreconditioner } PreconditionerType;

  extern	string	stringForPreconditionerType(PreconditionerType t);
  extern	string	shortStringForPreconditionerType(PreconditionerType t);
  extern	PreconditionerType	preconditionerTypeFromString(const string& t );


struct RestartMinimizer {};

  struct MinimizerCondition  {
  };

  struct MinimizerCondition_ExceededNumSteps : public MinimizerCondition {
    core::Symbol_sp    _Minimizer;
    size_t _NumberOfSteps;
  MinimizerCondition_ExceededNumSteps(core::Symbol_sp min, size_t steps)
    : _Minimizer(min),
      _NumberOfSteps(steps) {};
  };



  class   MinimizerCondition_Stuck
  {
  private:
    string _Message;
  public:
    string  message() { return "Minimizer is stuck!("+ this->_Message+")"; };
    MinimizerCondition_Stuck(const string& msg) {this->_Message = msg;};
    virtual ~MinimizerCondition_Stuck() throw() {};
  };

  SMART(Minimizer);
  class Minimizer_O : public core::CxxObject_O
  {
    LISP_CLASS(chem,ChemPkg,Minimizer_O,"Minimizer",core::CxxObject_O);
  public:
    static Minimizer_sp make(ScoringFunction_sp givenEnergyFunction);
  public:
    void initialize();
#if 0
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
#endif
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  private:
	// behavior
    double		_InitialLineSearchStep;

    double		_SteepestDescentTolerance;
    int			_NumberOfSteepestDescentSteps;

    int			_NumberOfConjugateGradientSteps;
    double		_ConjugateGradientTolerance;

    int			_NumberOfTruncatedNewtonSteps;
    double		_TruncatedNewtonTolerance;
    PreconditionerType	_TruncatedNewtonPreconditioner;

	// status
    bool		_DebugOn;
    MinimizerLog_sp	_Log;
    MinimizerStatus	_Status;
    PreconditionerType	_CurrentPreconditioner;
    int			_MinBracketSteps;
    size_t		_PrintIntermediateResults;
    int			_ReportEveryNSteps;
    ScoringFunction_sp	_ScoringFunction;
    int			_Iteration;
    dont_expose<std::chrono::time_point<std::chrono::steady_clock>>	_StartTime;
    bool		_ShowElapsedTime;
    double		_MinGradientMean;
    double		_RMSForce;
    NVector_sp		nvP1DSearchTemp1;
    NVector_sp		nvP1DSearchTemp2;
    NVector_sp		nvP1DSearchOrigin;
    NVector_sp		nvP1DSearchDirection;
        // Save coordinates and forces so that if we interrupt the minimization we can recover them
    gc::Nilable<NVector_sp> _Position;
    gc::Nilable<NVector_sp> _Force;

  private:	// Do not serialize
    core::T_sp          _StepCallback;

  public:
    void set_frozen(core::T_sp frozen);
    void set_initial_line_search_step(double step);
    
    void lineSearchInitialReport( StepReport_sp report,
                                  NVector_sp nvPos,
                                  core::T_sp energyScale,
                                  NVector_sp nvDir,
                                  NVector_sp nvForce,
                                  double xa, double xx, double xb,
                                  double fa, double fx, double fb,
                                  core::T_sp activeAtomMask );
    void lineSearchFinalReport( StepReport_sp report, double step, double fmin,
                                int energyEvals, int forceEvals, int dbrentSteps );
    void stepReport( StepReport_sp report, double energy,NVector_sp force,core::T_sp activeAtomMask );

    void	getPosition(NVector_sp nvResult, NVector_sp nvOrigin, NVector_sp nvDirection, double x, core::T_sp activeAtomMask);
    double dTotalEnergy(NVector_sp pos,core::T_sp energyScale,core::T_sp activeAtomMask);
    double dTotalEnergyForce( NVector_sp nvPos, core::T_sp energyScale, NVector_sp nvForce, core::T_sp activeAtomMask);
    double d1DTotalEnergy(double x, core::T_sp energyScale, core::T_sp activeAtomMask );
    double d1DTotalEnergyForce( double x, core::T_sp energyScale, double* fx, double* dfx, core::T_sp activeAtomMask );
    void   minBracket(	NVector_sp	nvOrigin,
                        core::T_sp      energyScale,
                        NVector_sp	nvDir,
                        double		*dPxa,
                        double		*dPxb,
                        double		*dPxc,
                        double		*dPfa,
                        double		*dPfb,
                        double		*dPfc,
                        core::T_sp activeAtomMask );
    double 	dbrent(	double ax, double bx, double cx,
                        core::T_sp energyScale,
			double tol, double& step,
			int&	energyEvals,
			int&	forceEvals,
			int&	dbrentSteps,
                        core::T_sp activeAtomMask );
    void	lineSearch(	double	*dPxnew,
				double	*dPfnew,
				NVector_sp nvOrigin,
                                core::T_sp energyScale,
				NVector_sp nvDirection,
				NVector_sp nvForce,
				NVector_sp nvTemp1,
				NVector_sp nvTemp2,
	       			int	iteration,
				StepReport_sp report,
                                core::T_sp activeAtomMask );
    void		define1DSearch( NVector_sp o, NVector_sp d,
					NVector_sp t1, NVector_sp t2 ) {
      this->nvP1DSearchOrigin = o;
      this->nvP1DSearchDirection = d;
      this->nvP1DSearchTemp1 = t1;
      this->nvP1DSearchTemp2 = t2; };

        /*! Return true on success */
    void	_steepestDescent( int numSteps,
				  NVector_sp p,
                                  core::T_sp energyScale,
				  double rmsGradientTol,
                                  core::T_sp activeAtomMask,
                                  core::T_sp callback );
        /*! Return true on success */
    void	_conjugateGradient( int numSteps,
				    NVector_sp p,
                                  core::T_sp energyScale,
				    double rmsGradientTol,
                                    core::T_sp activeAtomMask,
                                    core::T_sp callback );
    void	_truncatedNewton( int numSteps,
				  NVector_sp p,
                                  core::T_sp energyScale,
				  double rmsGradientTol,
                                  core::T_sp activeAtomMask,
                                  core::T_sp callback );


    void	_evaluateEnergyAndForceManyTimes( int numSteps,
                                                  NVector_sp nvPos,
                                                  core::T_sp energyScale,
                                                  core::T_sp activeAtomMask );
    void	validateForce(NVector_sp pos,
                              core::T_sp energyScale,
                              NVector_sp force,
                              core::T_sp activeAtomMask );

    void	debugBeginIteration(int i);
//	void	debugAdd(adapt::QDomNode_sp	node);
    void	debugAccumulateMessage( const char* str );
    void	debugEndIteration();


    bool _displayIntermediateMessage( NVector_sp        pos,
                                      double		step,
                                      double		fnew,
                                      double		forceRmsMag,
                                      double		cosAngle,
                                      bool		steepestDescent,
                                      bool              forcePrint = false,
                                      core::T_sp        activeAtomMask = nil<core::T_O>() );

  public:	// should be private
    void _truncatedNewtonInnerLoop(
                                   int				k,
                                   NVector_sp			xj,
                                   core::T_sp                   energyScale,
                                   SparseLargeSquareMatrix_sp	mprecon,
                                   SparseLargeSquareMatrix_sp	ldlt,
                                   NVector_sp 			force,
                                   double			rmsForceMag,
                                   NVector_sp			pj,
                                   NVector_sp			pjNext,
                                   NVector_sp			rj,
                                   NVector_sp			dj,
                                   NVector_sp			zj,
                                   NVector_sp			qj,
                                   bool& failedInnerLoop,
                                   core::T_sp activeAtomMask );




  public:


    void throwMinimizerExceededMaxSteps();
    void throwMinimizerStuck();
    void throwMinimizerError();

    CL_LISPIFY_NAME("getEnergyFunction");
    CL_DEFMETHOD 	ScoringFunction_sp	getEnergyFunction() {_OF(); ASSERTNOTNULL(this->_ScoringFunction); return this->_ScoringFunction; };
    void	useDefaultSettings();

    void	changeOptions(core::List_sp options);

    MinimizerStatus status()	{return this->_Status;};
    CL_LISPIFY_NAME("status");
    CL_DEFMETHOD 	int statusAsInt()	{return (int)(this->_Status);};

    CL_LISPIFY_NAME("getIteration");
    CL_DEFMETHOD 	int	getIteration()	{ return this->_Iteration; };
    string	statusAsString();
    string	statusAsShortString();

    CL_LISPIFY_NAME("getRMSForce");
    CL_DEFMETHOD 	double	getRMSForce() { return this->_RMSForce; };

    CL_DOCSTRING(R"dx(Set a callback that is called at every minimization step. The callback function
takes a single argument, the NVECTOR position of the atoms.)dx");
    CL_DEFMETHOD void	setStepCallback(core::T_sp stepCallback) { this->_StepCallback = stepCallback;};

    CL_LISPIFY_NAME("setMaximumNumberOfSteepestDescentSteps");
    CL_DEFMETHOD 	void	setMaximumNumberOfSteepestDescentSteps(int m) {this->_NumberOfSteepestDescentSteps = m;};
    CL_LISPIFY_NAME("setMaximumNumberOfConjugateGradientSteps");
    CL_DEFMETHOD 	void	setMaximumNumberOfConjugateGradientSteps(int m) {this->_NumberOfConjugateGradientSteps= m;};
    CL_LISPIFY_NAME("setMaximumNumberOfTruncatedNewtonSteps");
    CL_DEFMETHOD 	void	setMaximumNumberOfTruncatedNewtonSteps(int m) {this->_NumberOfTruncatedNewtonSteps= m;};
    CL_LISPIFY_NAME("setSteepestDescentTolerance");
    CL_DEFMETHOD 	void	setSteepestDescentTolerance(double m) {this->_SteepestDescentTolerance = m;};
    CL_LISPIFY_NAME("setConjugateGradientTolerance");
    CL_DEFMETHOD 	void	setConjugateGradientTolerance(double m) {this->_ConjugateGradientTolerance = m;};
    CL_LISPIFY_NAME("setTruncatedNewtonTolerance");
    CL_DEFMETHOD 	void	setTruncatedNewtonTolerance(double m) {this->_TruncatedNewtonTolerance = m;};
    void	setTruncatedNewtonPreconditioner(PreconditionerType m) {this->_TruncatedNewtonPreconditioner = m;};


    void	setEnergyFunction(ScoringFunction_sp ef);
    void	enablePrintIntermediateResults(size_t steps, size_t level);
    void	disablePrintIntermediateResults();
    void        setDebugOn(bool debugOn);
    void        setReportEveryNSteps(size_t steps);
    size_t      getReportEveryNSteps();
    
    string	configurationAsString();
    void	restart();

    void	minimizeSteepestDescent();
    void	minimizeConjugateGradient();
    void	resetAndMinimize(core::T_sp energyScale, core::T_sp activeAtomMask, core::T_sp callback );
    core::T_mv minimize(core::T_sp energyScale, core::T_sp activeAtomMask, core::T_sp callback, core::T_sp coordinates=nil<core::T_O>() );
        // If the minimization is aborted the intermediate results can be recovered
    void    writeIntermediateResultsToEnergyFunction();
    void	evaluateEnergyAndForceManyTimes(core::T_sp energyScale, int times,core::T_sp activeAtomMask );

    adapt::QDomNode_sp	asXml();

    Minimizer_O() :
        _PrintIntermediateResults(1),
        _StepCallback(nil<core::T_O>())
    {};
  };

  extern core::Symbol_sp& _sym__PLUS_minimizerStatusConverter_PLUS_;
  extern core::Symbol_sp& _sym__PLUS_preconditionerTypeConverter_PLUS_;

};

DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::MinimizerStatus,chem::_sym__PLUS_minimizerStatusConverter_PLUS_);
DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::PreconditionerType,chem::_sym__PLUS_preconditionerTypeConverter_PLUS_);


#endif
