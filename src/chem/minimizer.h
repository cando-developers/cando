       
       


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
#include "core/common.h"
//#include "bond.h"
#include "addon/vector3.h"
#include "matter.fwd.h"


#include "addon/quickDom.fwd.h"// minimizer.h wants QDomNode needs quickDom.fwd.h
#include "core/posixTime.fwd.h"// minimizer.h wants PosixTime needs posixTime.fwd.h
//#include "core/lispCallback.fwd.h"// minimizer.h wants LispCallback needs lispCallback.fwd.h

namespace       chem
{





//SMART(QDomNode);
    SMART(Atom);
    SMART(Residue);
    SMART(NVector);
    SMART(ForceField);
    SMART(MinimizerLog);
    SMART(EnergyFunction);
    SMART(SparseLargeSquareMatrix);
    SMART(StepReport);
    SMART(PosixTime);
    SMART(LispCallback);

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


    class	MinimizerCondition_ExceededNumSteps
    {
    public:
	string _Message;
	int	_NumberOfSteps;
	string	message() {stringstream ss; ss.str(""); ss<<"Exceeded max("<<this->_NumberOfSteps<<") steps in "<< this->_Message; return ss.str(); };

	MinimizerCondition_ExceededNumSteps() {};
	virtual ~MinimizerCondition_ExceededNumSteps() throw() {};
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
    class Minimizer_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,Minimizer_O,"Minimizer");
#if INIT_TO_FACTORIES
    public:
	static Minimizer_sp make(Matter_sp matter, ForceField_sp forceField, EnergyFunction_sp givenEnergyFunction);
#else
	DECLARE_INIT();
#endif
    public:
	static void lisp_initGlobals(core::Lisp_sp lisp);
    public:
	void initialize();
    public:
//	void	archiveBase(core::ArchiveP node);

    private:
	// behavior
	double			_InitialLineSearchStep;

	double			_SteepestDescentTolerance;
    	int			_NumberOfSteepestDescentSteps;

    	int			_NumberOfConjugateGradientSteps;
	double			_ConjugateGradientTolerance;

    	int			_NumberOfTruncatedNewtonSteps;
	double			_TruncatedNewtonTolerance;
	PreconditionerType	_TruncatedNewtonPreconditioner;

	// status
	bool			_DebugOn;
	MinimizerLog_sp		_Log;
	MinimizerStatus		_Status;
	PreconditionerType	_CurrentPreconditioner;
	stringstream		_Message;
    	int			_MinBracketSteps;
	stringstream		_IterationMessages;
	bool			_PrintIntermediateResults;
    	int			_ReportEverySteps;
    	EnergyFunction_sp	_EnergyFunction;
	int			_Iteration;
	core::PosixTime_sp		_StartTime;
	bool			_ShowElapsedTime;
	double			_MinGradientMean;
	double			_RMSForce;
	NVector_sp		nvP1DSearchTemp1;
	NVector_sp		nvP1DSearchTemp2;
	NVector_sp		nvP1DSearchOrigin;
	NVector_sp		nvP1DSearchDirection;

    private:	// Do not serialize
//	core::LispCallback_sp		_StepCallback;

    public:
	void lineSearchInitialReport( StepReport_sp report,
				      NVector_sp nvPos, NVector_sp nvDir,
				      NVector_sp nvForce,
				      double xa, double xx, double xb,
				      double fa, double fx, double fb);
	void lineSearchFinalReport( StepReport_sp report, double step, double fmin,
				    int energyEvals, int forceEvals, int dbrentSteps );
	void stepReport( StepReport_sp report, double energy,NVector_sp force );

	void	getPosition(NVector_sp nvResult, NVector_sp nvOrigin, NVector_sp nvDirection, double x);
	double dTotalEnergy(NVector_sp pos);
	double dTotalEnergyForce( NVector_sp nvPos, NVector_sp nvForce);
	double d1DTotalEnergy(double x);
	double d1DTotalEnergyForce( double x, double* fx, double* dfx );
	void   minBracket(	NVector_sp	nvOrigin,
				NVector_sp	nvDir,
				double		*dPxa,
				double		*dPxb,
				double		*dPxc,
				double		*dPfa,
				double		*dPfb,
				double		*dPfc );
	double 	dbrent(	double ax, double bx, double cx,
			double tol, double& step,
			int&	energyEvals,
			int&	forceEvals,
			int&	dbrentSteps );
	void	lineSearch(	double	*dPxnew,
				double	*dPfnew,
				NVector_sp nvOrigin,
				NVector_sp nvDirection,
				NVector_sp nvForce,
				NVector_sp nvTemp1,
				NVector_sp nvTemp2,
	       			int	iteration,
				StepReport_sp report );
	void		define1DSearch( NVector_sp o, NVector_sp d,
					NVector_sp t1, NVector_sp t2 ) {
	    this->nvP1DSearchOrigin = o;
	    this->nvP1DSearchDirection = d;
	    this->nvP1DSearchTemp1 = t1;
	    this->nvP1DSearchTemp2 = t2; };

	void	_steepestDescent( int numSteps,
				  NVector_sp p,
				  double rmsGradientTol );
	void	_conjugateGradient( int numSteps,
				    NVector_sp p,
				    double rmsGradientTol );
	void	_truncatedNewton( int numSteps,
				  NVector_sp p,
				  double rmsGradientTol );


	void	_evaluateEnergyAndForceManyTimes( int numSteps, NVector_sp nvPos );
	void	validateForce(NVector_sp pos, NVector_sp force);

	void	debugBeginIteration(int i);
//	void	debugAdd(addon::QDomNode_sp	node);
	void	debugAccumulateMessage( const char* str );
	void	debugEndIteration();


	void	_displayIntermediateMessage( double		step,
					     double		fnew,
					     double		forceMag,
					     double		forceRmsMag,
					     double		cosAngle,
					     bool		steepestDescent );

    public:	// should be private
	void _truncatedNewtonInnerLoop(
	    int				k,
	    NVector_sp			xj,
	    SparseLargeSquareMatrix_sp	mprecon,
	    SparseLargeSquareMatrix_sp	ldlt,
	    NVector_sp 			force,
	    double				rmsForceMag,
	    NVector_sp			pj,
	    NVector_sp			pjNext,
	    NVector_sp			rj,
	    NVector_sp			dj,
	    NVector_sp			zj,
	    NVector_sp			qj );




    public:


	void throwMinimizerExceededMaxSteps();
	void throwMinimizerStuck();
	void throwMinimizerError();

	EnergyFunction_sp	getEnergyFunction() {_OF(); ASSERTNOTNULL(this->_EnergyFunction); return this->_EnergyFunction; };
	void	useDefaultSettings();

	void	changeOptions(core::Cons_sp options);

	MinimizerStatus status()	{return this->_Status;};
	int statusAsInt()	{return (int)(this->_Status);};
	string		statusMessage() { return this->_Message.str();};

	int	getIteration()	{ return this->_Iteration; };
	string	statusAsString();
	string	statusAsShortString();

	double	getRMSForce() { return this->_RMSForce; };
#ifdef XML_ARCHIVE
	void	debugStart();
	void	debugStop(const string& fileName);
#endif

//	void	setStepCallback(core::LispCallback_sp stepCallback);

	void	setMaximumNumberOfSteepestDescentSteps(int m) {this->_NumberOfSteepestDescentSteps = m;};
	void	setMaximumNumberOfConjugateGradientSteps(int m) {this->_NumberOfConjugateGradientSteps= m;};
	void	setMaximumNumberOfTruncatedNewtonSteps(int m) {this->_NumberOfTruncatedNewtonSteps= m;};
	void	setSteepestDescentTolerance(double m) {this->_SteepestDescentTolerance = m;};
	void	setConjugateGradientTolerance(double m) {this->_ConjugateGradientTolerance = m;};
	void	setTruncatedNewtonTolerance(double m) {this->_TruncatedNewtonTolerance = m;};
	void	setTruncatedNewtonPreconditioner(PreconditionerType m) {this->_TruncatedNewtonPreconditioner = m;};


	void	setEnergyFunction(EnergyFunction_sp ef);
	void	enablePrintIntermediateResults();
	void	disablePrintIntermediateResults();

	string	configurationAsString();

	void	restart();

	void	minimizeSteepestDescent();
	void	minimizeConjugateGradient();
	void	resetAndMinimize();
	void	minimize();
	void	evaluateEnergyAndForceManyTimes(int times);

	addon::QDomNode_sp	asXml();

	DEFAULT_CTOR_DTOR(Minimizer_O);
    };

    extern core::Symbol_sp _sym__PLUS_minimizerStatusConverter_PLUS_;
    extern core::Symbol_sp _sym__PLUS_preconditionerTypeConverter_PLUS_;

};


TRANSLATE(chem::Minimizer_O);



DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::MinimizerStatus,chem::_sym__PLUS_minimizerStatusConverter_PLUS_);
DECLARE_ENUM_SYMBOL_TRANSLATOR(chem::PreconditionerType,chem::_sym__PLUS_preconditionerTypeConverter_PLUS_);


#endif
