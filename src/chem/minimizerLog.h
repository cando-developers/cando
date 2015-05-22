       
       
//
// (C) 2004 Christian E. Schafmeister
//


#ifndef MinimizerLog_H
#define	MinimizerLog_H

#include "core/common.h"
#include "numericalFunction.h"
#include "nVector.h"
#include "core/holder.h"

#include <ctime>

namespace       chem {



class	ReportBase_O;
    typedef	gctools::smart_ptr<ReportBase_O>	ReportBase_sp;
SMART(ReportBase);
class ReportBase_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ReportBase_O,"ReportBase");

public:
void initialize();
//void	archiveBase(core::ArchiveP node);
public:
    	int	_Time;
public:
	DEFAULT_CTOR_DTOR(ReportBase_O);
};


SMART(MessageReport);
class MessageReport_O : public ReportBase_O
{
    LISP_BASE1(ReportBase_O);
    LISP_CLASS(chem,ChemPkg,MessageReport_O,"MessageReport");

public:
//    void	archive(core::ArchiveP node);
public:
	string	_Message;
public:
	DEFAULT_CTOR_DTOR(MessageReport_O);
};



class	ForceMatchReport_O;
    typedef	gctools::smart_ptr<ForceMatchReport_O>	ForceMatchReport_sp;
SMART(ForceMatchReport);
class ForceMatchReport_O : public ReportBase_O
{
    LISP_BASE1(ReportBase_O);
    LISP_CLASS(chem,ChemPkg,ForceMatchReport_O,"ForceMatchReport");

public:
//    void	archive(core::ArchiveP node);
public:
	string	_Message;
	NVector_sp	_AnalyticalForce;
	NVector_sp	_NumericalForce;
public:
	DEFAULT_CTOR_DTOR(ForceMatchReport_O);
};






SMART(StepReport);
class StepReport_O : public ReportBase_O
{
    LISP_BASE1(ReportBase_O);
    LISP_CLASS(chem,ChemPkg,StepReport_O,"StepReport");

public:
	void	initialize();
//	void	archive(core::ArchiveP node);
public:
	int	_Iteration;
	double	_Xa;
	double	_Xb;
	double	_Xc;
	double	_Fa;
	double	_Fb;
	double	_Fc;
	int	_MinBracketSteps;
	string	_EnergyTermsEnabled;
	double	_TotalEnergy;
	double	_DirectionMagnitude;
	double	_ForceMagnitude;
	string	_MinimizerStatus;
	string	_Direction;
	NVector_sp	_LineSearchPosition;
	NVector_sp	_LineSearchDirection;
	NVector_sp	_LineSearchForce;
	NumericalFunction_sp _TotalEnergyFn;
	NumericalFunction_sp _StretchEnergyFn;
	NumericalFunction_sp _AngleEnergyFn;
	NumericalFunction_sp _DihedralEnergyFn;
	NumericalFunction_sp _NonbondEnergyFn;
	NumericalFunction_sp _ImproperEnergyFn;
	NumericalFunction_sp _ChiralRestraintEnergyFn;
	NumericalFunction_sp _AnchorRestraintEnergyFn;
	NumericalFunction_sp _ImproperRestraintEnergyFn;
	NumericalFunction_sp _FixedNonbondRestraintEnergyFn;
		// Final line search report stuff
	int	_DbrentSteps;
	int	_EnergyEvals;
	int	_ForceEvals;
	double	_Step;
	double	_FMin;
		// Step info
	double				_ForceMag;
	double				_StepEnergy;
	string				_IterationMessages;
	double	_AngleBetweenDirectionAndForceDeg;
		// Termination info
	bool	_PrematureTermination;
	string	_PrematureTerminationMessage;
public:
	void	prematureTermination(const string& msg)
	{
	    this->_PrematureTermination = true;
	    this->_PrematureTerminationMessage = msg;
	};
public:

	DEFAULT_CTOR_DTOR(StepReport_O);
};



    SMART(Minimizer);
SMART(EnergyFunction);



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
SMART(MinimizerLog);
class MinimizerLog_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,MinimizerLog_O,"MinimizerLog");

public:
//    void	archive(core::ArchiveP node);
public:
		// Permanent information
	Minimizer_sp		_Minimizer;
    gctools::Vec0<ReportBase_sp>	_Reports;
public:

    gctools::Vec0<ReportBase_sp>::iterator	begin_Reports() { return this->_Reports.begin(); };
    gctools::Vec0<ReportBase_sp>::iterator	end_Reports() { return this->_Reports.end(); };

	int numberOfReports() { return this->_Reports.size(); }

	void	addReport(ReportBase_sp msg) { this->_Reports.push_back(msg); };
	void	addMessage(const string& str)
	{
	    MessageReport_sp rep = MessageReport_O::create();
	    rep->_Message = str;
	    this->_Reports.push_back(rep);
	}

	DEFAULT_CTOR_DTOR(MinimizerLog_O);
};




};


TRANSLATE(chem::ReportBase_O);
TRANSLATE(chem::MessageReport_O);
TRANSLATE(chem::ForceMatchReport_O);
TRANSLATE(chem::StepReport_O);
TRANSLATE(chem::MinimizerLog_O);
#endif
