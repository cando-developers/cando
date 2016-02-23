/*
    File: minimizerLog.h
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


#ifndef MinimizerLog_H
#define	MinimizerLog_H

#include <clasp/core/common.h>
#include <cando/chem/numericalFunction.h>
#include <cando/chem/nVector.h>
#include <clasp/core/holder.h>

#include <ctime>

namespace       chem {



class	ReportBase_O;
    typedef	gctools::smart_ptr<ReportBase_O>	ReportBase_sp;
SMART(ReportBase);
class ReportBase_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ReportBase_O,"ReportBase",core::CxxObject_O);

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
    LISP_CLASS(chem,ChemPkg,MessageReport_O,"MessageReport",ReportBase_O);

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
    LISP_CLASS(chem,ChemPkg,ForceMatchReport_O,"ForceMatchReport",ReportBase_O);

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
    LISP_CLASS(chem,ChemPkg,StepReport_O,"StepReport",ReportBase_O);

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
class MinimizerLog_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,MinimizerLog_O,"MinimizerLog",core::CxxObject_O);

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
