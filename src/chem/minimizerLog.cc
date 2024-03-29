/*
    File: minimizerLog.cc
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


#include <cando/chem/minimizerLog.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/nVector.h>
#include <cando/chem/minimizer.h>
#include <clasp/core/wrappers.h>


namespace chem {


    void ReportBase_O::initialize()
    {
	this->Base::initialize();
	this->_Time = (unsigned long)std::time(0);
    }

#ifdef XML_ARCHIVE
    void	ReportBase_O::archiveBase(core::ArchiveP node)
{
    node->attribute("time",this->_Time);
}
#endif

#ifdef XML_ARCHIVE
    void	MessageReport_O::archive(core::ArchiveP node)
{
    this->archiveBase(node);
    node->archiveString("message",this->_Message);
}
#endif

#ifdef XML_ARCHIVE
    void	ForceMatchReport_O::archive(core::ArchiveP node)
{
    this->archiveBase(node);
    node->archiveString("message",this->_Message);
    node->attribute("analyticalForce",this->_AnalyticalForce);
    node->attribute("numericalForce",this->_NumericalForce);
}
#endif


void	StepReport_O::initialize()
{
    this->Base::initialize();
    this->_PrematureTermination = false;
    this->_LineSearchPosition = nil<NVector_O>();
    this->_LineSearchDirection = nil<NVector_O>();
    this->_LineSearchForce = nil<NVector_O>();
    this->_TotalEnergyFn = nil<NumericalFunction_O>();
    this->_StretchEnergyFn = nil<NumericalFunction_O>();
    this->_AngleEnergyFn = nil<NumericalFunction_O>();
    this->_DihedralEnergyFn = nil<NumericalFunction_O>();
    this->_NonbondEnergyFn = nil<NumericalFunction_O>();
    this->_ImproperEnergyFn = nil<NumericalFunction_O>();
    this->_ChiralRestraintEnergyFn = nil<NumericalFunction_O>();
    this->_AnchorRestraintEnergyFn = nil<NumericalFunction_O>();
    this->_PointToLineRestraintEnergyFn = nil<NumericalFunction_O>();
    this->_OutOfZPlaneEnergyFn = nil<NumericalFunction_O>();
    this->_ImproperRestraintEnergyFn = nil<NumericalFunction_O>();
}


#ifdef XML_ARCHIVE
    void	StepReport_O::archive(core::ArchiveP node)
{
    this->archiveBase(node);
    node->attribute("_Iteration",this->_Iteration);
    LOG("Status" );
    node->attribute("_Iteration", this->_Iteration);
    LOG("Status" );
    node->attribute("_Xa", this->_Xa);
    node->attribute("_Xb", this->_Xb);
    node->attribute("_Xc", this->_Xc);
    node->attribute("_Fa", this->_Fa);
    node->attribute("_Fb", this->_Fb);
    node->attribute("_Fc", this->_Fc);
    node->attribute("_MinBracketSteps", this->_MinBracketSteps);
    node->attribute("_EnergyTermsEnabled", this->_EnergyTermsEnabled);
    node->attribute("_TotalEnergy", this->_TotalEnergy);
    node->attribute("_DirectionMagnitude", this->_DirectionMagnitude);
    node->attribute("_ForceMagnitude", this->_ForceMagnitude);
    LOG("Status" );
    node->attribute("_MinimizerStatus", this->_MinimizerStatus);
    LOG("Status" );
    node->attribute("_Direction", this->_Direction);
    LOG("Finished attributes A" );
    LOG("Status" );
    node->attributeIfNotNil("_LineSearchPosition", this->_LineSearchPosition );
#if 0
    ANN(this->_LineSearchPosition);
    if ( this->_LineSearchPosition.notnilp() )
    {
	LOG("this->_LineSearchPosition size={}" , this->_LineSearchPosition->size()  );
    }
#endif
    LOG("Status" );
    node->attributeIfNotNil("_LineSearchDirection", this->_LineSearchDirection);
    LOG("Status" );
    ASSERTNOTNULL(this->_LineSearchForce);
    node->attributeIfNotNil("_LineSearchForce", this->_LineSearchForce);
    LOG("Status" );
    LOG("Finished NVectors A" );
    ASSERTNOTNULL(this->_TotalEnergyFn);
    node->attributeIfNotNil("_TotalEnergyFn",  this->_TotalEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_StretchEnergyFn);
    node->attributeIfNotNil("_StretchEnergyFn",  this->_StretchEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_AngleEnergyFn);
    node->attributeIfNotNil("_AngleEnergyFn",  this->_AngleEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_DihedralEnergyFn);
    node->attributeIfNotNil("_DihedralEnergyFn",  this->_DihedralEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_NonbondEnergyFn);
    node->attributeIfNotNil("_NonbondEnergyFn",  this->_NonbondEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_ImproperEnergyFn);
    node->attributeIfNotNil("_ImproperEnergyFn",  this->_ImproperEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_ChiralRestraintEnergyFn);
    node->attributeIfNotNil("_ChiralRestraintEnergyFn",  this->_ChiralRestraintEnergyFn);
    LOG("Status" );
    ASSERTNOTNULL(this->_AnchorRestraintEnergyFn);
    node->attributeIfNotNil("_AnchorRestraintEnergyFn",  this->_AnchorRestraintEnergyFn);
    ASSERTNOTNULL(this->_PointToLineRestraintEnergyFn);
    node->attributeIfNotNil("_PointToLineRestraintEnergyFn",  this->_PointToLineRestraintEnergyFn);
    ASSERTNOTNULL(this->_OutOfZPlaneEnergyFn);
    node->attributeIfNotNil("_OutOfZPlaneEnergyFn",  this->_OutOfZPlane!EnergyFn);
    ASSERTNOTNULL(this->_ImproperRestraintEnergyFn);
    node->attributeIfNotNil("_ImproperRestraintEnergyFn",  this->_ImproperRestraintEnergyFn);
    LOG("Finished NumFunctions A" );
    node->attribute("_DbrentSteps", this->_DbrentSteps);
    LOG("Status" );
    node->attribute("_EnergyEvals", this->_EnergyEvals);
    LOG("Status" );
    node->attribute("_ForceEvals", this->_ForceEvals);
    node->attribute("_Step", this->_Step);
    node->attribute("_FMin", this->_FMin);
    node->attribute("_ForceMag", this->_ForceMag);
    node->attribute("_StepEnergy", this->_StepEnergy);
    node->attribute("_IterationMessages", this->_IterationMessages);
    node->attribute("_AngleBetweenDirectionAndForceDeg", this->_AngleBetweenDirectionAndForceDeg);
    node->attribute("_PrematureTermination", this->_PrematureTermination);
    node->attribute("_PrematureTerminationMessage", this->_PrematureTerminationMessage);
    LOG("Finished Attributes B" );
    LOG("Status" );
}
#endif



#ifdef XML_ARCHIVE
    void	MinimizerLog_O::archive(core::ArchiveP node)
{
    node->attribute( "_Minimizer",this->_Minimizer );
    node->archiveVector0( "_Reports",this->_Reports );
}
#endif


};
