#define	DEBUG_LEVEL_FULL


#include <cando/chem/minimizerLog.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/nVector.h>
#include <cando/chem/minimizer.h>
#include <clasp/core/wrappers.h>


namespace chem {

    REGISTER_CLASS(chem,ReportBase_O);
    REGISTER_CLASS(chem,StepReport_O);
    REGISTER_CLASS(chem,MessageReport_O);
    REGISTER_CLASS(chem,ForceMatchReport_O);
    REGISTER_CLASS(chem,MinimizerLog_O);

    void ReportBase_O::initialize()
    {_OF();
	this->Base::initialize();
	this->_Time = (unsigned long)std::time(0);
    }

#ifdef XML_ARCHIVE
    void	ReportBase_O::archiveBase(core::ArchiveP node)
{_G();
    node->attribute("time",this->_Time);
}
#endif

#ifdef XML_ARCHIVE
    void	MessageReport_O::archive(core::ArchiveP node)
{_G();
    this->archiveBase(node);
    node->archiveString("message",this->_Message);
}
#endif

#ifdef XML_ARCHIVE
    void	ForceMatchReport_O::archive(core::ArchiveP node)
{_G();
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
    this->_LineSearchPosition = _Nil<NVector_O>();
    this->_LineSearchDirection = _Nil<NVector_O>();
    this->_LineSearchForce = _Nil<NVector_O>();
    this->_TotalEnergyFn = _Nil<NumericalFunction_O>();
    this->_StretchEnergyFn = _Nil<NumericalFunction_O>();
    this->_AngleEnergyFn = _Nil<NumericalFunction_O>();
    this->_DihedralEnergyFn = _Nil<NumericalFunction_O>();
    this->_NonbondEnergyFn = _Nil<NumericalFunction_O>();
    this->_ImproperEnergyFn = _Nil<NumericalFunction_O>();
    this->_ChiralRestraintEnergyFn = _Nil<NumericalFunction_O>();
    this->_AnchorRestraintEnergyFn = _Nil<NumericalFunction_O>();
    this->_ImproperRestraintEnergyFn = _Nil<NumericalFunction_O>();
}


#ifdef XML_ARCHIVE
    void	StepReport_O::archive(core::ArchiveP node)
{_G();
    this->archiveBase(node);
    node->attribute("_Iteration",this->_Iteration);
    LOG(BF("Status") );
    node->attribute("_Iteration", this->_Iteration);
    LOG(BF("Status") );
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
    LOG(BF("Status") );
    node->attribute("_MinimizerStatus", this->_MinimizerStatus);
    LOG(BF("Status") );
    node->attribute("_Direction", this->_Direction);
    LOG(BF("Finished attributes A") );
    LOG(BF("Status") );
    node->attributeIfNotNil("_LineSearchPosition", this->_LineSearchPosition );
#if 0
    ANN(this->_LineSearchPosition);
    if ( this->_LineSearchPosition.notnilp() )
    {
	LOG(BF("this->_LineSearchPosition size=%d") % this->_LineSearchPosition->size()  );
    }
#endif
    LOG(BF("Status") );
    node->attributeIfNotNil("_LineSearchDirection", this->_LineSearchDirection);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_LineSearchForce);
    node->attributeIfNotNil("_LineSearchForce", this->_LineSearchForce);
    LOG(BF("Status") );
    LOG(BF("Finished NVectors A") );
    ASSERTNOTNULL(this->_TotalEnergyFn);
    node->attributeIfNotNil("_TotalEnergyFn",  this->_TotalEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_StretchEnergyFn);
    node->attributeIfNotNil("_StretchEnergyFn",  this->_StretchEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_AngleEnergyFn);
    node->attributeIfNotNil("_AngleEnergyFn",  this->_AngleEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_DihedralEnergyFn);
    node->attributeIfNotNil("_DihedralEnergyFn",  this->_DihedralEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_NonbondEnergyFn);
    node->attributeIfNotNil("_NonbondEnergyFn",  this->_NonbondEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_ImproperEnergyFn);
    node->attributeIfNotNil("_ImproperEnergyFn",  this->_ImproperEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_ChiralRestraintEnergyFn);
    node->attributeIfNotNil("_ChiralRestraintEnergyFn",  this->_ChiralRestraintEnergyFn);
    LOG(BF("Status") );
    ASSERTNOTNULL(this->_AnchorRestraintEnergyFn);
    node->attributeIfNotNil("_AnchorRestraintEnergyFn",  this->_AnchorRestraintEnergyFn);
    ASSERTNOTNULL(this->_ImproperRestraintEnergyFn);
    node->attributeIfNotNil("_ImproperRestraintEnergyFn",  this->_ImproperRestraintEnergyFn);
    LOG(BF("Finished NumFunctions A") );
    node->attribute("_DbrentSteps", this->_DbrentSteps);
    LOG(BF("Status") );
    node->attribute("_EnergyEvals", this->_EnergyEvals);
    LOG(BF("Status") );
    node->attribute("_ForceEvals", this->_ForceEvals);
    node->attribute("_Step", this->_Step);
    node->attribute("_FMin", this->_FMin);
    node->attribute("_ForceMag", this->_ForceMag);
    node->attribute("_StepEnergy", this->_StepEnergy);
    node->attribute("_IterationMessages", this->_IterationMessages);
    node->attribute("_AngleBetweenDirectionAndForceDeg", this->_AngleBetweenDirectionAndForceDeg);
    node->attribute("_PrematureTermination", this->_PrematureTermination);
    node->attribute("_PrematureTerminationMessage", this->_PrematureTerminationMessage);
    LOG(BF("Finished Attributes B") );
    LOG(BF("Status") );
}
#endif



#ifdef XML_ARCHIVE
    void	MinimizerLog_O::archive(core::ArchiveP node)
{_G();
    node->attribute( "_Minimizer",this->_Minimizer );
    node->archiveVector0( "_Reports",this->_Reports );
}
#endif


};

