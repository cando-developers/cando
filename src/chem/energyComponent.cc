#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/chem/energyComponent.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#
// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace chem {








void	EnergyComponent_O::initialize()
{
    this->Base::initialize();
    this->enable();
    this->zeroEnergy();
    this->disableDebug();
    this->setScale(1.0);
    this->setDebug_NumberOfTermsToCalculate(-1);
}


void	EnergyComponent_O::zeroEnergy()
{
    this->_TotalEnergy = 0.0;
}

#ifdef XML_ARCHIVE
    void	EnergyComponent_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("_Enabled",this->_Enabled);
    node->attribute("_Scale",this->_Scale);
    node->attribute("_ErrorThreshold",this->_ErrorThreshold);
    node->attribute("_DebugEnergy",this->_DebugEnergy);
    node->attribute("_TotalEnergy",this->_TotalEnergy);
    node->attribute("_Debug_NumberOfTermsToCalculate",this->_Debug_NumberOfTermsToCalculate);
    node->attribute("_Scale",this->_Scale);
}
#endif


CL_LISPIFY_NAME("debugLogAsString");
CL_DEFMETHOD string EnergyComponent_O::debugLogAsString()
{
#if TURN_ENERGY_FUNCTION_DEBUG_ON == 1
	    return this->_DebugLog.str();
#else
	    return this->className()+"-no_log";
#endif
}

string	EnergyComponent_O::enabledAsString() 
{
stringstream	ss;
    ss.str("");
    if ( this->isEnabled() ) 
    {
	ss<<this->className();
#ifdef	DEBUG_CONTROL_THE_NUMBER_OF_TERMS_EVALAUTED
	if ( this->_Debug_NumberOfTermsToCalculate < 0 ) {
	    ss<<"(debug_All)";
	} else {
	    ss<<"(debug_";
	    ss<<this->_Debug_NumberOfTermsToCalculate<<")";
	}
#endif
	ss<<" ";
    }
    return ss.str();
};

string EnergyComponent_O::summarizeEnergyAsString() 
{
    stringstream ss;
    ss.str("");
    ss<< this->className();
    ss << " energy: " << this->getEnergy() << std::endl;
    return ss.str();
};


};




