#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/chem/energyComponent.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#
// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace chem {


void EnergyComponent_O::exposeCando(core::Lisp_sp e)
{
    core::class_<EnergyComponent_O>()
    .def("enable",&EnergyComponent_O::enable) 
    .def("disable",&EnergyComponent_O::disable) 
	.def("enableDebug",&EnergyComponent_O::enableDebug,"","","",false) 
	.def("disableDebug",&EnergyComponent_O::disableDebug,"","","",false) 
    .def("debugLogAsString",&EnergyComponent_O::debugLogAsString) 
    .def("getEnergy",&EnergyComponent_O::getEnergy) 
    .def("setDebug_NumberOfTermsToCalculate", 
	    &EnergyComponent_O::setDebug_NumberOfTermsToCalculate) 
    .def("setErrorThreshold",&EnergyComponent_O::setErrorThreshold) 
    .def("getErrorThreshold",&EnergyComponent_O::getErrorThreshold) 
    ;
}

    void EnergyComponent_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,EnergyComponent,"","",_lisp)
    ;
#endif //]
}




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


CL_NAME("debugLogAsString");
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


    EXPOSE_CLASS(chem,EnergyComponent_O);
};




