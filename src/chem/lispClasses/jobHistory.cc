       


#include "jobHistory.h"
//#include "core/archive.h"
//#include "core/archiveNode.h"
#include <ctime>
#include "core/wrappers.h"

namespace chem {

void TrainerHistoryEntry_O::initialize()
{
    this->Base::initialize();
    this->_Command = "";
    this->_Outcome = "";
    this->_Comment = "";
    this->_StartSeconds = 0;
    this->_StopSeconds = 0;

}


#ifdef XML_ARCHIVE
    void	TrainerHistoryEntry_O::archive(core::ArchiveP node)
{
    node->attribute("Command",this->_Command);
    node->attribute("StartSeconds",this->_StartSeconds);
    node->attribute("StopSeconds",this->_StopSeconds);
    node->attribute("Outcome",this->_Outcome);
    if ( node->loading() )
    {
        this->_Comment = node->characters();
    } else
    {
        node->setCharacters(this->_Comment);
    }
}
#endif


void	TrainerHistoryEntry_O::start(const string& command)
{
std::time_t	t1 = std::time(0);
int		seconds;
    seconds = (unsigned long)(t1);
    this->setStartSeconds(seconds);
    this->setCommand(command);
}


void	TrainerHistoryEntry_O::finish(const string& newStateName)
{
std::time_t	t1 = std::time(0);
int		seconds;
    seconds = (unsigned long)(t1);
    this->setStopSeconds(seconds);
    this->setState(newStateName);
}


#ifdef XML_ARCHIVE
    void	TrainerHistory_O::archive(core::ArchiveP node)
{
    node->archiveVector0<TrainerHistoryEntry_O>("entries",this->_Entries);
}
#endif








    void TrainerHistoryEntry_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<TrainerHistoryEntry_O>()
	.def("getElapsedSeconds",&TrainerHistoryEntry_O::getElapsedSeconds)
	.def("start",&TrainerHistoryEntry_O::start)
	.def("finish",&TrainerHistoryEntry_O::finish)
	.def("setComment",&TrainerHistoryEntry_O::setComment)
    ;

}
    void TrainerHistoryEntry_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,TrainerHistoryEntry,"","",_lisp)
// 	.def_readwrite("_Command",&TrainerHistoryEntry_O::_Command)
// 	.def_readwrite("_StartSeconds",&TrainerHistoryEntry_O::_StartSeconds)
// 	.def_readwrite("_StopSeconds",&TrainerHistoryEntry_O::_StopSeconds)
// 	.def_readwrite("_Outcome",&TrainerHistoryEntry_O::_Outcome)
// 	.def_readwrite("_Comment",&TrainerHistoryEntry_O::_Comment)
// 	.def_readwrite("_State",&TrainerHistoryEntry_O::_State)
	.def("getElapsedSeconds",&TrainerHistoryEntry_O::getElapsedSeconds)
	.def("start",&TrainerHistoryEntry_O::start)
	.def("finish",&TrainerHistoryEntry_O::finish)
	.def("setComment",&TrainerHistoryEntry_O::setComment)
	;
#endif //]
}




    void TrainerHistory_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<TrainerHistory_O>()
	.def("addToHistory", &TrainerHistory_O::add)
    ;

}
    void TrainerHistory_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,TrainerHistory,"","",_lisp)
//	.add_property("iterate_Entries",
//			boost::python::range(&TrainerHistory_O::begin_Entries,
//				&TrainerHistory_O::end_Entries))
	.def("addToHistory", &TrainerHistory_O::add)
	;
#endif //]
}


    EXPOSE_CLASS(chem,TrainerHistoryEntry_O);
    EXPOSE_CLASS(chem,TrainerHistory_O);
};



