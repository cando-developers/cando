       


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














};



