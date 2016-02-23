/*
    File: jobHistory.cc
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
