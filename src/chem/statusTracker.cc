/*
    File: statusTracker.cc
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

#include <cando/chem/statusTracker.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>



namespace chem {




//
// Constructor
//

void	StatusTracker_O::initialize()
{
    this->Base::initialize();
    this->_HasError = false;
    this->_ErrorStream.str("");
    this->_MessageStream.str("");
}




    StatusTracker_O::StatusTracker_O(const StatusTracker_O& ss ) : core::CxxObject_O(ss)
{
    this->_HasError = ss._HasError;
    this->_MessageStream.str(ss._MessageStream.str());
    this->_ErrorStream.str(ss._ErrorStream.str());
}




//
// Destructor
//





#ifdef XML_ARCHIVE
    void	StatusTracker_O::archive(core::ArchiveP node)
{
    node->attribute("hasError", this->_HasError );
    node->archiveStringStreamIfDefined("ErrorMessages",this->_ErrorStream);
    node->archiveStringStreamIfDefined("Messages",this->_MessageStream);
}
#endif





void	StatusTracker_O::addMessage(const string& msg)
{
    LOG(BF("Message: %s") % msg.c_str()  );
    this->_MessageStream << msg << std::endl;
}





void	StatusTracker_O::addError(const string& msg)
{
    LOG(BF("Message: %s") % msg.c_str()  );
    this->_HasError = true;
    this->_ErrorStream << msg << std::endl;
}






void	StatusTracker_O::reset()
{
    this->_HasError = false;
    this->_ErrorStream.str("");
    this->_MessageStream.str("");
}

string	StatusTracker_O::getStatus()
{
stringstream	stuff;
    stuff.str("");
    stuff << this->_MessageStream.str();
    if ( this->_HasError )
    {
        stuff << this->_ErrorStream.str();
    }
    LOG(BF("Status = |%s|") % stuff.str().c_str()  );
    return stuff.str();
}


StatusTracker_sp	StatusTracker_O::copy()
{
    GC_COPY(StatusTracker_O, stnew , *this); // = RP_Copy<StatusTracker_O>(this);
    return stnew;
}






};
