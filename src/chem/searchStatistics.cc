/*
    File: searchStatistics.cc
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
#define	DEBUG_LEVEL_FULL

#include <cando/chem/searchStatistics.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#
// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace chem {

#if INIT_TO_FACTORIES

#define ARGS_SearchStatistics_O_make "()"
#define DECL_SearchStatistics_O_make ""
#define DOCS_SearchStatistics_O_make "make SearchStatistics"
  SearchStatistics_sp SearchStatistics_O::make()
  {
    IMPLEMENT_ME();
  };

#else

    core::T_sp 	SearchStatistics_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{
    IMPLEMENT_ME();
}

#endif

void	SearchStatistics_O::initialize()
{
    this->Base::initialize();
    this->_ReportFrequency = UndefinedUnsignedInt;
    this->reset();
}

void	SearchStatistics_O::reset()
{
    this->_Builds = 0;
    this->_Rejects = 0;
    this->_Hits = 0;
}



#ifdef XML_ARCHIVE
    void	SearchStatistics_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("ReportFrequency",this->_ReportFrequency);
    node->attribute("Builds",this->_Builds);
    node->attribute("Rejects",this->_Rejects);
    node->attribute("Hits",this->_Hits);
}
#endif


void	SearchStatistics_O::setReportFrequency(uint i)
{
    this->_ReportFrequency = i;
}


CL_LISPIFY_NAME("incrementBuilds");
CL_DEFMETHOD void	SearchStatistics_O::incrementBuilds()
{_OF();
    this->_Builds++;
    if ( this->_ReportFrequency == UndefinedUnsignedInt ) return;
    if ( ( this->_Builds % this->_ReportFrequency ) == 0 )
    {
        _lisp->print(BF("SearchStatistics_O::incrementBuilds>>%s") % this->description().c_str() );
    }
}

CL_LISPIFY_NAME("getBuilds");
CL_DEFMETHOD Bignum SearchStatistics_O::getBuilds()
{
    return this->_Builds;
}


CL_LISPIFY_NAME("incrementRejects");
CL_DEFMETHOD void	SearchStatistics_O::incrementRejects()
{
    this->_Rejects++;
}

CL_LISPIFY_NAME("getRejects");
CL_DEFMETHOD Bignum SearchStatistics_O::getRejects()
{
    return this->_Rejects;
}


CL_LISPIFY_NAME("incrementHits");
CL_DEFMETHOD void	SearchStatistics_O::incrementHits()
{
    this->_Hits++;
}

CL_LISPIFY_NAME("getHits");
CL_DEFMETHOD Bignum SearchStatistics_O::getHits()
{
    return this->_Hits;
}

string	SearchStatistics_O::description() const
{
    boost::format fmt("Builds[%s] Rejects[%s] Hits[%s]");
    fmt % this->_Builds % this->_Rejects % this->_Hits;
    return fmt.str();
}









};
