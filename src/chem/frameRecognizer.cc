/*
    File: frameRecognizer.cc
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

#include <clasp/core/common.h>
#include <cando/chem/frameRecognizer.h>
#include <clasp/core/array.h>
#include <clasp/core/environment.h>
#include <clasp/core/wrappers.h>

namespace chem {


void	FrameRecognizer_O::initialize()
{_OF();
    this->_Name = _Nil<core::Symbol_O>();
    this->_GroupName = _Nil<core::Symbol_O>();
}

CL_LISPIFY_NAME("frame-recognizer-compile-smarts");
CL_DEFMETHOD void	FrameRecognizer_O::compileSmarts(const string& osm)
{
    this->_Smarts = osm;
    this->_Root = chem__compile_smarts(this->_Smarts,_Nil<core::T_O>());
    LOG(BF("Successfully compiled smarts code for atom O") );
//    LOG(BF("ChemInfo code = %s") % this->_ChemInfo->asXmlString());
}


#ifdef XML_ARCHIVE
    void	FrameRecognizer_O::archive(core::ArchiveP node)
{_OF();
#if PRODUCTION_CODE   // FIXME use "name" only and remove the test for "_key"
    node->attribute("name",this->_Name);
#else
    if ( node->loading() )
    {
	if ( node->hasAttribute("name") )
	{
	    node->attribute("name",this->_Name);
	} else
	{
	    node->attribute("_key",this->_Name);
	}
    } else
    {
	node->attribute("name",this->_Name);
    }
#endif
    node->attribute("smarts",this->_Smarts);
    node->attributeIfNotNil("groupName",this->_GroupName);
    	// Compile the smarts strings into ChemInfo code
    if ( node->loading() )
    {
        TRY()
	{
            this->compileSmarts(this->_Smarts);
	} catch ( core::Condition& err )
	{
	    SIMPLE_ERROR(BF("Could not parse Smarts code for %s %s")% this->_Smarts % err.message() );
	}
    }
};
#endif



CL_LISPIFY_NAME("recognizes");
CL_DEFMETHOD bool	FrameRecognizer_O::recognizes( Atom_sp o )
{
    ASSERTNOTNULL(this->_ChemInfo);
    ASSERT(this->_ChemInfo->compileSucceeded());
    core::T_mv match_mv = chem__chem_info_match(this->_Root,o);
    return match_mv.notnilp();
}


#if 0
CL_LISPIFY_NAME("getMatch");
CL_DEFMETHOD ChemInfoMatch_sp FrameRecognizer_O::getMatch()
{
  ChemInfoMatch_sp match = this->_ChemInfo->matches_atom(
    return this->_ChemInfo->getMatch();
}
#endif


    void FrameRecognizer_O::setRecognizerName(core::Symbol_sp fn)
{
    this->_Name = fn;
}

#if 0
CL_LISPIFY_NAME("depth");
CL_DEFMETHOD uint FrameRecognizer_O::depth()
{
    ASSERTNOTNULL(this->_ChemInfo);
    return this->_->depth();
}
#endif

CL_LISPIFY_NAME("getRecognizerName");
CL_DEFMETHOD     core::Symbol_sp FrameRecognizer_O::getRecognizerName()
{
    return this->_Name;
}

void FrameRecognizer_O::setGroupName(core::Symbol_sp fn)
{
    this->_GroupName = fn;
}

CL_LISPIFY_NAME("getGroupName");
CL_DEFMETHOD core::Symbol_sp FrameRecognizer_O::getGroupName()
{
    if ( this->_GroupName.nilp() )
	return this->_Name;
    return this->_GroupName;
}



string FrameRecognizer_O::description() const
{
    stringstream ss;
    ss << "( " << this->className() << " ";
    ss << _rep_(this->_Name);
    ss << " :groupName \"" << _rep_(this->_GroupName) << "\"";
    ss << " :smarts " << "\"" << this->_Smarts << "\" )";
    return ss.str();
}


#define ARGS_FrameRecognizer_O_make "(name smarts group_name)"
#define DECL_FrameRecognizer_O_make ""
#define DOCS_FrameRecognizer_O_make "make FrameRecognizer"
  FrameRecognizer_sp FrameRecognizer_O::make(core::Symbol_sp name, const string& smarts, core::Symbol_sp groupName)
  {
      GC_ALLOCATE(FrameRecognizer_O, me );
    me->_Name = name;
    me->_Smarts = smarts;
    me->_GroupName = groupName;
    me->_Root = chem__compile_smarts(me->_Smarts,_Nil<core::T_O>());
    return me;
  };




};
