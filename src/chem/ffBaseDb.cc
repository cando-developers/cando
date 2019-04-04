/*
    File: ffBaseDb.cc
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


//
// (C) 2004 Christian E. Schafmeister
//


/*!
 *	ffBase.cc
 *
 *	Maintain the base class for parameters
 */
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/forceField.h>
//#include "core/archiveNode.h"
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>

SYMBOL_EXPORT_SC_(ChemKwPkg,parameterized);
SYMBOL_EXPORT_SC_(ChemKwPkg,estimated);
SYMBOL_EXPORT_SC_(ChemKwPkg,rough);
SYMBOL_EXPORT_SC_(ChemKwPkg,unknown);

namespace chem {

CL_DEFMETHOD string	FFParameter_O::levelDescription()
{
  switch ( this->_Level ) {
  case parameterized:
      return "parameterized ";
  case estimated:
      return "estimated ";
  case rough:
      return "guessed ";
  case unknown:
      return "unknown ";
  }
  return "unknownParameterLevel ";
}


    core::NullTerminatedEnumAssociation paramLevelEnum[] = {
	{ "parameterized", parameterized },
	{ "estimated", estimated },
	{ "rough", rough },
	{ "", -1 }
};

void FFParameter_O::fields(core::Record_sp node)
{
  //this->Base::fields(node); // T_O
  node->field( INTERN_(kw,level), this->_Level);
}


void FFBaseDb_O::initialize()
{
      this->Base::initialize();
}
void FFBaseDb_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
}

void FFBaseDb_O::forceFieldMerge(FFBaseDb_sp other) {
  // inherited classes may overload this function but they need to call this
  // do nothing
}



void FFParameterBaseDb_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,parms),this->_Parameters);
  this->Base::fields(node);
}

void FFParameterBaseDb_O::initialize()
{
      this->_Parameters = core::HashTableEq_O::create_default();
      this->Base::initialize();
}


void FFParameterBaseDb_O::forceFieldMerge(FFBaseDb_sp other) {
  // inherited classes may overload this function but they need to call this
  // EXCEPT for FFNonbondDb_O - it uses _Parameters differently and we don't
  //        want to mess up what it does by stomping on this->_Parameters the way we
  //        are here.  FFNonbondDb_O will skip over this forceFieldMerge method and
  //        try to call this->Base::Base::forceFieldMerge(other);
  FFParameterBaseDb_sp parm_other = gc::As<FFParameterBaseDb_sp>(other);
  parm_other->_Parameters->maphash([this] (core::T_sp key, core::T_sp parm) {
      this->_Parameters->setf_gethash(key,parm);
    } );
  this->Base::forceFieldMerge(other);
}



string FFParameterBaseDb_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " ";
#ifdef USE_BOEHM
  ss << "@" << (void*)this;
#endif
  ss << ">";
  return ss.str();
}














};
