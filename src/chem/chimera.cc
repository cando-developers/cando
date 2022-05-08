/*
    File: chimera.cc
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


#include <clasp/core/foundation.h>
#include <clasp/core/common.h>
#include <cando/chem/chimera.h>
#include <string>
#include <clasp/core/wrappers.h>

namespace chem
{
void parseChimeraAtomSpecification(const string& spec, uint& sequenceNumber, string& chain, string& atomName)
{
  uint 	residueStart;
  size_t atomNameStart;
  atomNameStart = spec.find_first_of("@")+1;
  if ( atomNameStart == string::npos ) {
    SIMPLE_ERROR(("Illegal Chimera atom specification"));
  }
  residueStart = spec.find_first_of(":")+1;
  LOG("residueStart(%d)" , residueStart );
  string residueInfo = spec.substr(residueStart,atomNameStart-residueStart-1);
  LOG("residueInfo(%s)" , residueInfo.c_str() );
  atomName = spec.substr(atomNameStart,9999);
  LOG("atomName(%s)" , atomName.c_str() );
  size_t chainStart = residueInfo.find_first_of(".");
  chain = "";
  if ( chainStart != string::npos )
  {
    chain = residueInfo.substr(chainStart+1,9999);
  } else
  {
    chainStart = residueInfo.size();
  }
  LOG("chain(%s)" , chain.c_str() );
  LOG("chainStart(%d)" , chainStart );
  string seqNumStr = residueInfo.substr(0,chainStart);
  LOG("seqNumStr(%s)" , seqNumStr.c_str() );
  sequenceNumber = atoi(seqNumStr.c_str());
  LOG("sequenceNumber(%d)" , sequenceNumber );
}

/*! Split the string containing any number of Chimera atom specifications
that look like:   #1:1.A@H_289092  into a list of lists that have the form:
((residue-number chain-id atom-name) ...)
The chain-id and atom-name will be CKW keyword symbols.
*/
DOCGROUP(cando)
CL_DEFUN core::List_sp chem__parseChimeraAtomSpecs(const string& spec)
{
  vector<string> specs = core::split(spec," \n\t");
  core::List_sp result(nil<core::T_O>());
  uint sequenceNumber;
  string chainStr;
  string atomNameStr;
  for ( int i(specs.size()-1); i>= 0; --i ) {
    parseChimeraAtomSpecification(specs[i],sequenceNumber,chainStr,atomNameStr);
    core::Symbol_sp chain = chemkw_intern(chainStr);
    core::Symbol_sp atomName = chemkw_intern(atomNameStr);
    core::Cons_sp one = core::Cons_O::createList(core::clasp_make_fixnum(sequenceNumber),chain,atomName);
    result = core::Cons_O::create(one,result);
  }
  return result;
}

void initialize_chimera() {
//  core::af_def(ChemPkg,"parse-chimera-atom-specs",&parseChimeraAtomSpecs);
};
};
