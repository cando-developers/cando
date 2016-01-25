       
#define	DEBUG_LEVEL_NONE


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
    SIMPLE_ERROR(BF("Illegal Chimera atom specification"));
  }
  residueStart = spec.find_first_of(":")+1;
  LOG(BF("residueStart(%d)") % residueStart );
  string residueInfo = spec.substr(residueStart,atomNameStart-residueStart-1);
  LOG(BF("residueInfo(%s)") % residueInfo.c_str() );
  atomName = spec.substr(atomNameStart,9999);
  LOG(BF("atomName(%s)") % atomName.c_str() );
  size_t chainStart = residueInfo.find_first_of(".");
  chain = "";
  if ( chainStart != string::npos )
  {
    chain = residueInfo.substr(chainStart+1,9999);
  } else
  {
    chainStart = residueInfo.size();
  }
  LOG(BF("chain(%s)") % chain.c_str() );
  LOG(BF("chainStart(%d)") % chainStart );
  string seqNumStr = residueInfo.substr(0,chainStart);
  LOG(BF("seqNumStr(%s)") % seqNumStr.c_str() );
  sequenceNumber = atoi(seqNumStr.c_str());
  LOG(BF("sequenceNumber(%d)") % sequenceNumber );
}

/*! Split the string containing any number of Chimera atom specifications
that look like:   #1:1.A@H_289092  into a list of lists that have the form:
((residue-number chain-id atom-name) ...)
The chain-id and atom-name will be CKW keyword symbols.
*/
CL_DEFUN core::List_sp chem__parseChimeraAtomSpecs(const string& spec)
{
  vector<string> specs = core::split(spec," \n\t");
  core::List_sp result(_Nil<core::T_O>());
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
