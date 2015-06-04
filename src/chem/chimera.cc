       
#define	DEBUG_LEVEL_NONE


#include <clasp/core/common.h>
#include <cando/chem/chimera.h>
#include <string>

namespace chem
{
    void	parseChimeraAtomSpecification(const string& spec, uint& sequenceNumber, string& chain, string& atomName)
    {_G();
	uint 	residueStart;
	size_t atomNameStart;
	atomNameStart = spec.find_first_of("@")+1;
	if ( atomNameStart == string::npos )
	{
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

};
