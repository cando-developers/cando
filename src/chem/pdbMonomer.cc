/*
    File: pdbMonomer.cc
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
#include <clasp/core/lispStream.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/pdb.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/atom.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/residue.h>
#include <cando/adapt/stringList.h>
#include <clasp/core/array.h>
#include <clasp/core/package.h>
#include <cando/chem/bond.h>
#include <cando/chem/pdbMonomer.h>
#include <cando/chem/loop.h>
#include <clasp/core/wrappers.h>


namespace chem {




/*
Parse files that look like this:

.123456789.123456789.123456789.123456789.123456890
RESIDUE   000      9
.123456789.123456789.123456789.123456789.123456890
CONECT      C      3 O    OA   OXT
CONECT      O      1 C
CONECT      OA     2 C    CB
CONECT      CB     4 OA   HB   HBA  HBB
CONECT      OXT    2 C    HXT
CONECT      HB     1 CB
CONECT      HBA    1 CB
CONECT      HBB    1 CB
CONECT      HXT    1 OXT
END
HET    000              9
.123456789.123456789.123456789.123456789.123456890
HETNAM     000 methyl hydrogen carbonate
.123456789.123456789.123456789.123456789.123456890
FORMUL      000    C2 H4 O3
*/


struct	ResidueRec
{
    core::Symbol_sp	_PdbCode;
    int		        _NumAtoms;

    void parse(const string& line)
    {
	this->_PdbCode = chemkw_intern(core::trimWhiteSpace(line.substr(10,3)));
	this->_NumAtoms = atoi(core::trimWhiteSpace(line.substr(17,3)).c_str());
    };
};




struct 	ConnectMonomerRec
{
    core::Symbol_sp	_MainAtom;
    adapt::SymbolSet_sp _Bonded;

    void parse(const string& line)
    {
	vector<string>parts;
	parts = core::split(line," ");
	this->_MainAtom = chemkw_intern(parts[1]);
	this->_Bonded = adapt::SymbolSet_O::create();
	for (uint idx=3; idx<parts.size(); idx++ )
	{
          this->_Bonded->insert(chemkw_intern(parts[idx]));
	}
    }
};

struct	HetNameMonomerRec
{
    core::Symbol_sp	_PdbCode;
    core::Symbol_sp	_LongName;
    void parse(const string& line)
    {
	
	this->_PdbCode = chemkw_intern(core::trimWhiteSpace(line.substr(11,3)));
	this->_LongName = chemkw_intern(core::trimWhiteSpace(line.substr(15,99999)));
    };
};


struct FormulaMonomerRec
{
#define	MAX_FORMULA_TERMS 10
    core::Symbol_sp	_PdbCode;
    string	        _FormulaTerm[MAX_FORMULA_TERMS];
    void parse(const string& line)
    {
	vector<string>	parts;
	parts = core::split(line," ");
	this->_PdbCode = chemkw_intern(parts[1]);
	int eidx=0;
	for ( uint idx=2; idx<parts.size(); idx++ )
	{
	    this->_FormulaTerm[eidx] = parts[idx];
	    ASSERTF(eidx<MAX_FORMULA_TERMS,("Allocate more space for formula terms - need at least %d") , eidx );
	    eidx++;
	}
    }
};







DOCGROUP(cando);
CL_DEFUN PdbMonomerDatabase_sp chem__readPdbMonomerConnectivityDatabase(const string& fileName)
    {
	char buffer[1024];
        std::ifstream myfile( fileName.c_str() );
	PdbMonomerDatabase_sp monomerDatabase = PdbMonomerDatabase_O::create();
	if ( myfile.fail() )
	{
	    SIMPLE_ERROR(("Monomer connectivity database file not found: "+fileName));
	}
	{
	    while ( !myfile.eof() )
	    {
		myfile.getline(buffer,1023);
		string oneLine = buffer;
		LOG("PDB line |%s|" , oneLine.c_str()  );
		string recordName = core::trimWhiteSpace(oneLine.substr(0,7));
		if ( recordName == "RESIDUE" )
		{
		    ResidueRec residue;
		    residue.parse(oneLine);
		    PdbMonomerConnectivity_sp oneMonomer = PdbMonomerConnectivity_O::create();
		    oneMonomer->setPdbName(residue._PdbCode);
		    monomerDatabase->addMonomer(oneMonomer);
		    while ( !myfile.eof() )
		    {
			myfile.getline(buffer,1023);
			string oneLine = buffer;
			LOG("PDB line |%s|" , oneLine.c_str()  );
			string recordName = core::trimWhiteSpace(oneLine.substr(0,7));
			if ( recordName == "END" )
			{
			    break;
			} else if ( recordName == "CONECT")
			{
			    ConnectMonomerRec connect;
			    connect.parse(oneLine);
			    oneMonomer->addConnect(connect._MainAtom,connect._Bonded);
			} else
			{
			    SIMPLE_ERROR(("Illegal record(%s) in monomer connection database") , recordName );
			}
		    }
		} else if ( recordName == "HETNAM" )
		{
		    HetNameMonomerRec hetName;
		    hetName.parse(oneLine);
		    PdbMonomerConnectivity_sp oneMonomer = monomerDatabase->findMonomer(hetName._PdbCode);
		    oneMonomer->addLongName(hetName._LongName);
		}
	    }
	}
	return monomerDatabase;
    }








    void PdbMonomerConnectivity_O::initialize()
    {
	this->_LongNames = core::ComplexVector_T_O::make(8,nil<T_O>(),core::clasp_make_fixnum(0));
    }

#ifdef XML_ARCHIVE
    void PdbMonomerConnectivity_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("pdbName",this->_PdbName);
	node->attribute("longNames",this->_LongNames);
	if ( node->loading() )
	{
	    this->_BondAtomPairs.clear();
	    vector<string> data;
	    node->getDataAsVectorOfStrings(data);
	    for ( uint i = 0; i<data.size(); i = i + 2 )
	    {
		AtomPair one(data[i],data[i+1]);
		this->_BondAtomPairs.push_back(one);
	    }
	} else
	{
	    stringstream sout;
	    for ( vector<AtomPair>::iterator pi=this->_BondAtomPairs.begin(); pi!=this->_BondAtomPairs.end(); pi++ )
	    {
		sout << pi->first << " " << pi->second << std::endl;
	    }
	    node->setCharacters(sout.str());
	}
    }
#endif

uint PdbMonomerConnectivity_O::connectAtomsAndReturnUnconnected(Residue_sp res, gctools::Vec0<Atom_sp>& unconnected )
{
    uint unconCount = 0;
    uint unknown = 0;
    uint bondsMade = 0;
    for ( gctools::Vec0<AtomPair>::iterator pi=this->_BondAtomPairs.begin(); pi!=this->_BondAtomPairs.end(); pi++ )
    {
	MatterName s1 = pi->first;
	MatterName s2 = pi->second;
	Atom_sp a1 = res->contentWithNameOrNil(s1).as<Atom_O>();
	Atom_sp a2 = res->contentWithNameOrNil(s2).as<Atom_O>();
	if ( a1.nilp() && a2.notnilp() )
	{
	    LOG("Could not find atom(%s) - marking atom(%s) as unconnected" , pi->first , pi->second );
	    unconnected.push_back(a2);
	    unconCount++;
	    unknown++;
	}
	if ( a2.nilp() && a1.notnilp() ) 
	{
	    LOG("Could not find atom(%s) - marking atom(%s) as unconnected" , pi->second , pi->first );
	    unconnected.push_back(a1);
	    unconCount++;
	    unknown++;
	}
	if ( a1.nilp() && a2.nilp() )
	{
	    LOG("Could not find either atom(%s) or atom(%s)" , pi->first , pi->second );
	    unknown += 2;
	}
	if ( a1.notnilp() && a2.notnilp() ) 
	{
	    LOG("Connecting atoms %s - %s" , a1->description() , a2->description() )
	    a1->bondToSingle(a2);
	    bondsMade++;
	}
    }
    LOG("Total number of unconnected atoms(%d) and unknown atoms(%d)" , unconCount , unknown );
    return bondsMade;
}




    void PdbMonomerConnectivity_O::setPdbName(core::Symbol_sp pdbName)
{
    this->_PdbName = pdbName;
}


    void PdbMonomerConnectivity_O::addLongName(core::Symbol_sp longName)
    {
      this->_LongNames->vectorPushExtend(longName);
    }

void PdbMonomerConnectivity_O::addConnect(core::Symbol_sp mainAtom, adapt::SymbolSet_sp bonded )
{
    bonded->map( [this,&mainAtom] (core::T_sp key) {
            if ( !this->hasConnection(mainAtom,key) )
            {
                AtomPair two(mainAtom,key);
                this->_BondAtomPairs.push_back(two);
            }
        });
}

bool PdbMonomerConnectivity_O::hasConnection(MatterName a1, MatterName a2 )
{
    for ( iterator it=this->_BondAtomPairs.begin(); it!=this->_BondAtomPairs.end(); it++ )
    {
	if ( (it->first==a1 && it->second==a2) || (it->first==a2 && it->second==a1) )
	{
	    return true;
	}
    }
    return false;
}










#ifdef XML_ARCHIVE
    void PdbMonomerDatabase_O::archiveBase(core::ArchiveP node)
    {
	node->attributeStringMap("db",this->_Database);
    }
#endif

void PdbMonomerDatabase_O::addMonomer(PdbMonomerConnectivity_sp mon)
{
    this->_Database.set(mon->getPdbName(),mon);
}


    PdbMonomerConnectivity_sp PdbMonomerDatabase_O::findMonomer(core::Symbol_sp pdbName)
 {
    if ( this->_Database.contains(pdbName) )
    {
	return nil<PdbMonomerConnectivity_O>();
    }
    return this->_Database.get(pdbName);
}


uint PdbMonomerDatabase_O::connectVdwOverlappingUnconnectedAtoms(gctools::Vec0<Atom_sp>& atoms)
{
    uint connected = 0;
    LOG("There are %d unconnected atoms to connect together" , atoms.size()  );
    for ( gctools::Vec0<Atom_sp>::iterator i1 = atoms.begin(); i1!=atoms.end()-1; i1++ )
    {
	Vector3 v1 = (*i1)->getPosition();
//	double rad1 = vdwRadiusForElement((*i1)->getElement());
	for ( gctools::Vec0<Atom_sp>::iterator i2 = i1+1; i2!=atoms.end(); i2++ )
	{
	    Vector3 v2 = (*i2)->getPosition();
	    Vector3 delta = v2.sub(v1);
	    double len = delta.length();
	    double maxLen = maximumBondLengthBetweenElements((*i1)->getElement(),(*i2)->getElement());
	    if ( len < maxLen )
	    {
		if (!(*i1)->isBondedTo(*i2))
		{
		    LOG("Forming bond between atom(%s) and atom(%s) - they were %lf angstroms apart < max(%lf)" , (*i1)->description() , (*i2)->description() , len , maxLen );
		    (*i1)->bondToSingle(*i2);
		    connected++;
		}
	    }
	}
    }
    LOG("Connected %d bonds" , connected );
    return connected;
}







CL_LISPIFY_NAME("connectAtoms");
CL_DEFMETHOD     void PdbMonomerDatabase_O::connectAtoms(Matter_sp matter)
    {
        gctools::Vec0<Atom_sp>	unconnectedAtoms;
	Loop lResidues;
	uint numConnect = 0;
	lResidues.loopTopGoal(matter,RESIDUES);
	while ( lResidues.advance() )
	{
	    Residue_sp res = lResidues.getResidue();
	    if ( this->_Database.contains(res->getPdbName()) )
	    {
		LOG("Found residue(%s) in database connecting atoms" , res->getPdbName()   );
		PdbMonomerConnectivity_sp connector = this->_Database.get(res->getPdbName());
		numConnect += connector->connectAtomsAndReturnUnconnected(res,unconnectedAtoms);
	    } else
	    {
		LOG("Could not find residue(%s) pushing %d atoms into unconnectedAtoms" , res->getPdbName() , res->numberOfAtoms() );
		for ( Matter_O::contentIterator ci=res->begin_contents(); 
		      ci != res->end_contents(); ci++ )
		{
		    unconnectedAtoms.push_back((*ci).as<Atom_O>());
		}
	    }
	}
	uint numVdwConnected = this->connectVdwOverlappingUnconnectedAtoms(unconnectedAtoms);
	LOG("Number of atoms connected using connectivity tables: %d"  , numConnect );
	LOG("Number of atoms connected using vdw overlap: %d"  , numVdwConnected );
        core::writeln_bf_stream(fmt::sprintf("Number of unconnected atoms(%d)" , unconnectedAtoms.size() ));
        core::writeln_bf_stream(fmt::sprintf("Number of atoms connected by monomer connectivity(%d)" , numConnect ));
        core::writeln_bf_stream(fmt::sprintf("Number of atoms connected by vdw overlap(%d)" , numVdwConnected));
    }


   
};
