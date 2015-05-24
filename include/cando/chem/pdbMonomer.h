#ifndef	PdbMonomer_H //[
#define PdbMonomer_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/stringSet.fwd.h>

#include <cando/chem/chemPackage.h>


namespace chem {

SMART(Aggregate);
SMART(Matter);

struct ConnectMonomerRec;


    typedef std::pair<MatterName,MatterName> AtomPair;

SMART(PdbMonomerConnectivity );
class PdbMonomerConnectivity_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,PdbMonomerConnectivity_O,"PdbMonomerConnectivity");
public:
    void initialize();
//void archiveBase(core::ArchiveP node);
private:
    core::Symbol_sp           _PdbName;
    core::VectorObjectsWithFillPtr_sp     _LongNames;
    gctools::Vec0<AtomPair> _BondAtomPairs;
public:
    typedef gctools::Vec0<AtomPair>::iterator iterator;
public:
    core::Symbol_sp getPdbName() { return this->_PdbName;};
    void addLongName(MatterName longName);
    void addConnect(MatterName mainAtom, core::SymbolSet_sp neighborNames );
    void setPdbName(MatterName pdbName);
    bool hasConnection(MatterName name1, MatterName name2);
    uint connectAtomsAndReturnUnconnected(Residue_sp res, VectorAtom& unconnected);

    iterator begin_Connections() { return this->_BondAtomPairs.begin();};
    iterator end_Connections() {return this->_BondAtomPairs.end();};


    DEFAULT_CTOR_DTOR(PdbMonomerConnectivity_O);
};


    SMART(PdbMonomerDatabase);
    class PdbMonomerDatabase_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,PdbMonomerDatabase_O,"PdbMonomerDatabase");
public:
//    void archiveBase(core::ArchiveP node);

    private:
        core::SymbolMap<PdbMonomerConnectivity_O>   _Database; // was StringMap
public:
void addMonomer(PdbMonomerConnectivity_sp);
        PdbMonomerConnectivity_sp findMonomer(core::Symbol_sp pdbName);
uint connectVdwOverlappingUnconnectedAtoms(VectorAtom& atoms);

void connectAtoms(Matter_sp m);

DEFAULT_CTOR_DTOR(PdbMonomerDatabase_O);
    };


    };

TRANSLATE(chem::PdbMonomerConnectivity_O);
TRANSLATE(chem::PdbMonomerDatabase_O);
#endif
