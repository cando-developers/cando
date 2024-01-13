/*
    File: pdbMonomer.h
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
#ifndef	PdbMonomer_H //[
#define PdbMonomer_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.fwd.h>

#include <cando/chem/chemPackage.h>


namespace chem {

SMART(Aggregate);
SMART(Matter);

struct ConnectMonomerRec;


    typedef std::pair<MatterName,MatterName> AtomPair;

SMART(PdbMonomerConnectivity );
class PdbMonomerConnectivity_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,PdbMonomerConnectivity_O,"PdbMonomerConnectivity",core::CxxObject_O);
public:
    void initialize();
private:
    core::Symbol_sp           _PdbName;
    core::ComplexVector_T_sp     _LongNames;
    gctools::Vec0<AtomPair> _BondAtomPairs;
public:
    typedef gctools::Vec0<AtomPair>::iterator iterator;
public:
    core::Symbol_sp getPdbName() { return this->_PdbName;};
    void addLongName(MatterName longName);
    void addConnect(MatterName mainAtom, adapt::SymbolSet_sp neighborNames );
    void setPdbName(MatterName pdbName);
    bool hasConnection(MatterName name1, MatterName name2);
    uint connectAtomsAndReturnUnconnected(Residue_sp res, VectorAtom& unconnected);

    iterator begin_Connections() { return this->_BondAtomPairs.begin();};
    iterator end_Connections() {return this->_BondAtomPairs.end();};


    DEFAULT_CTOR_DTOR(PdbMonomerConnectivity_O);
};


    SMART(PdbMonomerDatabase);
    class PdbMonomerDatabase_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,PdbMonomerDatabase_O,"PdbMonomerDatabase",core::CxxObject_O);
public:

    private:
      gc::SmallMap<core::Symbol_sp, PdbMonomerConnectivity_sp>   _Database; // was StringMap
public:
void addMonomer(PdbMonomerConnectivity_sp);
        PdbMonomerConnectivity_sp findMonomer(core::Symbol_sp pdbName);
uint connectVdwOverlappingUnconnectedAtoms(VectorAtom& atoms);

void connectAtoms(Matter_sp m);

DEFAULT_CTOR_DTOR(PdbMonomerDatabase_O);
    };


    };

#endif
