/*
    File: stereoisomerAtoms.cc
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
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/matter.h>
#include <cando/chem/residue.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <clasp/core/wrappers.h>

namespace chem
{

void StereoisomerAtom_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,AtomName),this->_AtomName);
  node->field(INTERN_(kw,AtomIndex),this->_AtomIndex);
  node->field_if_not_default(INTERN_(kw,Charge),this->_Charge,0.0);
  node->field_if_not_nil(INTERN_(kw,AtomType),this->_AtomType);
};

  
StereoisomerAtom_sp StereoisomerAtom_O::create(ConstitutionAtom_sp constitutionAtom )
{
  if ( constitutionAtom.isA<ConstitutionVirtualAtom_O>() )
  {
    ConstitutionVirtualAtom_sp cva = constitutionAtom.as<ConstitutionVirtualAtom_O>();
    StereoisomerVirtualAtom_sp tva = StereoisomerVirtualAtom_O::create(cva);
    return tva;
  }
  GC_ALLOCATE_VARIADIC(StereoisomerAtom_O, ta, constitutionAtom->_AtomName,constitutionAtom->_Index,0.0,_Nil<core::Symbol_O>());
  return ta;
}

string StereoisomerAtom_O::__repr__() const
{
  stringstream ss;
  ss << "( " << this->className() << " ";
  ss << " :atomName \"" << this->_AtomName << "\"";
  ss << " :atomIndex " << this->_AtomIndex << " )";
  return ss.str();
}






;



StereoisomerVirtualAtom_sp StereoisomerVirtualAtom_O::create(ConstitutionVirtualAtom_sp atom )
{
  StereoisomerVirtualAtom_sp ta = StereoisomerVirtualAtom_O::create();
  ta->_AtomName = atom->_AtomName;
  return ta;
}



void StereoisomerAtoms_O::fields(core::Record_sp node) {
  node->field_if_not_empty(INTERN_(kw,Atoms),this->_Atoms);
};

;

StereoisomerAtoms_sp StereoisomerAtoms_O::create(ConstitutionAtoms_sp constitutionAtoms)
{
  GC_ALLOCATE(StereoisomerAtoms_O, topologyAtoms );
  for ( ConstitutionAtoms_O::iterator it=constitutionAtoms->_Atoms.begin();
        it!=constitutionAtoms->_Atoms.end(); it++)
  {
    ConstitutionAtom_sp ca = *it;
    StereoisomerAtom_sp ta = StereoisomerAtom_O::create(ca);
    topologyAtoms->_Atoms.push_back(ta);
  }
  return topologyAtoms;
}

CL_DEFUN StereoisomerAtoms_sp chem__make_stereoisomer_atoms(ConstitutionAtoms_sp constitutionAtoms)
{
  return StereoisomerAtoms_O::create(constitutionAtoms);
}


void StereoisomerAtoms_O::addStereoisomerVirtualAtom(StereoisomerVirtualAtom_sp ca)
{_OF();
  ConstitutionAtomIndex0N nextIndex = this->_Atoms.size();
  ca->_AtomIndex = nextIndex;
  this->_Atoms.push_back(ca);
};


StereoisomerAtom_sp& StereoisomerAtoms_O::operator[](ConstitutionAtomIndex0N idx)
{_OF();
  ASSERTF(idx>=0&&idx<(int)(this->_Atoms.size()),BF("index[%d] out of range") % idx );
  StereoisomerAtom_sp& ca = this->_Atoms[(uint)idx];
  ASSERTF(idx==ca->_AtomIndex,BF("A mismatch has occured between the index[%d] in the StereoisomerAtoms array and the StereoisomerAtom index[%d] - these have to match or we can't quickly look up atoms by their index") % idx % ca->_AtomIndex );
  return ca;
}

CL_LISPIFY_NAME("atomNamesAsSymbolSet");
CL_DEFMETHOD     adapt::SymbolSet_sp StereoisomerAtoms_O::atomNamesAsSymbolSet()
{_OF();
  adapt::SymbolSet_sp ss = adapt::SymbolSet_O::create();
  for ( gctools::Vec0<StereoisomerAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    ss->insert((*ci)->_AtomName);
  }
  return ss;
};


CL_LISPIFY_NAME("atomWithName");
CL_DEFMETHOD     StereoisomerAtom_sp StereoisomerAtoms_O::atomWithName(MatterName nm)
{_OF();
  for ( gctools::Vec0<StereoisomerAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    if ( (*ci)->_AtomName == nm ) return *ci;
  }
  SIMPLE_ERROR(BF("Could not find StereoisomerAtom with name[%s]") % nm );
}


CL_LISPIFY_NAME("atomWithId");
CL_DEFMETHOD     StereoisomerAtom_sp StereoisomerAtoms_O::atomWithId(ConstitutionAtomIndex0N id) const
{
  return this->_Atoms[id];
}

CL_LISPIFY_NAME("stereoisomer-atoms-index");
CL_DEFMETHOD     int StereoisomerAtoms_O::index(MatterName name) const
{_OF();
  for ( int idx=0; idx<(int)this->_Atoms.size(); idx++ )
  {
    const StereoisomerAtom_sp& atom = this->_Atoms[idx];
    if ( atom->_AtomName == name ) return idx;
  }
  SIMPLE_ERROR(BF("Unknown atom[%s]") % name );
}


};
