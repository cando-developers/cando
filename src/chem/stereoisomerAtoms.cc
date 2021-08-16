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
  node->field(INTERN_(kw,ConstitutionAtomIndex),this->_ConstitutionAtomIndex);
  node->field(INTERN_(kw,AtomCharge),this->_AtomCharge);
  node->field(INTERN_(kw,AtomType),this->_AtomType);
};

  
StereoisomerAtom_sp StereoisomerAtom_O::create(ConstitutionAtom_sp constitutionAtom )
{
  if ( constitutionAtom.isA<ConstitutionVirtualAtom_O>() )
  {
    ConstitutionVirtualAtom_sp cva = constitutionAtom.as<ConstitutionVirtualAtom_O>();
    StereoisomerVirtualAtom_sp tva = StereoisomerVirtualAtom_O::create(cva);
    return tva;
  }
  auto  ta = gctools::GC<StereoisomerAtom_O>::allocate( constitutionAtom->_AtomName,constitutionAtom->_Index,0.0,nil<core::Symbol_O>());
  return ta;
}

StereoisomerAtom_sp StereoisomerAtom_O::make(core::Symbol_sp atomName, core::Symbol_sp type, double charge, ConstitutionAtomIndex0N index )
{
  auto  ta = gctools::GC<StereoisomerAtom_O>::allocate( atomName, index, charge, type );
  return ta;
}

string StereoisomerAtom_O::__repr__() const
{
  stringstream ss;
  ss << "( " << this->className() << " ";
  ss << " :atomName \"" << this->_AtomName << "\"";
  ss << " :constitutionAtomIndex " << this->_ConstitutionAtomIndex << " )";
  return ss.str();
}






;



StereoisomerVirtualAtom_sp StereoisomerVirtualAtom_O::create(ConstitutionVirtualAtom_sp atom )
{
  StereoisomerVirtualAtom_sp ta = StereoisomerVirtualAtom_O::create();
  ta->_AtomName = atom->_AtomName;
  return ta;
}


StereoisomerAtoms_O::StereoisomerAtoms_O(core::Symbol_sp stereoisomerName) : _StereoisomerName(stereoisomerName)
{
  if (stereoisomerName.unboundp()) {
    SIMPLE_ERROR(BF("The stereoisomerName was unbound"));
  }
};

void StereoisomerAtoms_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,StereoisomerName),this->_StereoisomerName);
  node->field_if_not_empty(INTERN_(kw,Atoms),this->_Atoms);
};

;
StereoisomerAtoms_sp StereoisomerAtoms_O::create(core::Symbol_sp stereoisomerName )
{
  auto  stereoisomerAtoms = gctools::GC<StereoisomerAtoms_O>::allocate( stereoisomerName );
  return stereoisomerAtoms;
}

StereoisomerAtoms_sp StereoisomerAtoms_O::create(core::Symbol_sp stereoisomerName, ConstitutionAtoms_sp constitutionAtoms)
{
  auto  stereoisomerAtoms = gctools::GC<StereoisomerAtoms_O>::allocate( stereoisomerName );
  for ( ConstitutionAtoms_O::iterator it=constitutionAtoms->_Atoms.begin();
        it!=constitutionAtoms->_Atoms.end(); it++)
  {
    ConstitutionAtom_sp ca = *it;
    StereoisomerAtom_sp ta = StereoisomerAtom_O::create(ca);
    stereoisomerAtoms->_Atoms.push_back(ta);
  }
  return stereoisomerAtoms;
}

DOCGROUP(cando)
CL_DEFUN StereoisomerAtoms_sp chem__make_stereoisomer_atoms(core::Symbol_sp stereoisomerName, ConstitutionAtoms_sp constitutionAtoms)
{
  return StereoisomerAtoms_O::create(stereoisomerName,constitutionAtoms);
}


CL_DEFMETHOD core::Symbol_sp StereoisomerAtoms_O::getName() const {
  return this->_StereoisomerName;
}

CL_DEFMETHOD void StereoisomerAtoms_O::addStereoisomerAtom(StereoisomerAtom_sp ca)
{_OF();
  ConstitutionAtomIndex0N nextIndex = this->_Atoms.size();
  ca->_ConstitutionAtomIndex = nextIndex;
  this->_Atoms.push_back(ca);
};


StereoisomerAtom_sp& StereoisomerAtoms_O::operator[](ConstitutionAtomIndex0N idx)
{_OF();
  ASSERTF(idx>=0&&idx<(int)(this->_Atoms.size()),BF("index[%d] out of range") % idx );
  StereoisomerAtom_sp& ca = this->_Atoms[(uint)idx];
//  ASSERTF(idx==ca->_AtomIndex,BF("A mismatch has occured between the index[%d] in the StereoisomerAtoms array and the StereoisomerAtom index[%d] - these have to match or we can't quickly look up atoms by their index") % idx % ca->_AtomIndex );
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


CL_LAMBDA(name &optional (errorp t))
CL_LISPIFY_NAME("atomWithName");
CL_DEFMETHOD     core::T_mv StereoisomerAtoms_O::atomWithName(MatterName nm,bool errorp)
{_OF();
  for ( gctools::Vec0<StereoisomerAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    if ( (*ci)->_AtomName == nm ) return Values(*ci,_lisp->_true());
  }
  if (!errorp) {
    return Values(nil<core::T_O>(),nil<core::T_O>());
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
