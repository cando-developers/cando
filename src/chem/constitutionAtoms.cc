/*
    File: constitutionAtoms.cc
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
#include <clasp/core/bformat.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/matter.h>
#include <cando/chem/residue.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <clasp/core/wrappers.h>

namespace chem
{

void ConstitutionBond_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,toatom),this->_ToAtomIndex);
  node->field(INTERN_(kw,bond_order),this->_BondOrder);
  this->Base::fields(node);
}

void ConstitutionAtom_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,name),this->_AtomName);
  node->field(INTERN_(kw,element),this->_Element);
  node->field(INTERN_(kw,bonds),this->_Bonds);
  this->Base::fields(node);
}  

string ConstitutionAtom_O::__repr__() const
{
  stringstream ss;
  ss << "(" << this->className() << " ";
  ss << " :atomName \"" << this->_AtomName << "\"";
  ss << " :element " << symbolFromElement(this->_Element);
  ss << ")";
  return ss.str();
}




;



#ifdef XML_ARCHIVE
void ConstitutionVirtualAtom_O::archiveBase(core::ArchiveP node)
{_OF();
  this->Base::archiveBase(node);
}
#endif


  

CL_LISPIFY_NAME(makeConstitutionAtomsFromResidue);
CL_LAMBDA(residue &optional verbose);
CL_DEFUN ConstitutionAtoms_sp ConstitutionAtoms_O::makeConstitutionAtomsFromResidue(Residue_sp residue, bool verbose)
{
  residue->ensureAllAtomNamesAreUnique();
  GC_ALLOCATE(ConstitutionAtoms_O, catoms );
  Matter_O::contentIterator it;
  ConstitutionAtomIndex0N index = 0;
  MapAtomsToConstitutionAtomIndex0N atomToIndexMap;
  {_BLOCK_TRACEF(BF("First assign each atom a unique ConstitutionAtomIndex0N"));
    for ( it=residue->begin_atoms(); it!=residue->end_atoms(); it++, index++ )
    {
      Atom_sp atom = (*it).as<Atom_O>();
      ConstitutionAtom_sp catom = atom->asConstitutionAtom(index);
      catoms->_Atoms.push_back(catom);
      atomToIndexMap[atom] = index;
      if (verbose) BFORMAT_T(BF("Atom %s index: %d\n") % _rep_(atom) % index);
    }
  }
  {_BLOCK_TRACE("Then define the bonds for each atom using the indices");
    index = 0; // start indexing from 0 again
    for ( it=residue->begin_atoms(); it!=residue->end_atoms(); it++, index++ )
    {
      Atom_sp atom = (*it).as<Atom_O>();
      ConstitutionAtom_sp catom = catoms->_Atoms[index];
      atom->defineConstitutionAtomBonding(catom,atomToIndexMap);
    }
  }
  return catoms;
}


void ConstitutionAtoms_O::addConstitutionVirtualAtom(ConstitutionVirtualAtom_sp ca)
{_OF();
  ConstitutionAtomIndex0N nextIndex = this->_Atoms.size();
  this->_Atoms.push_back(ca);
};


ConstitutionAtom_sp& ConstitutionAtoms_O::operator[](ConstitutionAtomIndex0N idx)
{_OF();
  ASSERTF(idx>=0&&idx<(int)(this->_Atoms.size()),BF("index[%d] out of range") % idx );
  ConstitutionAtom_sp& ca = this->_Atoms[(uint)idx];
  ASSERTF(idx==ca->_AtomIndex,BF("A mismatch has occured between the index[%d] in the ConstitutionAtoms array and the ConstitutionAtom index[%d] - these have to match or we can't quickly look up atoms by their index") % idx % ca->_AtomIndex );
  return ca;
}

CL_LISPIFY_NAME("atomNamesAsSymbolSet");
CL_DEFMETHOD     adapt::SymbolSet_sp ConstitutionAtoms_O::atomNamesAsSymbolSet()
{_OF();
  adapt::SymbolSet_sp ss = adapt::SymbolSet_O::create();
  for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    ss->insert((*ci)->_AtomName);
  }
  return ss;
};

void ConstitutionAtoms_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,atoms),this->_Atoms);
  this->Base::fields(node);
}  


CL_LISPIFY_NAME("atomWithName");
CL_DEFMETHOD     ConstitutionAtom_sp ConstitutionAtoms_O::atomWithName(MatterName nm)
{_OF();
  for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    if ( (*ci)->_AtomName == nm ) return *ci;
  }
  SIMPLE_ERROR(BF("Could not find ConstitutionAtom with name[%s]") % nm );
}


CL_LISPIFY_NAME("atomWithId");
CL_DEFMETHOD     ConstitutionAtom_sp ConstitutionAtoms_O::atomWithId(ConstitutionAtomIndex0N id) const
{
  return this->_Atoms[id];
}

CL_LISPIFY_NAME("constitution-atoms-index");
CL_DEFMETHOD     int ConstitutionAtoms_O::index(MatterName name) const
{_OF();
  for ( int idx=0; idx<(int)this->_Atoms.size(); idx++ )
  {
    const ConstitutionAtom_sp& atom = this->_Atoms[idx];
    if ( atom->_AtomName == name ) return idx;
  }
  SIMPLE_ERROR(BF("Unknown atom[%s]") % name );
}


#if 0
CL_LISPIFY_NAME("buildResidue");
CL_DEFMETHOD     Residue_sp ConstitutionAtoms_O::buildResidue()
{
  Residue_sp res = Residue_O::create();
  gctools::Vec0<Atom_sp> atoms;
  atoms.resize(this->_Atoms.size());
  res->resizeContents(this->_Atoms.size());
  size_t idx = 0;
  for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ++ci, ++idx )
  {
    Atom_sp atom = Atom_O::create();
    atom->setName((*ci)->_AtomName);
    atom->setElement((*ci)->_Element);
    atoms[idx] = atom;
    res->putMatter((*ci)->_AtomIndex,atom);
  }
  for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    ASSERT_lt((*ci)->_AtomIndex,(int)atoms.size());
    Atom_sp fromAtom = atoms[(*ci)->_AtomIndex];
    for ( gctools::Vec0<ConstitutionBond_sp>::const_iterator bi=(*ci)->_Bonds.begin(); bi!=(*ci)->_Bonds.end(); bi++ )
    {
		// Avoid making bonds twice
      if ( (*bi)->_AtomIndex < (*ci)->_AtomIndex ) continue;
      ASSERT_lt((*bi)->_AtomIndex,(int)atoms.size());
      Atom_sp toAtom = atoms[(*bi)->_AtomIndex];
      BondOrder order = (*bi)->_BondOrder;
      fromAtom->bondTo(toAtom,order);
    }
  }
  return res;
}
#endif

// Stuff

};
