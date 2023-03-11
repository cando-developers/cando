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
#include <clasp/core/ql.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/lispStream.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/matter.h>
#include <cando/chem/residue.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/atom.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/bond.h>
#include <clasp/core/wrappers.h>

namespace chem
{

void ConstitutionBond_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,toatom),this->_ToAtomIndex);
  node->field(INTERN_(kw,bond_order),this->_BondOrder);
  // this->Base::fields(node);
}

string ConstitutionBond_O::__repr__() const
{
  stringstream ss;
  ss << "(" << this->className() << " ";
  ss << " :toAtom " << this->_ToAtomIndex;
  ss << " :bondOrder " << bondOrderToString(this->_BondOrder);
  ss << ">";
  return ss.str();
}

void ConstitutionAtom_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,name),this->_AtomName);
  node->field(INTERN_(kw,index),this->_Index);
  node->field(INTERN_(kw,element),this->_Element);
  node->field(INTERN_(kw,type),this->_AtomType);
  node->field_if_not_default(INTERN_(kw,stereochemType), this->_StereochemistryType, undefinedCenter);
  node->field(INTERN_(kw,bonds),this->_Bonds);
  node->field_if_not_nil(INTERN_(kw,properties),this->_Properties);
  this->Base::fields(node);
}  

CL_DOCSTRING(R"(Build a chem:atom from this chem:constitution-atom.)");
CL_DEFMETHOD Atom_sp ConstitutionAtom_O::buildAtom() const
{
  Atom_sp atom = Atom_O::create();
  atom->setName(this->_AtomName);
  atom->setElement(this->_Element);
  atom->setStereochemistryType(this->_StereochemistryType);
  atom->setProperties(core::cl__copy_list(this->_Properties));
  return atom;
}

CL_LISPIFY_NAME("atomName");
CL_DEFMETHOD 	MatterName ConstitutionAtom_O::atomName() const {
  return this->_AtomName;
};

CL_LISPIFY_NAME("atomType");
CL_DEFMETHOD 	core::T_sp ConstitutionAtom_O::atomType() const
{
  return this->_AtomType;
};


string ConstitutionAtom_O::__repr__() const
{
  stringstream ss;
  ss << "(" << this->className() << " ";
  ss << " :atomName " << this->_AtomName;
  ss << " :element " << symbolFromElement(this->_Element);
  if (this->_StereochemistryType != undefinedCenter) {
    ss << " :stereochemistry-type " << _rep_(gc::As<core::SymbolToEnumConverter_sp>(_sym__PLUS_stereochemistryTypeConverter_PLUS_->symbolValue())->symbolForEnum(this->_StereochemistryType));
  }
  ss << ")";
  return ss.str();
}

void ConstitutionVirtualAtom_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,calculate_position),this->_CalculatePositionCode);
  this->Base::fields(node);
}  


CL_DOCSTRING(R"(Build a chem:virtual-atom from this chem:constitution-virtual-atom.)");
CL_DEFMETHOD Atom_sp ConstitutionVirtualAtom_O::buildAtom() const
{
  VirtualAtom_sp atom = VirtualAtom_O::create(this->_AtomName,this->_CalculatePositionCode);
  atom->setElement(this->_Element);
  atom->setStereochemistryType(this->_StereochemistryType);
  atom->setProperties(core::cl__copy_list(this->_Properties));
  return atom;
}

  

CL_LISPIFY_NAME(makeConstitutionAtomsFromResidue);
CL_LAMBDA(residue &optional verbose);
DOCGROUP(cando);
CL_DEFUN ConstitutionAtoms_sp ConstitutionAtoms_O::makeConstitutionAtomsFromResidue(Residue_sp residue, bool verbose)
{
  residue->ensureAllAtomNamesAreUnique();
  auto  catoms  = gctools::GC<ConstitutionAtoms_O>::allocate_with_default_constructor();
  Matter_O::contentIterator it;
  ConstitutionAtomIndex0N index = 0;
  MapAtomsToConstitutionAtomIndex0N atomToIndexMap;
  {
    for ( it=residue->begin_atoms(); it!=residue->end_atoms(); it++, index++ )
    {
      Atom_sp atom = (*it).as<Atom_O>();
      ConstitutionAtom_sp catom = atom->asConstitutionAtom(index);
      catoms->_Atoms.push_back(catom);
      atomToIndexMap[atom] = index;
      if (verbose) core::write_bf_stream(fmt::sprintf("Atom %s index: %d\n" , _rep_(atom) , index));
    }
  }
  {
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
{
  this->_Atoms.push_back(ca);
};


ConstitutionAtom_sp& ConstitutionAtoms_O::operator[](ConstitutionAtomIndex0N idx)
{
  ASSERTF(idx>=0&&idx<(int)(this->_Atoms.size()),("index[%d] out of range") , idx );
  ConstitutionAtom_sp& ca = this->_Atoms[(uint)idx];
//  ASSERTF(idx==ca->_AtomIndex,("A mismatch has occured between the index[%d] in the ConstitutionAtoms array and the ConstitutionAtom index[%d] - these have to match or we can't quickly look up atoms by their index") , idx , ca->_AtomIndex );
  return ca;
}

CL_LISPIFY_NAME("atomNamesAsSymbolSet");
CL_DEFMETHOD     adapt::SymbolSet_sp ConstitutionAtoms_O::atomNamesAsSymbolSet()
{
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


CL_DOCSTRING(R"dx(Return (values constitution-atom T) if the atom with NAME is found.
If it is not found and ERRORP is NIL then return (values NIL NIL), otherwise signal an error.)dx")
CL_LAMBDA(atom-name &optional (errorp t));
CL_LISPIFY_NAME("atomWithName");
CL_DEFMETHOD core::T_mv ConstitutionAtoms_O::atomWithName(MatterName nm, bool errorp )
{
  for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++ )
  {
    if ( (*ci)->_AtomName == nm ) return Values(*ci,_lisp->_true());
  }
  if (!errorp) {
    return Values(nil<core::T_O>(),nil<core::T_O>());
  }
  SIMPLE_ERROR(("Could not find ConstitutionAtom with name[%s]") , core::_rep_(nm) );
}


CL_LISPIFY_NAME("atomWithId");
CL_DEFMETHOD     ConstitutionAtom_sp ConstitutionAtoms_O::atomWithId(ConstitutionAtomIndex0N id) const
{
  return this->_Atoms[id];
}

CL_NAME(CHEM:CONSTITUTION-ATOMS/INDEX);
CL_DEFMETHOD     int ConstitutionAtoms_O::index(MatterName name) const
{
  for ( int idx=0; idx<(int)this->_Atoms.size(); idx++ )
  {
    const ConstitutionAtom_sp& atom = this->_Atoms[idx];
    if ( atom->_AtomName == name ) return idx;
  }
  SIMPLE_ERROR(("Unknown atom[%s]") , core::_rep_(name) );
}

CL_DEFMETHOD core::List_sp ConstitutionAtoms_O::constitutionAtomsAsList() const {
  ql::list l;
  for ( int i=0; i<this->_Atoms.size(); ++i ) {
    l << this->_Atoms[i];
  }
  return l.cons();
}

CL_DOCSTRING(R"(Create a chem:residue that has all the atoms/bonds and all the necessary atom/bond
properties set properly for this chem:constitution-atoms.)")
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
    atom->setStereochemistryType((*ci)->_StereochemistryType);
    atom->setProperties(core::cl__copy_list((*ci)->_Properties));
    atoms[idx] = atom;
    res->putMatter(idx,atom);
  }
  idx = 0;
  for ( gctools::Vec0<ConstitutionAtom_sp>::const_iterator ci=this->_Atoms.begin();
        ci!=this->_Atoms.end(); ci++, idx++ )
  {
    ASSERT_lt(idx,(int)atoms.size());
    Atom_sp fromAtom = atoms[idx];
    for ( gctools::Vec0<ConstitutionBond_sp>::const_iterator bi=(*ci)->_Bonds.begin(); bi!=(*ci)->_Bonds.end(); bi++ )
    {
		// Avoid making bonds twice
      if ( (*bi)->_ToAtomIndex < idx ) continue;
      ASSERT_lt((*bi)->_ToAtomIndex,(int)atoms.size());
      Atom_sp toAtom = atoms[(*bi)->_ToAtomIndex];
      BondOrder order = (*bi)->_BondOrder;
      fromAtom->bondTo(toAtom,order);
    }
  }
  return res;
}

// Stuff

};
