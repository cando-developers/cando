/*
    File: residue.cc
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
//
// (C) 2004 Christian E. Schafmeister
//



#include <cando/chem/residue.h>

#include <iostream>
#include <vector>
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/matter.h>
#include <cando/chem/monomer.h>
//#include "core/archiveNode.h"
//#include <cando/chem/candoDatabaseReference.h>
//#include <cando/adapt/stringSet.h>
#include <cando/chem/loop.h>
//#include <cando/adapt/symbolSet.h>
#include <cando/chem/restraint.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/calculatePosition.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/constitution.h>
#include <cando/chem/bond.h>
#include <cando/chem/atom.h>
#include <cando/chem/atomIdMap.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>


namespace chem {


bool Residue_O::applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value ) {
  if ( this->Base::applyPropertyToSlot(prop,value) ) return true;
  return false;
}


#if 0
CL_DEFMETHOD bool	Residue_O::Residue_equal(core::T_sp obj) const
{
    if ( this->eq(obj) ) return true;
    if ( !obj.isA<Residue_O>() ) return false;
    Residue_sp other = obj.as<Residue_O>();
    if ( other->getName() != this->getName() ) return false;
    if ( other->_Contents.size() != this->_Contents.size() ) return false;
    Matter_O::const_contentIterator tit,oit;
    for ( tit=this->_Contents.begin(), oit=other->_Contents.begin();
	    tit!=this->_Contents.end(); tit++, oit++ )
    {
	if ( ! (*tit)->equal(*oit) ) return false;
    }
    return true;
}
#endif

void Residue_O::transferCoordinates(Matter_sp obj)
{
    if ( !obj.isA<Residue_O>() ) 
    {
	SIMPLE_ERROR("You can only transfer coordinates to a Residue from another Residue");
    }
    Residue_sp other = obj.as<Residue_O>();
    if ( other->_Contents.size() != this->_Contents.size() )
    {
	SIMPLE_ERROR("You can only transfer coordinates if the two residues have the same number of atoms");
    }
    Matter_O::contentIterator tit,oit;
    for ( tit=this->_Contents.begin(), oit=other->_Contents.begin();
	    tit!=this->_Contents.end(); tit++, oit++ )
    {
	(*tit)->transferCoordinates(*oit);
    }
}



CL_LISPIFY_NAME("addAtom");
CL_DEFMETHOD void	Residue_O::addAtom(Atom_sp a)
{
    this->addMatter(a);
    LOG("Added {} to {}" , a->description().c_str() , this->description().c_str() );
}

#if 0
bool Residue_O::recognizesMonomerAlias(core::Symbol_sp s)
{
    if ( this->_MonomerAliases.nilp() ) return false;
    return this->_MonomerAliases->contains(s);
}

void	Residue_O::setMonomerAliases(adapt::SymbolSet_sp s)
{
    this->_MonomerAliases = s;
}

CL_LISPIFY_NAME("addMonomerAlias");
CL_DEFMETHOD void Residue_O::addMonomerAlias(core::Symbol_sp s)
{
    if (this->_MonomerAliases.nilp() )
    {
	this->_MonomerAliases = adapt::SymbolSet_O::create();
    }
    this->_MonomerAliases->insert(s);
}



adapt::SymbolSet_sp Residue_O::getMonomerAliases()
{
    ASSERTNOTNULL(this->_MonomerAliases);
    return this->_MonomerAliases;
}
#endif


void	Residue_O::addVirtualAtom(MatterName name, CalculatePosition_sp proc)
{
    LOG("Residue_O::addVirtualAtom adding virtual atom with name({}) in residue({})" , _rep_(name) , _rep_(this->getName())  );
    if ( this->hasAtomWithName(name) )
    {
	stringstream ss;
	ss << "The " << this->description() << " already has the virtual atom: " << name;
	SIMPLE_ERROR("{}" , ss.str());
    }
    VirtualAtom_sp va = VirtualAtom_O::create(name,proc);
    va->setElement(element_Dummy);
    this->addMatter(va);
    va->calculatePosition(this->asSmartPtr());
}


void	Residue_O::fields( core::Record_sp node )
{
  LOG("Status" );
  node->field_if_not_unbound(INTERN_(kw,type),this->_Type);
  node->field_if_not_nil( INTERN_(kw,pdbName),this->_PdbName);
  node->field( INTERN_(kw,uniqueLabel),this->_UniqueLabel);
  node->/*pod_*/field_if_not_default( INTERN_(kw,NetCharge),this->_NetCharge,0);
  node->/*pod_*/field_if_not_default( INTERN_(kw,fileSeqNum),this->_FileSequenceNumber,UndefinedUnsignedInt);
  node->field_if_not_nil( INTERN_(kw,monomerAliases),this->_MonomerAliases);
  LOG("Status" );
  switch (node->stage()) {
  case core::Record_O::saving: {
	    // Accumulate intraresidue bonds into a vector
#if 1
    core::HashTable_sp atomToResidue = this->atomToResidueMap();
//    BondList_sp bondList = BondList_O::create();
    auto bondList  = gctools::GC<BondList_O>::allocate_with_default_constructor();
    { 
      for ( auto aa = this->begin_atoms(); aa != this->end_atoms(); ++aa ) {
        Atom_sp a = gc::As_unsafe<Atom_sp>(*aa);
        a->addUniqueIntraResidueBondCopiesToBondList(atomToResidue,bondList);
      }
    }
    node->field( INTERN_(kw,bl),bondList);
    ASSERTNOTNULL(bondList);
#endif    
  }
      break;
  case core::Record_O::initializing:
  case core::Record_O::loading: {
#if 1
    
    LOG("Creating the intraResidue bonds" );
	    // create the intraResidue bonds
    LOG("About to load bondList" );
    BondList_sp bondList = BondList_O::create();
    node->field( INTERN_(kw,bl),bondList);
    ASSERTNOTNULL(bondList);
    RECORD_LOG("residue bondList = {}" , _rep_(bondList));
    bondList->imposeYourself();
#endif
  }
  case core::Record_O::patching: {
    // Nothing should need to be done
  }
      break;
  }
  this->Base::fields(node);
}



//
// Constructor
//


void	Residue_O::initialize()
{
    this->Base::initialize();
//    this->_Constitution = _Nil<CandoDatabaseReference_O>();
    this->_FileSequenceNumber = UndefinedUnsignedInt;
    this->_NetCharge = 0;
    this->_MonomerAliases = nil<core::T_O>();
    this->_UniqueLabel = nil<core::T_O>();
    this->_PdbName = nil<core::Symbol_O>();
}

//
// Copy Constructor
//
Residue_O::Residue_O(const Residue_O& res)
	:Matter_O(res)
{
    this->_Selected = res._Selected;
    this->_TempInt = res._TempInt;
    this->_PdbName = res._PdbName;
//    ANN(res._Constitution);
//    this->_Constitution = res._Constitution;
    this->_FileSequenceNumber = res._FileSequenceNumber;
    this->_MonomerAliases = res._MonomerAliases;
    this->_UniqueLabel = res._UniqueLabel;
    this->_Type = res._Type;
}



//
// Constructor
//




void	Residue_O::duplicateFrom( const Residue_O* r )
{
    *this = *r;
}



Matter_sp	Residue_O::copyDontRedirectAtoms(core::T_sp new_to_old)
{
    contentIterator	a;
    Atom_sp				acopy;
    Atom_sp				aorig;

    auto  rPNew = gctools::GC<Residue_O>::copy( *this ); // = RP_Copy<Residue_O>(this); // _copy_Residue_sp(this);
    if (gc::IsA<core::HashTable_sp>(new_to_old)) {
      core::HashTable_sp new_to_old_ht = gc::As_unsafe<core::HashTable_sp>(new_to_old);
      new_to_old_ht->setf_gethash(rPNew,this->asSmartPtr());
    }
//    rPNew = Residue_sp(n_ew Residue_O(*this));
//    rPNew->duplicateFrom(this); //    *rPNew = *this;
    rPNew->eraseContents();
    for ( a=this->begin_atoms(); a!=this->end_atoms(); a++ ) 
    {
	aorig = (*a).as<Atom_O>();
	acopy = aorig->copyDontRedirectAtoms(new_to_old).as<Atom_O>();
        if (gc::IsA<core::HashTable_sp>(new_to_old)) {
          core::HashTable_sp new_to_old_ht = gc::As_unsafe<core::HashTable_sp>(new_to_old);
          new_to_old_ht->setf_gethash(acopy,aorig);
        }
	LOG("Copying atom({}) with {} bonds" , aorig->getName().c_str() , aorig->numberOfBonds()  );
	rPNew->addMatter((Matter_sp)(acopy));
	LOG("Completed copy for new {}" , acopy->description().c_str()  );
    }
    rPNew->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
    return rPNew;
}


void Residue_O::redirectAtoms()
{
    for ( contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ )
    {
	Atom_sp at = (*a).as<Atom_O>();
	LOG("Redirecting bonds for {} with {} bonds (only bonds where we are atom1)"
	    , at->description() , at->numberOfBonds() );
	at->redirectAtoms();
    }
    this->redirectRestraintAtoms();
}



CL_LISPIFY_NAME("copy");
CL_DEFMETHOD Matter_sp Residue_O::copy(core::T_sp new_to_old)
{ 
    Residue_sp newRes = this->copyDontRedirectAtoms(new_to_old).as<Residue_O>();
    Loop lbonds(this->asSmartPtr(),BONDS);
    while (lbonds.advanceLoopAndProcess()) {
      Bond_sp b = lbonds.getBond();
      Bond_sp bcopy = b->copyDontRedirectAtoms();
      if (chem__verbose(1)) {
        core::clasp_write_string(fmt::format("    b -> {}\n" , _rep_(b)));
        core::clasp_write_string(fmt::format("bcopy -> {}\n" , _rep_(bcopy)));
      }
      bcopy->addYourselfToCopiedAtoms();
    }
    newRes->redirectAtoms();
    return newRes;
}



CL_LISPIFY_NAME("removeAtomsWithNames");
CL_DEFMETHOD void Residue_O::removeAtomsWithNames(core::List_sp args)
{
  for ( auto c : args ) {
    MatterName atomName = oCar(c).as<MatterName::Type>();
    Atom_sp a = gc::As_unsafe<Atom_sp>(this->atomWithName(atomName));
    this->removeAtomDeleteBonds(a);
  }
}

//
//	removeAtomDeleteBonds
//
//	Remove the atom from the residue and delete its bonds
CL_LISPIFY_NAME("removeAtomDeleteBonds");
CL_DEFMETHOD void	Residue_O::removeAtomDeleteBonds(Atom_sp a)
{
  contentIterator	atom;
  Atom_sp		aTemp;
  LOG("Residue_O::removeAtomsDeleteBonds| Removing a:{} from r:{}" , &a , this  );
  for ( atom=this->begin_atoms(); atom != this->end_atoms(); atom++ ) {
    aTemp = (*atom).as<Atom_O>();
    if ( aTemp == a ) {
      this->eraseContent(atom);
      aTemp->removeAllBonds();
      LOG("Residue_O::removeAtomDeleteBonds setting atom {} parent to null" , &a  );
      return;
    }
  }
  SIMPLE_ERROR("Residue does not contain atom for removal");
}



CL_LISPIFY_NAME("containsAtom");
CL_DEFMETHOD bool	Residue_O::containsAtom(Atom_sp a)
{
  for (auto ii=this->begin_contents(); ii!=this->end_contents(); ++ii) {
    Atom_sp aa = gc::As_unsafe<Atom_sp>(*ii);
    if (aa == a) return true;
  }
  return false;
}



//
//	removeAtomDontDeleteBonds
//
//	Remove the atom from the residue and dont delete its bonds
//	This does not delete the atom!
CL_LISPIFY_NAME("removeAtomDontDeleteBonds");
CL_DEFMETHOD void	Residue_O::removeAtomDontDeleteBonds(Atom_sp a)
{
contentIterator	atom;
Atom_sp				aTemp;
    LOG("Residue_O::removeAtomsDontDeleteBonds| Removing a:{} from r:{}" , &a , this  );
    for ( atom=this->begin_atoms();
		atom != this->end_atoms(); atom++ ) {
      aTemp = (*atom).as<Atom_O>();
	LOG("Looking for atom name: {} in residue with name: {}" , a->getName().c_str() , aTemp->getName().c_str()  );
	if ( aTemp == a ) {
	LOG("     erasing" );
	    this->eraseContent(atom);
LOG("Residue_O::removeAtomDeleteBonds setting atom {} parent to null" , &a  );
//	    delete (aTemp);
	    return;
	}
    }
    SIMPLE_ERROR("You are asking to remove an atom from a residue that doesn't contain it");
}


size_t Residue_O::fillInImplicitHydrogens()
{
    gctools::Vec0<Atom_sp>	atoms;
contentIterator	aCur;
	/* fillInImplicitHydrogens will add more atoms to 
	 * the residue so we can't just iterate over the contents
	 * as we are adding to the contents.
	 * Pull out the atoms into a separate vector first and
	 * then iterate over that
	 */
    for ( aCur=this->_Contents.begin();aCur!=this->_Contents.end(); aCur++ ) 
    {
	atoms.push_back((*aCur).as<Atom_O>());
    }
    size_t count = 0;
    for (gctools::Vec0<Atom_sp>::iterator it=atoms.begin(); it!=atoms.end(); it++ )
    {
      count += (*it)->fillInImplicitHydrogensWithResidue(this->asSmartPtr());
    }
    return count;
}



void Residue_O::ensureAllAtomNamesAreUnique() const
{
  core::HashTableEq_sp names = core::HashTableEq_O::create_default();
  Loop lAtoms(this->asSmartPtr(),ATOMS);
  while (lAtoms.advanceLoopAndProcess()) {
    Atom_sp atom = lAtoms.getAtom();
    if (names->gethash(atom->getName()).notnilp()) {
      SIMPLE_ERROR("To define a constitution from a residue - all atom names must be unique - but the name {} was found to be duplicated" , _rep_(atom->getName()));
    }
    names->setf_gethash(atom->getName(),atom);
  }
}




#if 0
void CDFragment_O::uniqifyResidueAtomNames(Molecule_sp mol, bool verbose)
{
  Loop lResidues(mol,RESIDUES);
  size_t idx = 0;
  while (lResidues.advanceLoopAndProcess()) {
    Residue_sp res = lResidues.getResidue();
    Loop lAtoms(res,ATOMS);
    core::HashTableEq_sp namesToIndices = core::HashTableEq_O::create_default();
    while (lAtoms.advanceLoopAndProcess()) {
      Atom_sp atom = lAtoms.getAtom();
      if (namesToIndices->gethash(atom->getName()).notnilp()) {
        size_t nameTries = 0;
        core::Symbol_sp sym;
        do {
          stringstream ss;
          ss << atom->getName()->symbolNameAsString() << "_" << idx;
          sym = chemkw_intern(ss.str());
          if (namesToIndices->gethash(sym).nilp()) break;
        } while (++nameTries<1024);
        if (nameTries>=1024) {
          SIMPLE_ERROR("Tried to uniqify an atom name more than 1024 times");
        }
        atom->setName(sym);
        namesToIndices->setf_gethash(sym,sym);
      } else {
        namesToIndices->setf_gethash(atom->getName(),atom->getName());
      }
    }
  }
}
#endif
/*!
 * Make all the atom names unique by suffixing the non-unique ones with
 * a number so that the combination of name+number is unique.
 */
void Residue_O::makeAllAtomNamesInEachResidueUnique()
{
  Loop lAtoms(this->asSmartPtr(),ATOMS);
  core::HashTableEq_sp namesToIndices = core::HashTableEq_O::create_default();
  size_t index = 0;
  while (lAtoms.advanceLoopAndProcess()) {
    Atom_sp atom = lAtoms.getAtom();
    if (namesToIndices->gethash(atom->getName()).notnilp()) {
      size_t nameTries = 0;
      core::Symbol_sp sym;
      do {
        stringstream ss;
        ss << atom->getName()->symbolNameAsString() << '_' << index;
        sym = chemkw_intern(ss.str());
        if (namesToIndices->gethash(sym).nilp()) break;
        index++;
      } while (++nameTries<1024);
      if (nameTries>=1024) {
        SIMPLE_ERROR("Tried to uniqify an atom name more than 1024 times");
      }
      atom->setName(sym);
      namesToIndices->setf_gethash(sym,sym);
    } else {
      namesToIndices->setf_gethash(atom->getName(),atom->getName());
    }
  }
#if 0
    contentIterator ai;
    adapt::SymbolSet_sp allNames = adapt::SymbolSet_O::create();
    adapt::SymbolSet_sp allNamesAccumulate = adapt::SymbolSet_O::create();
    core::HashTableEq_sp atomsThatShareName = core::HashTableEq_O::create_default();
//    multimap<string,Atom_sp>	atomsThatShareName;
    for ( ai=this->begin_atoms();
          ai != this->end_atoms(); ai++ ) {
	Atom_sp atom = (*ai).as<Atom_O>();
	LOG("Looking at atom name({})" , _rep_(atom->getName())  );
        allNames->insert(atom->getName());
	allNamesAccumulate->insert(atom->getName());
        core::List_sp atoms = atomsThatShareName->gethash(atom->getName()
                                                          ,_Nil<core::T_O>()).as<core::Cons_O>();
        atoms = core::Cons_O::create(atom,atoms);
        atomsThatShareName->setf_gethash(atom->getName(),atoms);
    }
    allNames->map( [&allNames,&allNamesAccumulate,&atomsThatShareName] (core::Symbol_sp name) {
            core::List_sp atomsWithName = atomsThatShareName->gethash(name);
            uint numberSuffix = 1;
            for ( auto cur : atomsWithName ) {
                stringstream newName;
                uint fails = 0;
                while ( 1 ) 
                {
                    Atom_sp atom = oCar(cur).as<Atom_O>();
                    newName.str("");
                    newName << atom->getName();
                    newName << numberSuffix;
                    // Check if this is a unique name
                    MatterName symNewName = chemkw_intern(newName.str());
                    if ( allNamesAccumulate->contains(symNewName) )
                    {
			// We don't have that name so we can use it and
			// move on
                        atom->setName(symNewName);
                        allNamesAccumulate->insert(symNewName);
                        break;
                    } else fails++;
                    numberSuffix++;
                    if ( fails > 1000 )
                    {
                        SIMPLE_ERROR("There were more than 1000 fails trying to identify a unique atom name");
                    }
                }
            }
        } );
#endif
}


void Residue_O::setAliasesForAtoms(core::List_sp aliasAtoms, core::List_sp atomAliases)
{
  ASSERT_eq(core::cl__length(aliasAtoms),core::cl__length(atomAliases));
    while ( aliasAtoms.notnilp() )
    {
	MatterName aliasAtom = oCar(aliasAtoms).as<MatterName::Type>();
	MatterName atomAlias = oCar(atomAliases).as<MatterName::Type>();
	Atom_sp a = gc::As_unsafe<Atom_sp>(this->atomWithName(aliasAtom));
	a->setAlias(atomAlias);
	aliasAtoms = oCdr(aliasAtoms);
	atomAliases = oCdr(atomAliases);
    }
}


//
//	getAllUniqueAtomNames
//
//	Return a set of all unique atom names
#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("getAtomNamesAsSymbolSet");
CL_DEFMETHOD adapt::SymbolSet_sp	Residue_O::getAtomNamesAsSymbolSet()
{
contentIterator	atom;
Atom_sp				aTemp;
adapt::SymbolSet_sp unique = adapt::SymbolSet_O::create();
    for ( atom=this->begin_atoms();
		atom != this->end_atoms(); atom++ ) {
	aTemp = (*atom).as<Atom_O>();
	unique->insert(aTemp->getName());
    }
    return unique;
}
#endif

//
//	getAllUniqueAtomNames
//
//	Return a set of all unique atom names
#ifdef USE_TOPOLOGY
adapt::SymbolSet_sp	Residue_O::getAllUniqueAtomNames()
{
    adapt::SymbolSet_sp unique;
contentIterator	atom;
Atom_sp				aTemp;
    for ( atom=this->begin_atoms();
		atom != this->end_atoms(); atom++ ) {
	aTemp = (*atom).as<Atom_O>();
	unique->insert(aTemp->getName());
    }
    return unique;
}
#endif

//
//	getAtoms
//
//	Return a set of all atoms
VectorAtom	Residue_O::getAtoms()
{
VectorAtom			atoms;
Atom_sp				aTemp;
contentIterator	atom;
    LOG("Residue_O::allAtoms" );
    for ( atom=this->begin_atoms();
		atom != this->end_atoms(); atom++ ) {
	aTemp = (*atom).as<Atom_O>();
	atoms.push_back(aTemp);
    }
    return atoms;
}




core::HashTable_sp Residue_O::atomToResidueMap()
{
  core::HashTableEq_sp map = core::HashTableEq_O::create_default();
  for ( auto aa = this->begin_atoms(); aa!=this->end_atoms(); ++aa ) {
    map->setf_gethash(*aa,this->asSmartPtr());
  }
  return map;
}

//
//	getOutGoingBonds
//
//	Return a vector of bonds that go out of the residue
//
BondList_sp	Residue_O::getOutGoingBonds()
{
  BondList_sp			bonds = BondList_O::create();
  contentIterator	aPPCur;
  Atom_sp				a;
  core::HashTable_sp atomToResidue = this->asSmartPtr()->atomToResidueMap();
  for ( aPPCur=this->begin_atoms();
        aPPCur != this->end_atoms(); aPPCur++ ) {
    a = (*aPPCur).as<Atom_O>();
    a->addInterResidueBondsToBondList(atomToResidue,bonds);
  }
  return bonds;
}


Atom_sp Residue_O::atomWithAliasOrNil(core::Symbol_sp alias)
{
contentIterator	aCur;
    for ( aCur=this->_Contents.begin();aCur!=this->_Contents.end(); aCur++ )
    {
	Atom_sp a = (*aCur).as<Atom_O>();
	LOG("Looking at({}) for alias({})" , (*aCur)->getName().c_str() , _rep_(alias)  );
        if ( a->getAlias() == alias)
	{
            return a;
        }
    }
    return nil<Atom_O>();
}

Atom_sp Residue_O::atomWithAlias(core::Symbol_sp alias)
{
    Atom_sp a = this->atomWithAliasOrNil(alias);
    if ( a.notnilp() ) return a;
    LOG("Matter({}) with {} contents does not contain content with alias({})" , this->name.c_str() , this->_Contents.size() , _rep_(alias)  );
    SIMPLE_ERROR("residue ({}) does not contain atom with alias({})", this->_Name, alias);
}


bool	Residue_O::invalid()
{
contentIterator	aPPCur;
Atom_sp				a;
    for ( aPPCur=this->begin_atoms();
		aPPCur != this->end_atoms(); aPPCur++ ) {
      a = (*aPPCur).as<Atom_O>();
	if ( a->invalid() ) return true;
    }
    return false;
}

CL_LAMBDA((residue !)atom-name &optional (errorp t));
CL_DEFMETHOD
core::T_mv Residue_O::atomWithName(MatterName name, bool errorp)
{
  contentIterator	aCur;
  for ( aCur=this->_Contents.begin();aCur!=this->_Contents.end(); aCur++ ) {
    LOG("Looking at({}) for({})" , (*aCur)->getName().c_str() , sName.c_str()  );
    if ( (*aCur)->getName() == name ) {
      return Values(*aCur,_lisp->_true());
    }
  }
  if (!errorp) {
    return Values(nil<core::T_O>(),nil<core::T_O>());
  }
  SIMPLE_ERROR("{} ({}) does not contain name({})", this->className(), this->_Name, name);
}

Vector3 Residue_O::positionOfAtomWithName(MatterName name)
{
  Atom_sp a = gc::As_unsafe<Atom_sp>(this->atomWithName(name));
    return a->getPosition();
}
void	Residue_O::failIfInvalid()
{
    if ( this->invalid() ) {
      SIMPLE_ERROR("INVALID {}" , this->description());
    }
}


CL_LISPIFY_NAME("useAtomCoordinatesToDefineAnchors");
CL_DEFMETHOD void	Residue_O::useAtomCoordinatesToDefineAnchors()
{
    IMPLEMENT_ME();
}

uint Residue_O::numberOfAtoms() const
{
    return this->_Contents.size();
}


    AtomIdMap_sp Residue_O::buildAtomIdMap() const
    {
	AtomIdMap_sp atomIdMap = AtomIdMap_O::create();
	atomIdMap->resizeAggregate(1);
	int mid = 0;
	atomIdMap->resizeMolecule(mid,1);
	int rid = 0;
	int numAtoms = this->_Contents[mid]->_Contents[rid]->_Contents.size();
	atomIdMap->resizeResidue(mid,rid,numAtoms);
	for ( int aid=0; aid<numAtoms; aid++ )
	{
	    AtomId atomId(mid,rid,aid);
	    atomIdMap->set(atomId,this->_Contents[mid]->_Contents[rid]->_Contents[aid]);
	}
	return atomIdMap;
    }



    Atom_sp Residue_O::atomWithAtomId(const AtomId& atomId) const
    {
	int aid = atomId.atomId();
	if ( aid >=0 && aid <=(int)this->_Contents.size() )
	{
	    Atom_sp atom = this->_Contents[aid].as<Atom_O>();
	    return atom;
	}
	SIMPLE_ERROR("Illegal atomId[{}] must be less than {}" , aid , this->_Contents.size() );
    }







CL_LISPIFY_NAME(make-residue);
DOCGROUP(cando);
CL_LAMBDA(&optional (name nil) atoms)
CL_DEFUN Residue_sp Residue_O::make(core::Symbol_sp name, core::List_sp atoms)
{
  auto me = gctools::GC<Residue_O>::allocate_with_default_constructor();
  me->setName(name);
  for ( auto cur : atoms ) {
    auto atm = gc::As<Atom_sp>(CONS_CAR(cur));
    me->addMatter(atm);
  }
    return me;
};





};
