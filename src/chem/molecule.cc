/*
    File: molecule.cc
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



#include <clasp/core/common.h>
#include <clasp/core/record.h>
#include <clasp/core/bformat.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/molecule.h>
#include <cando/chem/loop.h>
#include <clasp/core/numerics.h>
//#include "core/serialize.h"
#include <cando/chem/restraint.h>
#include <cando/chem/bond.h>
#include <clasp/core/translators.h>
#include <clasp/core/hashTable.h>
#include <cando/chem/atomIdMap.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>



namespace chem {

void	Molecule_O::initialize()
{
  this->Base::initialize();
}

bool Molecule_O::applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value ) {
  if ( this->Base::applyPropertyToSlot(prop,value) ) return true;
  return false;
}

string Molecule_O::__repr__() const
{
  return this->Base::__repr__();
}

void Molecule_O::fields(core::Record_sp node)
{
  node->field_if_not_unbound(INTERN_(kw,force_field_name),this->_ForceFieldName);
  node->field_if_not_unbound(INTERN_(kw,type),this->_Type);
  switch (node->stage()) {
  case core::Record_O::saving: {
	    // Accumulate intraresidue bonds into a vector
#if 1
    core::HashTable_sp atomToResidue = this->atomToResidueMap();
    BondList_sp bondList = BondList_O::create();
    { 
      Loop lmol(this->asSmartPtr(),ATOMS);
      while (lmol.advance()) {
        Atom_sp a = lmol.getAtom();
        a->addUniqueInterResidueBondCopiesToBondList(atomToResidue,bondList);
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
    RECORD_LOG("residue bondList = %s" , _rep_(bondList));
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

Molecule_O::Molecule_O(const Molecule_O& mol)
  :Matter_O(mol), _ForceFieldName(mol._ForceFieldName)
{
  this->_Type = mol._Type;
}

//
// Destructor
//

//
//	addResidue
//
//	Add the residue to the end of the molecule
//
void Molecule_O::addResidue( Matter_sp r )
{
  this->addMatter( r );
  LOG("Added %s to %s" , r->description().c_str() , this->description().c_str()  );
}

//
//	addResidueRetainId
//
//	Add the residue to the end of the molecule
//
void Molecule_O::addResidueRetainId( Matter_sp r )
{
  LOG("Adding %s to %s" , r->description().c_str() , this->description().c_str()  );
  this->addMatterRetainId( r );
}





//
//	removeResidue
//
//	Remove the residue
//
CL_LISPIFY_NAME("removeResidue");
CL_DEFMETHOD void Molecule_O::removeResidue( Matter_sp a )
{
  for ( auto it=this->begin_residues(); it!= this->end_residues(); it++ ) {
    if ( (*it).as<Residue_O>() == a ) {
      this->eraseContent(it);
      return;
    }
  }
  SIMPLE_ERROR(("removeResidue: Molecule does not contain residue: %s") , _rep_(a->getName()) );
}

#if 0

bool	Molecule_O::equal(core::T_sp obj) const
    {
      if ( this->eq(obj) ) return true;
      if ( Molecule_sp other = obj.asOrNull<Molecule_O>() ) {
	if ( other->getName() != this->getName() ) return false;
	if ( other->_Contents.size() != this->_Contents.size() ) return false;
	Matter_O::const_contentIterator tit = this->_Contents.begin();
	Matter_O::const_contentIterator oit = this->_Contents.begin();
	for ( ; tit!=this->_Contents.end(); tit++, oit++ )
	{
          if ( ! (*tit)->equal(*oit) ) return false;
	}
	return true;
      }
      return false;
    }
#endif

void Molecule_O::transferCoordinates(Matter_sp obj)
    {
	if ( !obj.isA<Molecule_O>() ) 
	{
	    SIMPLE_ERROR(("You can only transfer coordinates to a Molecule from another Molecule"));
	}
	Molecule_sp other = obj.as<Molecule_O>();
	if ( other->length() != this->length() )
	{
	    SIMPLE_ERROR(("You can only transfer coordinates if the two Molecules have the same number of contents"));
	}
	Matter_O::contentIterator tit,oit;
	for ( tit=this->_Contents.begin(), oit=other->_Contents.begin();
	      tit!=this->_Contents.end(); tit++, oit++ )
	{
	    (*tit)->transferCoordinates(*oit);
	}
    }
CL_LISPIFY_NAME("moveAllAtomsIntoFirstResidue");
CL_DEFMETHOD     void	Molecule_O::moveAllAtomsIntoFirstResidue()
{
  SIMPLE_ERROR(("Reimplement moveAllAtomsIntoFirstResidue()"));
#if 0
  contentIterator	a;
  contentIterator	r;
  contentIterator	rHead;
  contentIterator	rRest;
  rHead = this->begin_residues();
  rRest = rHead+1;
  for ( r=rRest;r!=this->end_residues(); ) {
    for ( a = (*r)->begin_contents(); a!= (*r)->end_contents();) {
      (*a)->reparent(*rHead);
      a = (*r)->eraseContent(a);
    }
//	delete (*r);
    r = this->eraseContent(r);
  }
#endif
}



    void	Molecule_O::duplicate(const Molecule_O* a )
    {
	*this = *a;
    }


Matter_sp	Molecule_O::copyDontRedirectAtoms(core::T_sp new_to_old)
{
  Residue_sp			res;
  auto  newMol  = gctools::GC<Molecule_O>::copy( *this); // = RP_Copy<Molecule_O>(this);
  if (gc::IsA<core::HashTable_sp>(new_to_old)) {
    core::HashTable_sp new_to_old_ht = gc::As_unsafe<core::HashTable_sp>(new_to_old);
    new_to_old_ht->setf_gethash(newMol,this->asSmartPtr());
  }
  for ( const_contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ ) {
    res = (*a).as<Residue_O>();
    newMol->addMatter(res->copyDontRedirectAtoms(new_to_old));
  }
  newMol->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
  return newMol;
}

void	Molecule_O::redirectAtoms()
    {_OF();
	LOG("Molecule_O::redirectAtoms START" );
	for ( contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ )
	{
          Residue_sp res = (*a).as<Residue_O>();
	    res->redirectAtoms();
	}
	this->redirectRestraintAtoms();
	LOG("Molecule_O::redirectAtoms DONE" );
    }



CL_LISPIFY_NAME("copy");
CL_DEFMETHOD     Matter_sp Molecule_O::copy(core::T_sp new_to_old)
    {
	Molecule_sp	newMol;
	newMol = this->copyDontRedirectAtoms(new_to_old).as<Molecule_O>();
        Loop lbonds(this->asSmartPtr(),BONDS);
        while (lbonds.advanceLoopAndProcess()) {
          Bond_sp b = lbonds.getBond();
          // printf("%s:%d Molecule_O::copy copy bond %s\n", __FILE__, __LINE__, _rep_(b).c_str() );
          Bond_sp bcopy = b->copyDontRedirectAtoms();
          bcopy->addYourselfToCopiedAtoms();
        }
	newMol->redirectAtoms();
	return newMol;
    }
#if 0 //[

    void	Molecule_O::dump()
    {
	contentIterator	rit;
	Residue_sp			res;
	printf( "Molecule: %s at: 0x%lx contains %d residues\n", this->getName().c_str(), this, this->getContents().size() );
	for ( rit=this->getContents().begin();
	      rit!=this->getContents().end(); rit++ ) {
          res = (*rit).as<Residue_O>();
	    res->dump();
	}
    }
#endif //]

VectorResidue	Molecule_O::getResiduesWithName(MatterName name ) {
	VectorResidue	result;
	Loop		lr;
	lr.loopTopMoleculeGoal( this->sharedThis<Molecule_O>(), RESIDUES );
	while ( lr.advanceLoopAndProcess() ) {
	    if ( lr.getResidue()->getName() == name ) {
		result.push_back(lr.getResidue());
	    }
	}
	return result;

    }


CL_LISPIFY_NAME("numberOfResiduesWithName");
CL_DEFMETHOD     int	Molecule_O::numberOfResiduesWithName( MatterName name )
    {
	VectorResidue	residues;
	residues = this->getResiduesWithName(name);
	return residues.size();
    }

CL_LISPIFY_NAME("getFirstResidueWithName");
CL_DEFMETHOD     Residue_sp	Molecule_O::getFirstResidueWithName(MatterName name)
    {
	VectorResidue residues = this->getResiduesWithName(name);
	if ( residues.size() > 0 ) {
	    return *(residues.begin());
	}
	SIMPLE_ERROR(("getFirstResidueWithName: Molecule does not contain residues with name: %s") , _rep_(name) );
    }



#ifdef RENDER
    geom::Render_sp	Molecule_O::rendered(core::List_sp opts)
    {
	GrPickableMatter_sp	rend;
	rend = GrPickableMatter_O::create();
	rend->setFromMatter(this->sharedThis<Molecule_O>());
	return rend;
    }
#endif

    uint	Molecule_O::numberOfAtoms()
    {
	Loop				lb;
	Vector3				v,vd;
	int				numberOfAtoms;
	numberOfAtoms = 0;
	lb.loopTopMoleculeGoal( this->sharedThis<Molecule_O>(), RESIDUES);
	while ( lb.advanceLoopAndProcess() ) 
	{
	    numberOfAtoms += lb.getResidue()->numberOfAtoms();
	}
	return numberOfAtoms;
    }


AtomIdMap_sp Molecule_O::buildAtomIdMap() const
{
  AtomIdMap_sp atomIdMap = AtomIdMap_O::create();
  atomIdMap->resizeAggregate(1);
  int mid = 0;
  int numResidues = this->_Contents.size();
  atomIdMap->resizeMolecule(mid,numResidues);
  for ( int rid =0; rid<numResidues; rid++ )
  {
    int numAtoms = this->_Contents[rid]->_Contents.size();
    core::write_bf_stream(fmt::sprintf("%s:%d rid %d of %d  numAtoms-> %d\n" , __FILE__ , __LINE__ , rid , numResidues , numAtoms ));
    atomIdMap->resizeResidue(mid,rid,numAtoms);
    for ( int aid=0; aid<numAtoms; aid++ )
    {
      AtomId atomId(mid,rid,aid);
      Atom_sp atom = this->_Contents[rid]->_Contents[aid].as<Atom_O>();
      core::write_bf_stream(fmt::sprintf("%s:%d Adding %d %d %d -> %s\n" , __FILE__ , __LINE__ , mid , rid , aid , _rep_(atom)));
      atomIdMap->set(atomId,atom);
    }
  }
  return atomIdMap;
}

    Atom_sp Molecule_O::atomWithAtomId(const AtomId& atomId) const
    {_OF();
	int resId = atomId.residueId();
	if ( resId >=0 && resId <=(int)this->_Contents.size() )
	{
	    Residue_sp residue = this->_Contents[resId].as<Residue_O>();
	    return residue->atomWithAtomId(atomId);
	}
	SIMPLE_ERROR(("Illegal residueId[%d] must be less than %d") , resId , this->_Contents.size() );
    }




#define ARGS_Molecule_O_make "(&key (name \"\"))"
#define DECL_Molecule_O_make ""
#define DOCS_Molecule_O_make "make Molecule args: &key name"
CL_LAMBDA(&optional (name nil))
CL_LISPIFY_NAME(make-molecule);
DOCGROUP(cando)
CL_DEFUN Molecule_sp Molecule_O::make(core::Symbol_sp name)
{
  auto me = gctools::GC<Molecule_O>::allocate_with_default_constructor();
    me->setName(name);
    return me;
};



CL_DEFMETHOD void Molecule_O::setf_force_field_name(core::T_sp name) {
  this->_ForceFieldName = name;
}

CL_DEFMETHOD core::T_sp Molecule_O::force_field_name() const {
  return this->_ForceFieldName;
}



}; // namespace chem
