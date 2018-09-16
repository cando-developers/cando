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
#include <cando/chem/molecule.h>
#include <cando/chem/loop.h>
#include <clasp/core/numerics.h>
//#include "core/serialize.h"
#include <cando/chem/restraint.h>
#include <cando/chem/bond.h>
#include <clasp/core/translators.h>
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
  Bond_sp	bond;
  Atom_sp		a;
  Atom_sp	a1, a2;
  BondOrder	o;
  Loop lb;
  switch (node->stage()) {
  case core::Record_O::saving: {
            // Accumulate intraresidue bonds into a vector
    _BLOCK_TRACE("Xmling inter-residue bonds");
    BondList_sp bondList = BondList_O::create();
    { _BLOCK_TRACE("Accumulating inter residue bond copies");
      LOG(BF("Serializing InterResidueBonds") );
      lb.loopTopMoleculeGoal( this->sharedThis<Molecule_O>(), BONDS );
      while ( lb.advanceLoopAndProcess() ) {
        a1 = lb.getBondA1();
        a2 = lb.getBondA2();
		    // LOG(BF("bond with atom1=%s") % a1->description().c_str()  );
		    // LOG(BF("bond with atom2=%s") % a2->description().c_str()  );
        o  = lb.getBondOrder();
        bond = Bond_O::create(a1,a2,o);
        if ( bond->isInterResidueBond() )
        {
          LOG(BF("Adding bond between atom1(%s) and atom2(%s)")
              % a1->description()
              % a2->description() );
          bondList->addBond(bond);
        }
      }
    }
    node->field( INTERN_(kw,bl),bondList);
  }
      break;
  case core::Record_O::initializing:
  case core::Record_O::loading: {
    	    // create the intraResidue bonds
    BondList_sp bondList;
    node->field( INTERN_(kw,bl),bondList);
    ASSERTNOTNULL(bondList);
    bondList->imposeYourself();
  }
      break;
  case core::Record_O::patching: {
  };
      break;
  }
  node->field_if_not_unbound(INTERN_(kw,type),this->_Type);
  this->Base::fields(node);
}

Molecule_O::Molecule_O(const Molecule_O& mol)
  :Matter_O(mol)
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
  LOG(BF("Added %s to %s") % r->description().c_str() % this->description().c_str()  );
}

//
//	addResidueRetainId
//
//	Add the residue to the end of the molecule
//
void Molecule_O::addResidueRetainId( Matter_sp r )
{
  LOG(BF("Adding %s to %s") % r->description().c_str() % this->description().c_str()  );
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
  gctools::Vec0<Matter_sp>::iterator	it;
  for ( it=this->getContents().begin(); it!= this->getContents().end(); it++ ) {
    if ( (*it).as<Residue_O>() == a ) {
      this->eraseContent(it);
      return;
    }
  }
  SIMPLE_ERROR(BF("removeResidue: Molecule does not contain residue: %s") % _rep_(a->getName()) );
}

    bool	Molecule_O::equal(core::T_sp obj) const
    {
      if ( this->eq(obj) ) return true;
      if ( Molecule_sp other = obj.asOrNull<Molecule_O>() ) {
	if ( other->getName() != this->getName() ) return false;
	if ( other->_contents.size() != this->_contents.size() ) return false;
	Matter_O::const_contentIterator tit = this->_contents.begin();
	Matter_O::const_contentIterator oit = this->_contents.begin();
	for ( ; tit!=this->_contents.end(); tit++, oit++ )
	{
          if ( ! (*tit)->equal(*oit) ) return false;
	}
	return true;
      }
      return false;
    }
    void Molecule_O::transferCoordinates(Matter_sp obj)
    {
	if ( !obj.isA<Molecule_O>() ) 
	{
	    SIMPLE_ERROR(BF("You can only transfer coordinates to a Molecule from another Molecule"));
	}
	Molecule_sp other = obj.as<Molecule_O>();
	if ( other->_contents.size() != this->_contents.size() )
	{
	    SIMPLE_ERROR(BF("You can only transfer coordinates if the two Molecules have the same number of contents"));
	}
	Matter_O::contentIterator tit,oit;
	for ( tit=this->_contents.begin(), oit=other->_contents.begin();
	      tit!=this->_contents.end(); tit++, oit++ )
	{
	    (*tit)->transferCoordinates(*oit);
	}
    }
CL_LISPIFY_NAME("moveAllAtomsIntoFirstResidue");
CL_DEFMETHOD     void	Molecule_O::moveAllAtomsIntoFirstResidue()
    {
	contentIterator	a;
	contentIterator	r;
	contentIterator	rHead;
	contentIterator	rRest;
	rHead = this->getContents().begin();
	rRest = rHead+1;
	for ( r=rRest;r!=this->getContents().end(); ) {
	    for ( a = (*r)->getContents().begin(); a!= (*r)->getContents().end();) {
		(*a)->reparent(*rHead);
		a = (*r)->eraseContent(a);
	    }
//	delete (*r);
	    r = this->eraseContent(r);
	}
    }



    void	Molecule_O::duplicate(const Molecule_O* a )
    {
	*this = *a;
    }


Matter_sp	Molecule_O::copyDontRedirectAtoms()
{
  Residue_sp			res;
  GC_COPY(Molecule_O, newMol , *this); // = RP_Copy<Molecule_O>(this);
  for ( const_contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ ) {
    res = (*a).as<Residue_O>();
    newMol->addMatter(res->copyDontRedirectAtoms());
  }
  newMol->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
  return newMol;
}

void	Molecule_O::redirectAtoms()
    {_OF();
	LOG(BF("Molecule_O::redirectAtoms START") );
	for ( contentIterator a=this->getContents().begin(); a!=this->getContents().end(); a++ )
	{
          Residue_sp res = (*a).as<Residue_O>();
	    res->redirectAtoms();
	}
	this->redirectRestraintAtoms();
	LOG(BF("Molecule_O::redirectAtoms DONE") );
    }



CL_LISPIFY_NAME("copy");
CL_DEFMETHOD     Matter_sp Molecule_O::copy()
    {
	Molecule_sp	newMol;
	newMol = this->copyDontRedirectAtoms().as<Molecule_O>();
        Loop lbonds(this->asSmartPtr(),BONDS);
        while (lbonds.advanceLoopAndProcess()) {
          Bond_sp b = lbonds.getBond();
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
	SIMPLE_ERROR(BF("getFirstResidueWithName: Molecule does not contain residues with name: %s")% _rep_(name) );
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


AtomIdToAtomMap_sp Molecule_O::buildAtomIdMap() const
{
  AtomIdToAtomMap_sp atomIdMap = AtomIdToAtomMap_O::create();
  atomIdMap->resize(1);
  int mid = 0;
  int numResidues = this->_contents.size();
  atomIdMap->resize(mid,numResidues);
  for ( int rid =0; rid<numResidues; rid++ )
  {
    int numAtoms = this->_contents[rid]->_contents.size();
    BFORMAT_T(BF("%s:%d rid %d of %d  numAtoms-> %d\n") % __FILE__ % __LINE__ % rid % numResidues % numAtoms );
    atomIdMap->resize(mid,rid,numAtoms);
    for ( int aid=0; aid<numAtoms; aid++ )
    {
      AtomId atomId(mid,rid,aid);
      Atom_sp atom = this->_contents[rid]->_contents[aid].as<Atom_O>();
      BFORMAT_T(BF("%s:%d Adding %d %d %d -> %s\n") % __FILE__ % __LINE__ % mid % rid % aid % _rep_(atom));
      atomIdMap->set(atomId,atom);
    }
  }
  return atomIdMap;
}

    Atom_sp Molecule_O::atomWithAtomId(const AtomId& atomId) const
    {_OF();
	int resId = atomId.residueId();
	if ( resId >=0 && resId <=(int)this->_contents.size() )
	{
	    Residue_sp residue = this->_contents[resId].as<Residue_O>();
	    return residue->atomWithAtomId(atomId);
	}
	SIMPLE_ERROR(BF("Illegal residueId[%d] must be less than %d") % resId % this->_contents.size() );
    }




#define ARGS_Molecule_O_make "(&key (name \"\"))"
#define DECL_Molecule_O_make ""
#define DOCS_Molecule_O_make "make Molecule args: &key name"
CL_LAMBDA(&optional (name nil));
CL_LISPIFY_NAME(make-molecule);
CL_DEFUN Molecule_sp Molecule_O::make(core::Symbol_sp name)
{
    GC_ALLOCATE(Molecule_O,me);
    me->setName(name);
    return me;
};






}; // namespace chem
