/*
    File: loop.cc
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





#include <stdio.h>
#include <clasp/core/common.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/designators.h>
#include <cando/chem/loop.h>
#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <clasp/core/wrappers.h>

namespace chem {

//
// pushSubLoop
//
// Push a subloop onto the subloop stack
void	Loop::pushSubLoop( Matter_sp c )
{
  LOG(BF("pushSubLoop: Starting a subloop over: %s") % c->description().c_str()  );
  this->curSubLoop++;
  this->subLoopTop[this->curSubLoop] = c;
  this->subLoopIteratorCur[this->curSubLoop] = c->_contents.begin();
  this->subLoopIteratorEnd[this->curSubLoop] = c->_contents.end();
  LOG(BF("pushSubLoop:  number of children count = %d") % c->_contents.size()  );
}


bool	Loop::advanceLoopAndProcess()
{
  if ( !this->done) {
    LOG(BF("Loop may still have values to generate") );
    this->advanceLoop();
    if ( this->done )
    {
      LOG(BF("The loop is done,returning false") );
      return false;
    }
    LOG(BF("The loop continues,returning true") );
    return true;
  }
  LOG(BF("The loop was done in the last step,returning false") );
  return false;
}




/*
 *      bLoopAtomVisible
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the atom flags satisfy the
 *      Visible/Invisible requirements of the loop.
 *      This means that the atom has ALL of the same flags
 *      set as set in fVisibleFlagsOn and
 *      ALL the same flags reset as SET in fVisibleFlagsOff.
 */
bool	Loop:: bLoopAtomVisible(Atom_sp aAtom )
{
  AtomFlags		fFlags;

  fFlags = aAtom->getFlags();
  if ( (fFlags&this->fVisibleFlagsOn)!=this->fVisibleFlagsOn )
    return(false);

  if ( (fFlags|(~(this->fVisibleFlagsOff)))!=(~(this->fVisibleFlagsOff)) )
    return(false);

			/* If the TempInt field is being used to determine */
			/* visibility, then check if the ATOM_CLASS has the */
			/* proper TempInt field */

  if ( this->fVisibilityFlags & TEMPINTUSED ) {
    if ( this->fVisibilityFlags & TEMPINTINVISIBLE ) {
      if ( aAtom->getTempInt() == this->iTempInt ) {
        return(false);
      }
    } else {
      if ( aAtom->getTempInt() != this->iTempInt ) {
        return(false);
      }
    }
  }

		/* Check if the ATOM_CLASS is the one invisible ATOM_CLASS for the loop */
		/* Used to construct spanning trees that do not pass a */
		/* certian ATOM_CLASS, good for looping over side chains */

  if ( this->aInvisibleAtom == aAtom ) return(false);

  return(true);
}




/*
 *      bSpanAtomVisible
 *
 *	Author:	Christian Schafmeister (1991)
 *
 *      Return TRUE if the atom is visible to a SPANNINGTREE loop.
 *      This means that the atom has not been seen before ( Set
 *      bPSeenBefore to TRUE if it has ) and that the atom
 *      has ALL of the same flags set as set in fVisibleFlagsOn and
 *      ALL the same flags reset as SET in fVisibleFlagsOff.
 */
bool	Loop::bSpanAtomVisible( Atom_sp aAtom, bool *bPSeenBefore )
{
  *bPSeenBefore = false;
  if ( aAtom->getSeenId() == this->iSeenId ) {
    *bPSeenBefore = true;
    return(false);
  }

		/* If the ATOM_CLASS in the spanning tree is too far */
		/* from the root ATOM_CLASS then say it is invisible */

  if ( this->iMaxDistanceFromRoot >= 0 ) {
    if ( this->iMaxDistanceFromRoot <
         this->aCurSpan->getBackCount()+1 ) {
      return(false);
    }
  }
  return(this->bLoopAtomVisible( aAtom ));
}





//
// Constructor
//
Loop::Loop()
{
}


Loop::Loop(Matter_sp cont, int over)
{
  this->loopTopGoal(cont,over);
}


void	Loop::clearAtomIndices()
{
  this->iIndex0 = 0;
  this->iIndex1 = 0;
  this->iIndex2 = 0;
  this->iIndex3 = 0;
  this->_Impropers.clear();
  this->_ImproperIndex = -1;
}


void Loop::buildListOfImpropersCenteredOn(Atom_sp a)
{
  BondList_sp		blrest0;
  gctools::Vec0<Bond_sp>::iterator	bl1, bl2, bl3;
  this->_Impropers.clear();
  blrest0 = a->getBondList();
  LOG(BF("Attempting to build a list of impropers centered on %s") % a->description() );
  LOG(BF("It has the bonds: %s") % blrest0->description() );
  if ( blrest0->size() >= 3 )
  {
    LOG(BF("This atom has more than three bonded neighbors"));
    for ( bl1=blrest0->begin(); bl1!=blrest0->end(); bl1++ )
    {
      LOG(BF("Loop 1") );
      GC_ALLOCATE(BondList_O, accumulate1 );
      accumulate1->append(*bl1);
      GC_COPY(BondList_O, blrest1 , *blrest0); // = RP_Copy<BondList_O>(blrest0);
      LOG(BF("Copy blrest1 = %s") % blrest1->description() );
      ASSERTF(blrest1->size() == blrest0->size(),
              BF("The blrest1 copy doesn't have the same number of elements as the original"));
      blrest1->removeBond(*bl1);
      for ( bl2 = blrest1->begin(); bl2!=blrest1->end(); bl2++ )
      {
        LOG(BF("Loop 2") );
        GC_COPY(BondList_O, accumulate2 , *accumulate1); // = RP_Copy<BondList_O>(accumulate1);
        accumulate2->append(*bl2);
        GC_COPY(BondList_O, blrest2 , *blrest1); // = RP_Copy<BondList_O>(blrest1);
        blrest2->removeBond(*bl2);
        for ( bl3=blrest2->begin(); bl3!=blrest2->end(); bl3++ ) {
          LOG(BF("Loop 3") );
          GC_COPY(BondList_O, accumulate3 , *accumulate2); // = RP_Copy<BondList_O>(accumulate2);
          accumulate3->append(*bl3);
          this->_Impropers.push_back(accumulate3);
        }
      }
    }
  }
  LOG(BF("Number of impropers = %d") % this->_Impropers.size() );
}


void	Loop::loopTopGoal( Matter_sp c, int goal )
{

  if ( c.nilp() ) {
    SIMPLE_ERROR(BF("Attempt to loop over null object"));
  }

  this->goal = 	goal;
  this->top = 	c;
  this->initialized = false;
  this->done	= 	false;

                /* Initialize the LOOP memory fields to the default */
                /* of no memory. */

  this->iTempInt = 0;
  this->curSubLoop = -1;
  this->clearAtomIndices();
  this->fVisibilityFlags = 0;

		// Now setup the loop stack
  LOG(BF("Loop::loopTopGoal") );
  if ( !((goal & WAYONLY)==SPANNINGTREE) ) {
    LOG(BF("Loop::loopTopGoal pushing subloop object=%s") % c->description().c_str()  );
    this->pushSubLoop( c );
  }
}


/*!
 * Identify the next object(bond/angle/proper/improper) in the atom
 * and put it in the atom array this->atoms[] and return true
 * If there are no more objects then return false
 */
bool	Loop::nextObjectInAtom()
{
  Atom_sp		top;
  Atom_sp		aAtom1;
  Atom_sp		aAtom2;
  Atom_sp		aAtom3;
  bool		allowDuplicates;
  bool		done;
  void*		topOrder;

  top = (this->getCurSubLoopTop()).as<Atom_O>();
  allowDuplicates = ( (this->goal&ALLOWDUPLICATES) != 0 );
  LOG(BF("Loop::nextObjectInAtom atom name: %s  address=0x%08x") % top->getName().c_str() % top.get() );
  switch ( this->goal&GOALONLY )
  {
	/* When LOOPing over BONDS, use iIndex0 as      */
	/* the bond index for the next Atom to check    */
	/* keep looking until a suitable bond is found  */
	/* Eg: where id(Atom) < id(BONDED_ATOM)         */
	/* and place the two atoms in oaObj[0] and      */
	/* oaObj[1] of the Top Loop                     */
	/* then return the atom we are looping over     */
	/* or NULL if there are no more bonds.          */
	/* The atoms in the bonds can then be retrieved */
	/* Using the LoopGetBond routine                */
  case BONDS:
      LOG(BF("Looking for bond") );
      while ( this->iIndex0 < top->coordination() )
      {
        LOG(BF("status") );
        topOrder = top.get();
        ASSERT(this->iIndex0<top->numberOfBonds());
        Atom_sp bondedTo = top->bondedNeighbor(this->iIndex0);
        Atom_O* bondedOrder = bondedTo.get();
        if ( ( topOrder < bondedOrder ) || allowDuplicates ) {
          this->atoms[0] = top;
          this->atoms[1] = top->bondedNeighbor(this->iIndex0);
          this->bondOrder = top->bondedOrder(this->iIndex0);
          this->_Bond = top->getBonds()[this->iIndex0];
          this->iIndex0++;
          return(true);
        }
        this->iIndex0++;
      }
      break;

	/* When LOOPing over ANGLES, use iIndex0 and    */
	/* iIndex1 as the bond angle indices.           */
	/* An outer Atom is the current Atom.           */
	/* iIndex0 ranges through the current atoms     */
	/* coordination and iIndex1, ranges from 0      */
	/* to the current atoms coordination.           */
	/* When an angle is found place the three atoms */
	/* in order in oaObj, and return the current    */
	/* atom.  Otherwise return NULL                 */
	/* The atoms of the angle are retrieved using   */
	/*the LoopGetAngle routine                      */
  case ANGLES:
      done = false;
      while ( ! done ) {
	    /* First test if all the indices are valid */
        if ( this->iIndex0>=top->coordination() ) goto ANONE;
        aAtom1 = top->bondedNeighbor( this->iIndex0 );
        if ( this->iIndex1>=aAtom1->coordination() ) goto AINC1;
        aAtom2 = aAtom1->bondedNeighbor( this->iIndex1 );
	    /* Test if the atoms are unique */
        if ( top == aAtom2 ) goto AINC1;

	    /* Test if the id's are in the right order */

        if ( !top->atLowerUniqueAtomOrderThan(aAtom2) && !allowDuplicates )
        {
          goto AINC1;
        }

	    /* If it passed all these tests then    */
	    /* it is a valid angle                 */

        done = true;

	    /* Increment everything */
      AINC1:
        this->iIndex1++;
        if ( this->iIndex1 >= aAtom1->coordination() ) {
//AINC0:
          this->iIndex0++;
          this->iIndex1 = 0;
        }
      }
	/* An angle was found, place it in the oaObj array */
      this->atoms[0] = top;
      this->atoms[1] = aAtom1;
      this->atoms[2] = aAtom2;
      return(true);
  ANONE:
      break;




	/* When LOOPing over TORSIONs, use iIndex0,     */
	/* iIndex1, and iIndex2 as indices for the      */
	/* next torsion.                                */
	/* iIndex0 gives the bond index from the current*/
	/* atom to the second atom in the torsion.      */
	/* iIndex1 gives the bond index from the second */
	/* to the third and iIndex2 is from the third   */
	/* to the fourth                                */
	/* Keep searching for torsions until there      */
	/* are no more bonds to search on on the first  */
	/* atom.  Return only torsions where            */
	/* id(first) < id(last) and where every atom in */
	/* the torsion is unique                        */
	/* The search is done using a test/increment    */
	/* sequence.                                    */

  case PROPERS:
      done = false;
      while ( ! done ) {
	    /* First test if all the indices are valid */
        if ( this->iIndex0>=top->coordination() ) goto TNONE;
        aAtom1 = top->bondedNeighbor(this->iIndex0);
        if ( this->iIndex1>=aAtom1->coordination() ) goto TINC1;
        aAtom2 = aAtom1->bondedNeighbor( this->iIndex1 );
        if ( this->iIndex2>=aAtom2->coordination() ) goto TINC2;
        aAtom3 = aAtom2->bondedNeighbor( this->iIndex2 );
	    /* Test if the id's are in the right order */
        if ( !top->atLowerUniqueAtomOrderThan(aAtom3) && !allowDuplicates )  goto TINC2;
	    /* Test if the atoms are unique */
        if ( top == aAtom2 ) goto TINC1;
        if ( top == aAtom3 ) goto TINC2;
        if ( aAtom1 == aAtom3 ) goto TINC2;
	    /* If it passed all these tests then    */
	    /* it is a valid torsion                */
        done = true;
	    /* Increment everything */
      TINC2:
        this->iIndex2++;
        if ( this->iIndex2 >= aAtom2->coordination() ) {
        TINC1:
          this->iIndex1++;
          this->iIndex2 = 0;
          if ( this->iIndex1 >= aAtom1->coordination() ) {
//TINC0:
            this->iIndex0++;
            this->iIndex1 = 0;
          }
        }
      }
	/* A torsion was found, place it in the oaObj array */
      this->atoms[0] = top;
      this->atoms[1] = aAtom1;
      this->atoms[2] = aAtom2;
      this->atoms[3] = aAtom3;
      return(true);
  TNONE:      break;

	/* When LOOPing over IMPROPERS, use iIndex0 and */
	/* iIndex1,iIndex2 as the improper indices.     */
	/* The central atom of the improper is the current */
	/* improper being looped over.                  */
	/* iIndex0 ranges through the current atoms     */
	/* coordination and iIndex1, ranges from iIndex0*/
	/* +1 to the current atoms coordination.        */
	/* iIndex2 ranges form iIndex1+1 to the current */
	/* atoms coordination                           */
	/* When an angle is found place the three atoms */
	/* in order in oaObj, and return the current    */
	/* atom.  Otherwise return NULL                 */
	/* The atoms of the angle are retrieved using   */
	/*the LoopGetImproper routine                   */
	/* The central atom is stored in oaObj[2], the  */
	/* third atom in the list of four.  This is to  */
	/* maintain consistancy with the way impropers  */
	/* are defined in PARM89A.DAT                   */
  case IMPROPERS:
    {_BLOCK_TRACE("Goal is IMPROPERS");
        if ( this->_Impropers.size() == 0 )
        {
          this->buildListOfImpropersCenteredOn(top);
          if ( this->_Impropers.size() == 0 ) goto FAIL;
          this->_ImproperIndex = 0;
        }
        if ( this->_ImproperIndex < this->_Impropers.size() )
        {
          BondList_sp	bl;
          this->atoms[2] = top;
          bl = this->_Impropers[this->_ImproperIndex];
          this->atoms[3] = bl->atIndex(0)->getOtherAtom(top);
          this->atoms[0] = bl->atIndex(1)->getOtherAtom(top);
          this->atoms[1] = bl->atIndex(2)->getOtherAtom(top);
          LOG(BF("Pulled an improper off of the list: %s-%s-%s-%s")
              % this->atoms[0]
              % this->atoms[1]
              % this->atoms[2]
              % this->atoms[3] );
          this->_ImproperIndex++;
          return true;
        }
        goto FAIL;
#if 0
	    /* Anticipate that the indexes all start at 0 */
        if ( this->iIndex1 <= this->iIndex0 )
          this->iIndex1 = this->iIndex0+1;
        if ( this->iIndex2 <= this->iIndex1 )
          this->iIndex2 = this->iIndex1+1;
        else this->iIndex2++;
        if ( this->iIndex2 >= top->coordination() ) {
          this->iIndex1++;
          this->iIndex2 = this->iIndex1 + 1;
          if ( this->iIndex2 >= top->coordination()) {
            this->iIndex0++;
            this->iIndex1 = this->iIndex0 + 1;
            this->iIndex2 = this->iIndex1 + 1;
          }
        }
        if ( this->iIndex2 < (top->coordination()) ) {
          this->atoms[2] = top;
          this->atoms[0] = top->bondedNeighbor(this->iIndex0);
          this->atoms[1] = top->bondedNeighbor(this->iIndex1);
          this->atoms[3] = top->bondedNeighbor(this->iIndex2);
          return(true);
        }
#endif
    }
    break;

  }
 FAIL:
  this->clearAtomIndices();
  return(false);
}



//
//	nextHierarchyMatter
//
//	Private method for determining the next container in the loop
//	that is being iterated through the hierarchy of containres
//
Matter_sp	Loop::nextHierarchyMatter()
{
  Matter_sp	retVal;
  bool		hitEnd;
  int		goalOnly;

  goalOnly = this->goal & GOALONLY;
  hitEnd = false;
  while (1) 
  {
    LOG(BF("running") );
	    // Get the current object in the current subloop.
	    // If the current object is an atom and the goal is
	    //     internal coordinates then do special processing.
	    // Otherwise advance the subloop.
	    // Check if the object satisfies the goal.
	    // If not then create a subloop over the object
    ANN(this->getCurSubLoopTop());
    if ( (this->getCurSubLoopTop()->getMatterType()==ATOM_CLASS) &&
         (this->goal&NONCONTAINERLOOP)!=0 ){
		    // Loop over internal coordinates
      if ( this->nextObjectInAtom() ) {
        LOG(BF("nextObjectInAtom success") );
        retVal = this->getCurSubLoopTop();
        goto DONE;
      } else {
        LOG(BF("nextObjectInAtom failed") );
        hitEnd = true;
      }
    } else {
      LOG(BF("this->curSubLoop=%d") % this->curSubLoop  );
		    // get the next Matter* in the list
		    //
      LOG(BF("running") );
      hitEnd = this->curSubLoopIteratorDone();
      LOG(BF("hitEnd=%d") % hitEnd  );
      if ( !hitEnd ) {
        LOG(BF("running") );
        retVal = this->getCurSubLoopIteratorCur();
        LOG(BF("getCurSubLoopIteratorCur = %s") % retVal->description().c_str()  );
        this->curSubLoopAdvanceIteratorCur();
      }
      LOG(BF("running") );
    }

	    // If we hit the end of the subloop then pop the subloop
	    // and if there are no more subloops then return that
	    // the loop is done
	    // If there are still subloops to be looped over continue
    if (hitEnd) {
      LOG(BF("running") );
      this->popSubLoop();
      LOG(BF("running") );
      if ( this->noSubLoops() ) {
        LOG(BF("returning an empty Matter_sp") );
        retVal = _Nil<Matter_O>();
        this->done = true;
        goto DONE;
      }
      continue;
    }


	    // Now check if the current object satisfies
	    // the goal, if not, create a subloop over it

    LOG(BF("goalOnly = %d") % goalOnly  );
    LOG(BF("retVal->getMatterType() = %d") % retVal->getMatterType()  );
    if ( goalOnly == retVal->getMatterType() ) {
      LOG(BF("goalOnly == retVal->getMatterType() ") );
      if ( goalOnly == (ALL_MATTERS) )
      {
        LOG(BF("goalOnly == ALL_MATTERS") );
        if ( retVal->getMatterType()!=ATOM_CLASS ) 
        {
          LOG(BF("retVal->getMatterType() !=ATOM_CLASS") );
          this->pushSubLoop(retVal);
        }
      }

      if ( retVal->getMatterType() == ATOM_CLASS ) {
		    // Check if the atom is invisible and continue
		    // if it is
      }
      goto DONE;
    }

		    // The current subloop does not satisfy the goals
		    // try another subloop

    LOG(BF("running") );
    this->pushSubLoop(retVal);
  }
 DONE:
  LOG(BF("returning") );
  return retVal;
}


//
// advanceLoop
//
// At the end of this function either this->currentObject is the next object
// or this->done == true, in which case this->currentObject is invalid.
//
void	Loop::advanceLoop()
{
  Matter_sp	retVal;
  int		goalOnly;

  LOG(BF("Loop::next start") );
  if ( this->done ) {
    LOG(BF("Loop done,returning immediatly") );
    return;
  }
  LOG(BF("Loop::next running") );

  goalOnly = this->goal & GOALONLY;

  LOG(BF("Loop::next running") );
  retVal = this->nextHierarchyMatter();
//DONE:
  this->currentObject = retVal;
#ifdef	DEBUG_ON
  ANN(retVal);
  if ( retVal.notnilp() )
  {
    LOG(BF("Loop::next returning with object: %s") % retVal->description().c_str()  );
  }
#endif
  return;
}


#define ARGS_map_molecules "(result-type func matter)"
#define DECL_map_molecules ""
#define DOCS_map_molecules "Loop over molecules, invoke function for each molecule"
CL_DEFUN core::T_sp chem__map_molecules(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m )
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  Loop l(m,MOLECULES);
  Molecule_sp a;
  if ( result_type.nilp() ) {
    while (l.advanceLoopAndProcess()) {
      a = l.getMolecule();
      core::eval::funcall(func,a);
    }
    return _Nil<core::T_O>(); 
  }
  if ( result_type == cl::_sym_list ) {
    ql::list res;
    while (l.advanceLoopAndProcess()) {
      a = l.getMolecule();
      res << core::eval::funcall(func,a);
    }
    return res.cons();
  }
  if ( result_type == cl::_sym_vector ) {
    core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
    while (l.advanceLoopAndProcess()) {
      a = l.getMolecule();
      vo->vectorPushExtend(core::eval::funcall(func,a));
    }
    return vo;
  }
  SIMPLE_ERROR(BF("Illegal return type: %s") % _rep_(result_type));
};

#define ARGS_map_residues "(result-type func matter)"
#define DECL_map_residues ""
#define DOCS_map_residues "Loop over residues, invoke function for each residue"
CL_DEFUN core::T_sp chem__map_residues(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  Loop l(m,RESIDUES);
  Residue_sp a;
  if ( result_type.nilp() ) {
    while (l.advanceLoopAndProcess()) {
      a = l.getResidue();
      core::eval::funcall(func,a);
    }
    return _Nil<core::T_O>(); 
  }
  if ( result_type == cl::_sym_list ) {
    ql::list res;
    while (l.advanceLoopAndProcess()) {
      a = l.getResidue();
      res << core::eval::funcall(func,a);
    }
    return res.cons();
  }
  if ( result_type == cl::_sym_vector ) {
    core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
    while (l.advanceLoopAndProcess()) {
      a = l.getResidue();
      vo->vectorPushExtend(core::eval::funcall(func,a));
    }
    return vo;
  }
  SIMPLE_ERROR(BF("Illegal return type: %s") % _rep_(result_type));
};


#define ARGS_map_atoms "(result-type func matter)"
#define DECL_map_atoms ""
#define DOCS_map_atoms "Loop over atoms, invoke function for each atom"
CL_DEFUN core::T_sp chem__map_atoms(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  Loop l(m,ATOMS);
  Atom_sp a;
  if ( result_type.nilp() ) {
    while (l.advanceLoopAndProcess()) {
      a = l.getAtom();
      core::eval::funcall(func,a);
    }
    return _Nil<core::T_O>(); 
  }
  if ( result_type == cl::_sym_list ) {
    ql::list res;
    while (l.advanceLoopAndProcess()) {
      a = l.getAtom();
      res << core::eval::funcall(func,a);
    }
    return res.cons();
  }
  if ( result_type == cl::_sym_vector ) {
    core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
    while (l.advanceLoopAndProcess()) {
      a = l.getAtom();
      vo->vectorPushExtend(core::eval::funcall(func,a));
    }
    return vo;
  }
  SIMPLE_ERROR(BF("Illegal return type: %s") % _rep_(result_type));
};

#define DOCS_map_atoms "Loop over bonds, invoke function for each bond passing the two atoms and the bond order"
CL_DEFUN core::T_sp chem__map_bonds(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  Loop l(m,BONDS);
  Atom_sp a1, a2;
  BondOrder o;
  if ( result_type.nilp() ) {
    while (l.advanceLoopAndProcess()) {
      a1 = l.getBondA1();
      a2 = l.getBondA2();
      o = l.getBondOrder();
      core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o));
    }
    return _Nil<core::T_O>(); 
  }
  if ( result_type == cl::_sym_list ) {
    ql::list res;
    while (l.advanceLoopAndProcess()) {
      a1 = l.getBondA1();
      a2 = l.getBondA2();
      o = l.getBondOrder();
      res << core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o));
    }
    return res.cons();
  }
  if ( result_type == cl::_sym_vector ) {
    core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
    while (l.advanceLoopAndProcess()) {
      a1 = l.getBondA1();
      a2 = l.getBondA2();
      o = l.getBondOrder();
      vo->vectorPushExtend(core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o)));
    }
    return vo;
  }
  SIMPLE_ERROR(BF("Illegal return type: %s") % _rep_(result_type));
};

void validate_result_type(core::Symbol_sp result_type)
{
  if (result_type.nilp()) return;
  if (result_type == cl::_sym_list) return;
  if (result_type == cl::_sym_vector) return;
  SIMPLE_ERROR(BF("Illegal result type - must be nil, list or vector"));
}

CL_DOCSTRING("Loop over angles, invoke function for each bond passing the three atoms");
CL_DEFUN core::T_sp chem__map_angles(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  validate_result_type(result_type);
  Loop l(m,ANGLES);
  Atom_sp a1, a2, a3;
  ql::list list_result;
  core::ComplexVector_T_sp vector_result;
  if ( result_type == cl::_sym_vector ) {
    vector_result = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
  }
  while (l.advanceLoopAndProcess()) {
    a1 = l.getAtom1();
    a2 = l.getAtom2();
    a3 = l.getAtom3();
    core::T_sp val = core::eval::funcall(func,a1,a2,a3);
    if(result_type.notnilp()) {
      if (result_type == cl::_sym_list) {
        list_result << val;
      } else {
        vector_result->vectorPushExtend(val);
      }
    }
  }
  if(result_type.nilp()) {
    return _Nil<core::T_O>();
  } else if (result_type == cl::_sym_list) {
    return list_result.cons();
  } else {
    return vector_result;
  }
};

CL_DOCSTRING("Loop over dihedrals, invoke function for each bond passing the four atoms");
CL_DEFUN core::T_sp chem__map_dihedrals(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  validate_result_type(result_type);
  Loop l(m,DIHEDRALS);
  Atom_sp a1, a2, a3, a4;
  ql::list list_result;
  core::ComplexVector_T_sp vector_result;
  if ( result_type == cl::_sym_vector ) {
    vector_result = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
  }
  while (l.advanceLoopAndProcess()) {
    a1 = l.getAtom1();
    a2 = l.getAtom2();
    a3 = l.getAtom3();
    a4 = l.getAtom4();
    core::T_sp val = core::eval::funcall(func,a1,a2,a3,a4);
    if(result_type.notnilp()) {
      if (result_type == cl::_sym_list) {
        list_result << val;
      } else {
        vector_result->vectorPushExtend(val);
      }
    }
  }
  if(result_type.nilp()) {
    return _Nil<core::T_O>();
  } else if (result_type == cl::_sym_list) {
    return list_result.cons();
  } else {
    return vector_result;
  }
};

CL_DOCSTRING("Loop over impropers, invoke function for each bond passing the four atoms");
CL_DEFUN core::T_sp chem__map_impropers(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  validate_result_type(result_type);
  Loop l(m,IMPROPERS);
  Atom_sp a1, a2, a3, a4;
  ql::list list_result;
  core::ComplexVector_T_sp vector_result;
  if ( result_type == cl::_sym_vector ) {
    vector_result = core::ComplexVector_T_O::make(16,_Nil<core::T_O>(),core::clasp_make_fixnum(0));
  }
  while (l.advanceLoopAndProcess()) {
    a1 = l.getAtom1();
    a2 = l.getAtom2();
    a3 = l.getAtom3();
    a4 = l.getAtom4();
    core::T_sp val = core::eval::funcall(func,a1,a2,a3,a4);
    if(result_type.notnilp()) {
      if (result_type == cl::_sym_list) {
        list_result << val;
      } else {
        vector_result->vectorPushExtend(val);
      }
    }
  }
  if(result_type.nilp()) {
    return _Nil<core::T_O>();
  } else if (result_type == cl::_sym_list) {
    return list_result.cons();
  } else {
    return vector_result;
  }
};



};
