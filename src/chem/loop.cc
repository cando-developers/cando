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
  LOG("pushSubLoop: Starting a subloop over: {}" , c->description().c_str()  );
  this->curSubLoop++;
  this->subLoopTop[this->curSubLoop] = c;
  this->subLoopIteratorCur[this->curSubLoop] = c->_Contents.begin();
  this->subLoopIteratorEnd[this->curSubLoop] = c->_Contents.end();
  LOG("pushSubLoop:  number of children count = {}" , c->_Contents.size()  );
}


bool	Loop::advanceLoopAndProcess()
{
  if ( !this->done) {
    LOG("Loop may still have values to generate" );
    this->advanceLoop();
    if ( this->done )
    {
      LOG("The loop is done,returning false" );
      return false;
    }
    LOG("The loop continues,returning true" );
    return true;
  }
  LOG("The loop was done in the last step,returning false" );
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
    SIMPLE_ERROR("THis code won't work with new way to handle atoms");
#if 0
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
#endif
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


void	Loop::clearAtomIndexes()
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
  LOG("Attempting to build a list of impropers centered on {}" , a->description() );
  LOG("It has the bonds: {}" , blrest0->description() );
  if ( blrest0->size() >= 3 )
  {
    LOG("This atom has more than three bonded neighbors");
    for ( bl1=blrest0->begin(); bl1!=blrest0->end(); bl1++ )
    {
      LOG("Loop 1" );
      auto  accumulate1  = gctools::GC<BondList_O>::allocate_with_default_constructor();
      accumulate1->append(*bl1);
      BondList_sp blrest1 = blrest0->deepishCopy(); 
      LOG("Copy blrest1 = {}" , blrest1->description() );
      ASSERTF(blrest1->size() == blrest0->size(),
              ("The blrest1 copy doesn't have the same number of elements as the original"));
      blrest1->removeBond(*bl1);
      for ( bl2 = blrest1->begin(); bl2!=blrest1->end(); bl2++ )
      {
        LOG("Loop 2" );
        BondList_sp accumulate2 = accumulate1->deepishCopy();
        accumulate2->append(*bl2);
        BondList_sp blrest2 = blrest1->deepishCopy();
        blrest2->removeBond(*bl2);
        for ( bl3=blrest2->begin(); bl3!=blrest2->end(); bl3++ ) {
          LOG("Loop 3" );
          BondList_sp accumulate3 = accumulate2->deepishCopy();
          accumulate3->append(*bl3);
          this->_Impropers.push_back(accumulate3);
        }
      }
    }
  }
  LOG("Number of impropers = {}" , this->_Impropers.size() );
}


void	Loop::loopTopGoal( Matter_sp c, int goal )
{

  if ( c.nilp() ) {
    SIMPLE_ERROR("Attempt to loop over null object");
  }

  this->goal = 	goal;
  this->top = 	c;
  this->initialized = false;
  this->done	= 	false;

                /* Initialize the LOOP memory fields to the default */
                /* of no memory. */

  this->iTempInt = 0;
  this->curSubLoop = -1;
  this->clearAtomIndexes();
  this->fVisibilityFlags = 0;

		// Now setup the loop stack
  LOG("Loop::loopTopGoal" );
  if ( !((goal & WAYONLY)==SPANNINGTREE) ) {
    LOG("Loop::loopTopGoal pushing subloop object={}" , c->description().c_str()  );
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
  LOG("Loop::nextObjectInAtom atom name: {}  address=0x%08x" , top->getName().c_str() , top.get() );
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
      LOG("Looking for bond" );
      while ( this->iIndex0 < top->coordination() )
      {
        LOG("status" );
        topOrder = &*top;
        ASSERT(this->iIndex0<top->numberOfBonds());
        Atom_sp bondedTo = top->bondedNeighbor(this->iIndex0);
        Atom_O* bondedOrder = &*bondedTo;
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
	/* iIndex1 as the bond angle indexes.           */
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
	    /* First test if all the indexes are valid */
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
	/* iIndex1, and iIndex2 as indexes for the      */
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
	    /* First test if all the indexes are valid */
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
	/* iIndex1,iIndex2 as the improper indexes.     */
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
    {
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
          LOG("Pulled an improper off of the list: {}-{}-{}-{}"
              , this->atoms[0]
              , this->atoms[1]
              , this->atoms[2]
              , this->atoms[3] );
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
  this->clearAtomIndexes();
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
    LOG("running" );
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
        LOG("nextObjectInAtom success" );
        retVal = this->getCurSubLoopTop();
        goto DONE;
      } else {
        LOG("nextObjectInAtom failed" );
        hitEnd = true;
      }
    } else {
      LOG("this->curSubLoop={}" , this->curSubLoop  );
		    // get the next Matter* in the list
		    //
      LOG("running" );
      hitEnd = this->curSubLoopIteratorDone();
      LOG("hitEnd={}" , hitEnd  );
      if ( !hitEnd ) {
        LOG("running" );
        retVal = this->getCurSubLoopIteratorCur();
        LOG("getCurSubLoopIteratorCur = {}" , retVal->description().c_str()  );
        this->curSubLoopAdvanceIteratorCur();
      }
      LOG("running" );
    }

	    // If we hit the end of the subloop then pop the subloop
	    // and if there are no more subloops then return that
	    // the loop is done
	    // If there are still subloops to be looped over continue
    if (hitEnd) {
      LOG("running" );
      this->popSubLoop();
      LOG("running" );
      if ( this->noSubLoops() ) {
        LOG("returning an empty Matter_sp" );
        retVal = nil<Matter_O>();
        this->done = true;
        goto DONE;
      }
      continue;
    }


	    // Now check if the current object satisfies
	    // the goal, if not, create a subloop over it

    LOG("goalOnly = {}" , goalOnly  );
    LOG("retVal->getMatterType() = {}" , retVal->getMatterType()  );
    if ( goalOnly == retVal->getMatterType() ) {
      LOG("goalOnly == retVal->getMatterType() " );
      if ( goalOnly == (ALL_MATTERS) )
      {
        LOG("goalOnly == ALL_MATTERS" );
        if ( retVal->getMatterType()!=ATOM_CLASS ) 
        {
          LOG("retVal->getMatterType() !=ATOM_CLASS" );
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

    LOG("running" );
    this->pushSubLoop(retVal);
  }
 DONE:
  LOG("returning" );
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

  LOG("Loop::next start" );
  if ( this->done ) {
    LOG("Loop done,returning immediatly" );
    return;
  }
  LOG("Loop::next running" );

  goalOnly = this->goal & GOALONLY;

  LOG("Loop::next running" );
  retVal = this->nextHierarchyMatter();
//DONE:
  this->currentObject = retVal;
#ifdef	DEBUG_ON
  ANN(retVal);
  if ( retVal.notnilp() )
  {
    LOG("Loop::next returning with object: {}" , retVal->description().c_str()  );
  }
#endif
  return;
}




CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive molecules in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_molecules(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m )
{
  if (gc::IsA<Molecule_sp>(m)) {
    core::Function_sp func = core::coerce::functionDesignator(funcDesig);
    if ( result_type.nilp() ) {
      core::eval::funcall(func,m);
      return nil<core::T_O>(); 
    } else if ( result_type == cl::_sym_list ) {
      ql::list res;
      res << core::eval::funcall(func,m);
      return res.cons();
    } else if ( result_type == cl::_sym_vector ) {
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      vo->vectorPushExtend(core::eval::funcall(func,m));
      return vo;
    }
    SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
  } else if (gc::IsA<Aggregate_sp>(m)) {
    core::Function_sp func = core::coerce::functionDesignator(funcDesig);
    Loop l(m,MOLECULES);
    Molecule_sp a;
    if ( result_type.nilp() ) {
      while (l.advanceLoopAndProcess()) {
        a = l.getMolecule();
        core::eval::funcall(func,a);
      }
      return nil<core::T_O>(); 
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
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      while (l.advanceLoopAndProcess()) {
        a = l.getMolecule();
        vo->vectorPushExtend(core::eval::funcall(func,a));
      }
      return vo;
    }
    SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
  }
  TYPE_ERROR(m,_sym_Aggregate_O);
};

CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive residues in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_residues(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  if (gc::IsA<Residue_sp>(m)) {
    core::Function_sp func = core::coerce::functionDesignator(funcDesig);
    if ( result_type.nilp() ) {
      core::eval::funcall(func,m);
      return nil<core::T_O>(); 
    } else if ( result_type == cl::_sym_list ) {
      ql::list res;
      res << core::eval::funcall(func,m);
      return res.cons();
    } else if ( result_type == cl::_sym_vector ) {
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      vo->vectorPushExtend(core::eval::funcall(func,m));
      return vo;
    }
    SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
  } else if (gc::IsA<Aggregate_sp>(m)||gc::IsA<Molecule_sp>(m)) {
    core::Function_sp func = core::coerce::functionDesignator(funcDesig);
    Loop l(m,RESIDUES);
    Residue_sp a;
    if ( result_type.nilp() ) {
      while (l.advanceLoopAndProcess()) {
        a = l.getResidue();
        core::eval::funcall(func,a);
      }
      return nil<core::T_O>(); 
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
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      while (l.advanceLoopAndProcess()) {
        a = l.getResidue();
        vo->vectorPushExtend(core::eval::funcall(func,a));
      }
      return vo;
    }
    SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
  }
  TYPE_ERROR(m,core::Cons_O::createList(cl::_sym_or,_sym_Aggregate_O,_sym_Molecule_O));
};


CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive atoms in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_atoms(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  if (gc::IsA<Atom_sp>(m)) {
    if ( result_type.nilp() ) {
      core::eval::funcall(func,m);
      return nil<core::T_O>(); 
    } else if ( result_type == cl::_sym_list ) {
      ql::list res;
      res << core::eval::funcall(func,m);
      return res.cons();
    } else if ( result_type == cl::_sym_vector ) {
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      vo->vectorPushExtend(core::eval::funcall(func,m));
      return vo;
    }
    SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
  } else if ((gc::IsA<Aggregate_sp>(m)||gc::IsA<Molecule_sp>(m)||gc::IsA<Residue_sp>(m))) {
    Loop l(m,ATOMS);
    Atom_sp a;
    if ( result_type.nilp() ) {
      while (l.advanceLoopAndProcess()) {
        a = l.getAtom();
        core::eval::funcall(func,a);
      }
      return nil<core::T_O>(); 
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
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      while (l.advanceLoopAndProcess()) {
        a = l.getAtom();
        vo->vectorPushExtend(core::eval::funcall(func,a));
      }
      return vo;
    }
    SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
  };
  TYPE_ERROR(m,core::Cons_O::createList(cl::_sym_or,_sym_Aggregate_O,_sym_Molecule_O,_sym_Residue_O));
}

CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive bonds in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_bonds(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  if (gc::IsA<Atom_sp>(m)) {
    Atom_sp am = gc::As_unsafe<Atom_sp>(m);
    if ( result_type.nilp() ) {
      for ( auto b=am->_Bonds.begin();b!=am->_Bonds.end() ; b++ ) {
        core::eval::funcall( func, (*b)->getAtom1(), (*b)->getAtom2(),
                             translate::to_object<BondOrder>::convert((*b)->_DirectionalOrder), (*b));
      }
      return nil<core::T_O>();
    } else if ( result_type == cl::_sym_list ) {
      ql::list l;
      for ( auto b=am->_Bonds.begin();b!=am->_Bonds.end() ; b++ ) {
        l << core::eval::funcall( func, (*b)->getAtom1(), (*b)->getAtom2(),
                                  translate::to_object<BondOrder>::convert((*b)->_DirectionalOrder), (*b));
      }
      return l.cons();
    } else if ( result_type == cl::_sym_array ) {
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      for ( auto b=am->_Bonds.begin();b!=am->_Bonds.end() ; b++ ) {
        core::T_sp obj = core::eval::funcall( func, (*b)->getAtom1(), (*b)->getAtom2(),
                                              translate::to_object<BondOrder>::convert((*b)->_DirectionalOrder), (*b));
        vo->vectorPushExtend(obj);
      }
      return vo;
    }
  } else {
    Loop l(m,BONDS);
    Atom_sp a1, a2;
    BondOrder o;
    if ( result_type.nilp() ) {
      while (l.advanceLoopAndProcess()) {
        Bond_sp bond = l.getBond();
        a1 = bond->getAtom1();
        a2 = bond->getAtom2();
        o = bond->getOrderFromAtom(a1);
        core::eval::funcall( func, a1, a2, translate::to_object<BondOrder>::convert(o), bond);
      }
      return nil<core::T_O>(); 
    }
    if ( result_type == cl::_sym_list ) {
      ql::list res;
      while (l.advanceLoopAndProcess()) {
        Bond_sp bond = l.getBond();
        a1 = bond->getAtom1();
        a2 = bond->getAtom2();
        o = bond->getOrderFromAtom(a1);
        res << core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o), bond );
      }
      return res.cons();
    }
    if ( result_type == cl::_sym_vector ) {
      core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
      while (l.advanceLoopAndProcess()) {
        Bond_sp bond = l.getBond();
        a1 = bond->getAtom1();
        a2 = bond->getAtom2();
        o = bond->getOrderFromAtom(a1);
        vo->vectorPushExtend(core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o), bond));
      }
      return vo;
    }
  }
  SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
};

CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive bonds in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_bond_objects(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  Loop l(m,BONDS);
  Atom_sp a1, a2;
  BondOrder o;
  if ( result_type.nilp() ) {
    while (l.advanceLoopAndProcess()) {
      Bond_sp bond = l.getBond();
      a1 = bond->getAtom1();
      a2 = bond->getAtom2();
      o = bond->getOrderFromAtom(a1);
      core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o));
    }
    return nil<core::T_O>(); 
  }
  if ( result_type == cl::_sym_list ) {
    ql::list res;
    while (l.advanceLoopAndProcess()) {
      Bond_sp bond = l.getBond();
      a1 = bond->getAtom1();
      a2 = bond->getAtom2();
      o = bond->getOrderFromAtom(a1);
      res << core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o));
    }
    return res.cons();
  }
  if ( result_type == cl::_sym_vector ) {
    core::ComplexVector_T_sp vo = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
    while (l.advanceLoopAndProcess()) {
      Bond_sp bond = l.getBond();
      a1 = bond->getAtom1();
      a2 = bond->getAtom2();
      o = bond->getOrderFromAtom(a1);
      vo->vectorPushExtend(core::eval::funcall(func,a1,a2,translate::to_object<BondOrder>::convert(o)));
    }
    return vo;
  }
  SIMPLE_ERROR("Illegal return type: {}" , _rep_(result_type));
};

void validate_result_type(core::Symbol_sp result_type)
{
  if (result_type.nilp()) return;
  if (result_type == cl::_sym_list) return;
  if (result_type == cl::_sym_vector) return;
  SIMPLE_ERROR("Illegal result type - must be nil, list or vector");
}

CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive angles in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_angles(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  validate_result_type(result_type);
  Loop l(m,ANGLES);
  Atom_sp a1, a2, a3;
  ql::list list_result;
  core::ComplexVector_T_sp vector_result;
  if ( result_type == cl::_sym_vector ) {
    vector_result = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
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
    return nil<core::T_O>();
  } else if (result_type == cl::_sym_list) {
    return list_result.cons();
  } else {
    return vector_result;
  }
};

CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive dihedrals in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_dihedrals(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  validate_result_type(result_type);
  Loop l(m,DIHEDRALS);
  Atom_sp a1, a2, a3, a4;
  ql::list list_result;
  core::ComplexVector_T_sp vector_result;
  if ( result_type == cl::_sym_vector ) {
    vector_result = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
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
    return nil<core::T_O>();
  } else if (result_type == cl::_sym_list) {
    return list_result.cons();
  } else {
    return vector_result;
  }
};

CL_LAMBDA(result-type function matter);
CL_DOCSTRING(R"dx(Applies function to successive impropers in matter. The result-type (nil, 'list, 'vector) specifies the type of the resulting sequence.)dx");
DOCGROUP(cando);
CL_DEFUN core::T_sp chem__map_impropers(core::Symbol_sp result_type, core::T_sp funcDesig, Matter_sp m)
{
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  validate_result_type(result_type);
  Loop l(m,IMPROPERS);
  Atom_sp a1, a2, a3, a4;
  ql::list list_result;
  core::ComplexVector_T_sp vector_result;
  if ( result_type == cl::_sym_vector ) {
    vector_result = core::ComplexVector_T_O::make(16,nil<core::T_O>(),core::clasp_make_fixnum(0));
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
    return nil<core::T_O>();
  } else if (result_type == cl::_sym_list) {
    return list_result.cons();
  } else {
    return vector_result;
  }
};



};
