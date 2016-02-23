/*
    File: spanningLoop.cc
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
#define	DEBUG_LEVEL_FULL
//
// (C) 2004 Christian E. Schafmeister
//



#include <clasp/core/common.h>

#include <stdio.h>
#include <clasp/core/designators.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/loop.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/matter.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <clasp/core/environment.h>
#include <cando/chem/atom.h>
#include <cando/chem/bond.h>
#include <clasp/core/wrappers.h>

namespace chem
{
    static		int	SiNextSeenId = 0;

#define ARGS_SpanningLoop_O_make "(file_name)"
#define DECL_SpanningLoop_O_make ""
#define DOCS_SpanningLoop_O_make "make SpanningLoop"
CL_NAME(CHEM:MAKE-SPANNING-LOOP);
CL_DEFUN SpanningLoop_sp SpanningLoop_O::make(Atom_sp root)
  {
      GC_ALLOCATE(SpanningLoop_O, me );
      me->setTop(root);
      return me;
  };
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
    bool	SpanningLoop_O:: bLoopAtomVisible(Atom_sp aAtom )
    {
	AtomFlags		fFlags;
	LOG(BF("status") );
	fFlags = aAtom->getFlags();
	if ( (fFlags&this->fVisibleFlagsOn)!=this->fVisibleFlagsOn )
	    return(false);

	LOG(BF("status") );
	if ( (fFlags|(~(this->fVisibleFlagsOff)))!=(~(this->fVisibleFlagsOff)) )
	    return(false);

	LOG(BF("status") );
	/* If the TempInt field is being used to determine */
	/* visibility, then check if the ATOM_CLASS has the */
	/* proper TempInt field */

	LOG(BF("status") );
	if ( this->fVisibilityFlags & TEMPINTUSED ) {
	    if ( this->fVisibilityFlags & TEMPINTINVISIBLE ) {
		if ( aAtom->getTempInt() == this->iTempInt ) {
		    LOG(BF("status") );
		    return(false);
		}
	    } else {
		if ( aAtom->getTempInt() != this->iTempInt ) {
		    LOG(BF("status") );
		    return(false);
		}
	    }
	}
	LOG(BF("status") );

	/* Check if the ATOM_CLASS is the one invisible ATOM_CLASS for the loop */
	/* Used to construct spanning trees that do not pass a */
	/* certian ATOM_CLASS, good for looping over side chains */

        if ( this->aInvisibleAtom == aAtom ) return false;
#if 0        
	if ( this->aInvisibleAtom.otherp() )
	{
	    if ( this->aInvisibleAtom == aAtom ) 
	    {
		return(false);
	    }
	}
#endif
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
    bool	SpanningLoop_O::bSpanAtomVisible( Atom_sp aAtom, BondOrder order, bool *bPSeenBefore )
    {
	*bPSeenBefore = false;
	if ( this->_OnlyFollowRealBonds )
	{
	    if ( order != singleBond 
		 && order != doubleBond 
		 && order != tripleBond 
		) return false;
	}
	if ( aAtom->getSeenId() == this->iSeenId ) {
	    *bPSeenBefore = true;
	    return(false);
	}

	/* If the ATOM_CLASS in the spanning tree is too far */
	/* from the root ATOM_CLASS then say it is invisible */

	if ( this->iMaxDistanceFromRoot >=0 ) {
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




CL_LISPIFY_NAME("setTop");
CL_DEFMETHOD     void SpanningLoop_O::setTop( Atom_sp c )
    {
	this->top = 	c;
	this->initialized = false;
	this->done	= 	false;

	/* Initialize the LOOP memory fields to the default */
	/* of no memory. */

	this->iTempInt = 0;
	this->fVisibilityFlags = 0;

	// Now setup the loop stack


	LOG(BF("SpanningLoop_O::setTop atom moeIndex =%d") % (c->getMoeIndex() ) );
	this->iMaxDistanceFromRoot = -1;
	this->aCurSpan = c;
	this->aLastSpan = c;
	SiNextSeenId++;
	if ( SiNextSeenId == 0 ) SiNextSeenId = 1;
	this->iSeenId = SiNextSeenId;
	this->iInvisibleCollisions = 0;
	this->aLastCollision = Atom_sp();
	this->fVisibilityFlags = 0;
	this->fVisibleFlagsOn = 0;
	this->fVisibleFlagsOff = 0;
    }



CL_LISPIFY_NAME("getAtom");
CL_DEFMETHOD     Atom_sp	SpanningLoop_O::getAtom()
    {
	return downcast<Atom_O>(this->currentObject);
    }

CL_LISPIFY_NAME("advanceLoopAndProcess");
CL_DEFMETHOD     bool	SpanningLoop_O::advanceLoopAndProcess()
    {
	if ( !this->done) {
	    this->advanceLoop();
	    return true;
	} else {
	    return false;
	}
    }

    void SpanningLoop_O::initialize()
    {
	this->Base::initialize();
	this->initialized = false;
	this->_OnlyFollowRealBonds = false;
    }


//
//	nextSpanningAtom
//
//	Private method for advancing to the next atom in a
//	spanning tree.
Atom_sp	SpanningLoop_O::nextSpanningAtom(std::function<bool (Atom_sp fromAtom, Bond_sp b)> bondTester )
    {
	Atom_sp		oObject;
	Atom_sp		aBond;
	Atom_sp		aPrev;
	int		i;
	bool		bSeenBefore;

	/* The first time the LOOP is entered, set up */
	/* some stuff on the first ATOM_CLASS */

	if ( !this->initialized ) {
          ASSERT(this->top);
          gc::As<Atom_sp>(this->top)->invalidateBackSpan();
          gc::As<Atom_sp>(this->top)->setBackCount(0);
          gc::As<Atom_sp>(this->top)->setSeenId(this->iSeenId);
          this->initialized = true;
	}

	LOG(BF("--- Entered nextSpanningAtom aCurSpan = %d; nextSpan valid?=%d") 
	    % this->aCurSpan->getMoeIndex() % this->aCurSpan->isNextSpanValid() );
	ASSERT(this->aCurSpan.notnilp());
//    if (IsNull(this->aCurSpan)) {
//	LOG(BF("this->aCurSpan is NULL") );
//	oObject = Atom();	// Return a dummy atom
//	this->done = true;
//	goto DONE;
	//   }

	/* Connect all visible atoms bonded to the current atom */
	/* into the spanning tree */

	aPrev = this->aLastSpan;
	aPrev->invalidateNextSpan();
	LOG(BF("--- Invalidated nextSpan for atom: %d") % (aPrev->getMoeIndex() ) );
	LOG(BF("--- Currently on atom: %d") % (this->aCurSpan->getMoeIndex() ) );
//        printf("%s:%d nextSpanningAtom: %s coordination: %d\n", __FILE__, __LINE__, _rep_(this->aCurSpan).c_str(), this->aCurSpan->coordination());
	for ( i=0; i<this->aCurSpan->coordination(); i++ ) {
          Bond_sp bond = this->aCurSpan->bondAtIndex(i);
          aBond = bond->getOtherAtom(this->aCurSpan);
          BondOrder order = bond->getOrder();
          bool followBond = bondTester(this->aCurSpan,bond);
	    /* If the atom is visible then add it */
          if ( followBond && this->bSpanAtomVisible( aBond, order, &bSeenBefore ) ) {
            LOG(BF("--- looking at bonded atom: %d") % (aBond->getMoeIndex() ) );
            aBond->setSeenId( this->iSeenId );
            aBond->setBackCount( this->aCurSpan->getBackCount()+1 );
            aBond->setBackSpan( this->aCurSpan );
            aPrev->setNextSpan( aBond );
            LOG(BF("--- Setting atom: %d nextSpan to: %d") % (aPrev->getMoeIndex()) % (aBond->getMoeIndex() ) );

            aBond->invalidateNextSpan();
            LOG(BF("--- Is atom: %d nextSpan valid? ==> %d") % (aPrev->getMoeIndex()) % (aPrev->isNextSpanValid() ) );
            aPrev = aBond;
          } else {
		LOG(BF("--- NOT Visible") );

		/* If the atom is invisible, but has not been seen */
		/* before then it counts as an invisible collision */
		/* Increment the collision count and save the atom */
		/* which caused the collision                      */
		if ( !bSeenBefore ) {
		    LOG(BF("--- COLLISION!") );
		    this->iInvisibleCollisions++;
		    this->aLastCollision = aBond;
		}
	    }
	}

	/* Advance to the next atom and return the current one */
	LOG(BF("Advancing spanning tree") );
	this->aLastSpan = aPrev;
	oObject = this->aCurSpan;
	if ( this->aCurSpan->isNextSpanValid() ) {
	    this->aCurSpan = this->aCurSpan->getNextSpan();
	    LOG(BF("--- Setting this->aCurSpan to: %d") % (this->aCurSpan->getMoeIndex() ) );
	    LOG(BF("--- Is atom: %d nextSpan valid? ==> %d") % (this->aCurSpan->getMoeIndex()) % (this->aCurSpan->isNextSpanValid() ) );
	} else {
	    LOG(BF("atom %d does not have a next span atom, loop done") % (this->aCurSpan->getMoeIndex() ) );
	    this->done = true;
	    goto DONE;
	}

	LOG(BF("About to return") );

	/* Return with oObject */
    DONE:
	return(oObject);
    }



//
// advanceLoop
//
// At the end of this function either this->currentObject is the next object
// or this->done == true, in which case this->currentObject is invalid.
//
void	SpanningLoop_O::advanceLoop()
{
  if ( this->done ) return;
  Atom_sp retVal = this->nextSpanningAtom( [] (Atom_sp a, Bond_sp b) -> bool { return true; });
  this->currentObject = retVal;
  return;
}

bool SpanningLoop_O::advanceLoopAndProcessWhenTestTrue(std::function<bool (Atom_sp, Bond_sp)> followBondIfTrue )
{
  if ( this->done ) return false;
  Atom_sp retVal = this->nextSpanningAtom(followBondIfTrue);
  this->currentObject = retVal;
  return true;
}


CL_LISPIFY_NAME("next");
CL_DEFMETHOD core::T_sp SpanningLoop_O::next(core::T_sp funcDesig) {
  if ( this->done ) return _Nil<core::T_O>();
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  return this->nextSpanningAtom([&func] (Atom_sp a, Bond_sp b) -> bool {
      core::T_sp res = core::eval::funcall(func,a,b);
      return res.notnilp();
    } );
}


CL_LISPIFY_NAME("allAtoms");
CL_DEFMETHOD     core::List_sp	SpanningLoop_O::allAtoms()
    {
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
	core::Cons_sp cur = first;
	Atom_sp theTop = this->top;
	this->setTop(theTop);
	while ( this->advanceLoopAndProcess() )
	{
	    core::Cons_sp one = core::Cons_O::create(this->getAtom(),_Nil<core::T_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	return first->cdr();
    }








};
