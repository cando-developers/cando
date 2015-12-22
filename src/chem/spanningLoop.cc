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

#if INIT_TO_FACTORIES

#define ARGS_SpanningLoop_O_make "(file_name)"
#define DECL_SpanningLoop_O_make ""
#define DOCS_SpanningLoop_O_make "make SpanningLoop"
  SpanningLoop_sp SpanningLoop_O::make(Atom_sp root)
  {_G();
      GC_ALLOCATE(SpanningLoop_O, me );
      me->setTop(root);
      return me;
  };

#else

    core::T_sp SpanningLoop_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	Atom_sp root = translate::from_object<Atom_O>::convert(env->lookup(ChemPkg,"root"));
	this->setTop(root);
	return _Nil<core::T_O>();
    }

#endif

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
    {_G();
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
    {_G();
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




CL_NAME("setTop");
CL_DEFMETHOD     void SpanningLoop_O::setTop( Atom_sp c )
    {_G();
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



CL_NAME("getAtom");
CL_DEFMETHOD     Atom_sp	SpanningLoop_O::getAtom()
    {_G();
	return downcast<Atom_O>(this->currentObject);
    }

CL_NAME("advanceLoopAndProcess");
CL_DEFMETHOD     bool	SpanningLoop_O::advanceLoopAndProcess()
    {_G();
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
    {_G();
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


CL_NAME("next");
CL_DEFMETHOD core::T_sp SpanningLoop_O::next(core::T_sp funcDesig) {
  if ( this->done ) return _Nil<core::T_O>();
  core::Function_sp func = core::coerce::functionDesignator(funcDesig);
  return this->nextSpanningAtom([&func] (Atom_sp a, Bond_sp b) -> bool {
      core::T_sp res = core::eval::funcall(func,a,b);
      return res.notnilp();
    } );
}


CL_NAME("allAtoms");
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


void SpanningLoop_O::exposeCando(core::Lisp_sp lisp)
{
  core::class_<SpanningLoop_O>()
//	    .def_raw("core:__init__",&SpanningLoop_O::__init__,"(self &key root)")
    .def("setTop",&SpanningLoop_O::setTop)
    .def("next",&SpanningLoop_O::next)
    .def("advance",&SpanningLoop_O::advance)
    .def("advanceLoopAndProcess",&SpanningLoop_O::advanceLoopAndProcess)
    .def("getAtom",&SpanningLoop_O::getAtom)
    .def("allAtoms",&SpanningLoop_O::allAtoms)
    ;
}

    void SpanningLoop_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,SpanningLoop,"","",_lisp)
	    .def("setTop",&SpanningLoop_O::setTop)
	    .def("advance",&SpanningLoop_O::advance)
	    .def("advanceLoopAndProcess",&SpanningLoop_O::advanceLoopAndProcess)
	    .def("getAtom",&SpanningLoop_O::getAtom)
	    .def("allAtoms",&SpanningLoop_O::allAtoms)
	    ;
#endif
    }



    EXPOSE_CLASS(chem,SpanningLoop_O);


};
