       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	spanningLoop.h
 *
 *	Define an object that allows looping over a spanning
 *	tree of connected atoms.
 */

#ifndef	SPANNINGLOOP_H
#define	SPANNINGLOOP_H

#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/molecule.h>
#include <cando/chem/aggregate.h>


#define	MAXLOOPS	5

#include <cando/chem/chemPackage.h>


namespace chem {


    SMART(SpanningLoop );
    class SpanningLoop_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,SpanningLoop_O,"SpanningLoop");
#if INIT_TO_FACTORIES
    public:
	static SpanningLoop_sp make(Atom_sp root);
#else
	DECLARE_INIT();
#endif
public:
	void initialize();
private:
	int				goal;
	bool				initialized;
	bool				done;
	bool				_OnlyFollowRealBonds;
	Atom_sp				top;
	Atom_sp				currentObject;

			// Spanning tree stuff
	Atom_sp				aCurSpan;
	Atom_sp				aLastSpan;
	int				iSeenId;
	int				iMaxDistanceFromRoot;
	int				iInvisibleCollisions;
	Atom_sp				aLastCollision;
	Atom_sp				aInvisibleAtom;
	AtomFlags			fVisibilityFlags;
	AtomFlags			fVisibleFlagsOn;
	AtomFlags			fVisibleFlagsOff;
	int				iTempInt;


	void	clearAtomIndices();
	bool	nextObjectInAtom();
	bool		bLoopAtomVisible(Atom_sp a);
	bool		bSpanAtomVisible(Atom_sp a, BondOrder order, bool* b);
	Atom_sp		nextSpanningAtom();

public:
	static SpanningLoop_sp create(Atom_sp t) {
	    GC_ALLOCATE(SpanningLoop_O, s );
	    s->setTop(t);
	    return s;
	};
public:

	void		setOnlyFollowRealBonds(bool b) { this->_OnlyFollowRealBonds = b; };
	void		setTop(Atom_sp a);
	void		advanceLoop();
			// Return TRUE as long as there is something to return
	bool		advanceLoopAndProcess();
	bool		advance() { return this->advanceLoopAndProcess(); };
	Atom_sp		getAtom();

	core::List_sp		allAtoms();



	DEFAULT_CTOR_DTOR(SpanningLoop_O);
    };

};

TRANSLATE(chem::SpanningLoop_O);
#endif
