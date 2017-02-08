/*
    File: spanningLoop.h
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
    class SpanningLoop_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,SpanningLoop_O,"SpanningLoop",core::CxxObject_O);
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
        std::function<bool (Bond_sp)>*  _FollowBondIfTrueCallback;
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
        core::T_sp next(core::T_sp funcDesig);
	bool	nextObjectInAtom();
	bool		bLoopAtomVisible(Atom_sp a);
	bool		bSpanAtomVisible(Atom_sp a, BondOrder order, bool* b);
	Atom_sp		nextSpanningAtom(std::function<bool (Atom_sp fromAtom, Bond_sp b)> bondTester);

public:
	static SpanningLoop_sp create(Atom_sp t) {
	    GC_ALLOCATE(SpanningLoop_O, s );
	    s->setTop(t);
	    return s;
	};
public:

        void setFollowBondsWhereTrue(std::function<bool (Bond_sp)>* fn);
	void		setOnlyFollowRealBonds(bool b) { this->_OnlyFollowRealBonds = b; };
	void		setTop(Atom_sp a);
	void		advanceLoop();
	bool		advanceLoopAndProcessWhenTestTrue(std::function<bool (Atom_sp, Bond_sp)>);
			// Return TRUE as long as there is something to return
	bool		advanceLoopAndProcess();
CL_LISPIFY_NAME("advance");
CL_DEFMETHOD 	bool		advance() { return this->advanceLoopAndProcess(); };
	Atom_sp		getAtom();

	core::List_sp		allAtoms();



	DEFAULT_CTOR_DTOR(SpanningLoop_O);
    };

};

TRANSLATE(chem::SpanningLoop_O);
#endif
