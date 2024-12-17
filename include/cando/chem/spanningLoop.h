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
#include <clasp/core/hashTable.fwd.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/molecule.h>
#include <cando/chem/aggregate.h>


#define	MAXLOOPS	5

#include <cando/chem/chemPackage.h>


namespace chem {

    SMART(SpanningInfo );
    class SpanningInfo_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,SpanningInfo_O,"SpanningInfo",core::CxxObject_O);
    public:
        int        _Distance;
        core::T_sp _ToRoot;
        core::T_sp _Next;
        SpanningInfo_O(int dist=0, core::T_sp toRoot=nil<core::T_O>(), core::T_sp toNext=nil<core::T_O>() ) : _Distance(dist), _ToRoot(toRoot), _Next(toNext) {};
        std::string __repr__() const;
    };


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
	Atom_sp			aCurSpan;
        Atom_sp  		aLastSpan;
        SpanningInfo_sp         siLastSpan;
	int				iSeenId;
	int				iMaxDistanceFromRoot;
	int				iInvisibleCollisions;
	Atom_sp				aLastCollision;
	Atom_sp				aInvisibleAtom;
	AtomFlags			fVisibilityFlags;
	AtomFlags			fVisibleFlagsOn;
	AtomFlags			fVisibleFlagsOff;
	int				iTempInt;
        core::HashTable_sp            _BackSpan;

	void	clearAtomIndexes();
        core::T_sp next(core::T_sp funcDesig);
	bool	nextObjectInAtom();
	bool		bLoopAtomVisible(Atom_sp a);
	bool		bSpanAtomVisible(Atom_sp a, BondOrder order, bool* b);
	Atom_sp		nextSpanningAtom(std::function<bool (Atom_sp fromAtom, Bond_sp b)> bondTester);

    public:
        bool isBackSpanValid(Atom_sp a);
        core::T_sp getBackSpan(Atom_sp a);
        bool isNextSpanValid(Atom_sp a);
        core::T_sp getNextSpan(Atom_sp a);
        int getBackCount(Atom_sp a);
        bool lookupSpanningInfo(Atom_sp a, SpanningInfo_sp& info);
        SpanningInfo_sp getSpanningInfo(Atom_sp a);
        SpanningInfo_sp storeSpanningInfo(Atom_sp key, int distance, core::T_sp toRoot=nil<core::T_O>(), core::T_sp next=nil<core::T_O>());

public:
	static SpanningLoop_sp create(Atom_sp t) {
          auto  s  = gctools::GC<SpanningLoop_O>::allocate_with_default_constructor();
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

#endif
