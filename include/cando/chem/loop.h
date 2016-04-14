/*
    File: loop.h
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
 *	loop.h
 *
 *	Define an object that allows looping over bonds/atoms in
 *	an Atom/Residue/Molecule/Aggregate
 */

#ifndef	LOOP_H
#define	LOOP_H

#include <clasp/core/common.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/molecule.h>
#include <cando/chem/aggregate.h>


#define GOALONLY		0x0000FFFF
#define ALL_MATTERS		(0x000000FF&MATTER_CLASS)
#define	AGGREGATES		(0x000000FF&AGGREGATE_CLASS)
#define MOLECULES		(0x000000FF&MOLECULE_CLASS)
#define RESIDUES		(0x000000FF&RESIDUE_CLASS)
#define ATOMS			(0x000000FF&ATOM_CLASS)
#define DIRECTCONTENTS		0x00000100
#define DIRECTCONTENTSBYSEQNUM	0x00000200

#define NONCONTAINERLOOP        0x00001000
#define BONDS			0x00001001
#define ANGLES			0x00001002
#define PROPERS			0x00001004
#define TORSIONS		PROPERS
#define DIHEDRALS		PROPERS
#define IMPROPERS		0x00001008


#define	TEMPINTVISIBLE		0x00000001
#define	TEMPINTINVISIBLE	0x00000002
#define	TEMPINTUSED		0x00000004

        /* LOOP ways */

                /* ALLOWDUPLICATES allows the loops over BONDS, */
                /* ANGLES, TORSIONS to return duplicates.  Eg:  */
                /* A bond between atoms A-B can be returned from */
                /* both A and B */

#define WAYONLY                 0xFFFF0000
#define DEFAULT                 0x00000000
#define SPANNINGTREE            0x00010000
#define ALLOWDUPLICATES         0x00020000



#define	MAXLOOPS	5

#include <cando/chem/chemPackage.h>


namespace chem {



typedef	enum {
goalOnly = GOALONLY,
allMatters = ALL_MATTERS,
aggregates = AGGREGATES,
molecules = MOLECULES,
residues = RESIDUES,
atoms = ATOMS,
directContents = DIRECTCONTENTS,
bonds = BONDS,
angles = ANGLES,
propers = PROPERS,
torsions = TORSIONS,
dihedrals = DIHEDRALS,
impropers = IMPROPERS
} LoopEnum;

class Loop {
public:
	int				goal;
	bool				initialized;
	bool				done;

	Matter_sp			top;
	Matter_sp			currentObject;
			// Maintain loop stack

  int				curSubLoop;
  Matter_sp			subLoopTop[MAXLOOPS];
  Matter_O::contentIterator	subLoopIteratorCur[MAXLOOPS];
  Matter_O::contentIterator	subLoopIteratorEnd[MAXLOOPS];
			// Maintain indices for bonds/angles/dihedrals/impropers
  int				iIndex0;
  int				iIndex1;
  int				iIndex2;
  int				iIndex3;
  Atom_sp				atoms[4];
  BondOrder			bondOrder;
  Bond_sp				_Bond;

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

  gctools::Vec0<BondList_sp>		_Impropers;
  uint				_ImproperIndex;

  void	clearAtomIndices();
  bool	nextObjectInAtom();

  void	pushSubLoop(Matter_sp c);
  void	popSubLoop()	{this->curSubLoop--;};
  bool	noSubLoops()	{return(this->curSubLoop<0);};
  Matter_sp getCurSubLoopTop()	{return(this->subLoopTop[this->curSubLoop]);};
  Matter_sp getCurSubLoopIteratorCur(){return(*(this->subLoopIteratorCur[this->curSubLoop]));};
  void	curSubLoopAdvanceIteratorCur(){this->subLoopIteratorCur[this->curSubLoop]++;};
  bool	curSubLoopIteratorDone(){return(this->subLoopIteratorCur[this->curSubLoop]>=this->subLoopIteratorEnd[this->curSubLoop]);};

  bool		bLoopAtomVisible(Atom_sp a);
  bool		bSpanAtomVisible(Atom_sp a, bool* b);
  Atom_sp		nextSpanningAtom();
  Matter_sp  nextHierarchyMatter();

 public:

  Loop();
  Loop( Matter_sp cont, int over );

  void	loopTopGoal( Matter_sp cont, int over );

  void 	buildListOfImpropersCenteredOn(Atom_sp a);

  void	loopTopAtomGoal( Atom_sp cont, int goal ) { loopTopGoal( cont, goal ); }
  void	loopTopResidueGoal( Residue_sp cont, int goal ) { loopTopGoal( cont, goal ); }
  void	loopTopMoleculeGoal( Molecule_sp cont, int goal ) { loopTopGoal( cont, goal ); }
  void	loopTopAggregateGoal( Aggregate_sp cont, int goal ) { loopTopGoal( cont, goal ); }

  void		advanceLoop();
			// Return TRUE as long as there is something to return
  bool		advanceLoopAndProcess();
			// Return TRUE as long as there is something to return
  bool		advance() { return this->advanceLoopAndProcess(); };

  Matter_sp	getMatter() { return this->currentObject; }
  Bond_sp		getBond() const { return this->_Bond;};
  Atom_sp		getAtom() { return this->currentObject.as<Atom_O>(); }
  Residue_sp	getResidue() { return this->currentObject.as<Residue_O>(); }
  Molecule_sp	getMolecule() { return this->currentObject.as<Molecule_O>(); }
  void		getBond( Atom_sp *a1, Atom_sp *a2, BondOrder* o ) {*a1=this->atoms[0];*a2=this->atoms[1];*o=this->bondOrder;}
  Atom_sp		getBondA1() { return this->atoms[0]; };
  Atom_sp		getBondA2() { return this->atoms[1]; };
  Atom_sp		getAtom1() { return this->atoms[0]; };
  Atom_sp		getAtom2() { return this->atoms[1]; };
  Atom_sp		getAtom3() { return this->atoms[2]; };
  Atom_sp		getAtom4() { return this->atoms[3]; };
  BondOrder	getBondOrder() { return this->bondOrder; };
  void		getAngle( Atom_sp* a1, Atom_sp* a2, Atom_sp* a3 ) {*a1=this->atoms[0];*a2=this->atoms[1];*a3=this->atoms[2];};
  void		getProper( Atom_sp* a1, Atom_sp* a2, Atom_sp* a3, Atom_sp* a4 ) {*a1=this->atoms[0];*a2=this->atoms[1];*a3=this->atoms[2];*a4=this->atoms[3];};
  void		getImproper( Atom_sp* a1, Atom_sp* a2, Atom_sp* a3, Atom_sp* a4 ) {*a1=this->atoms[0];*a2=this->atoms[1];*a3=this->atoms[2];*a4=this->atoms[3];};


};

};


#endif
