/*
    File: delayedBondedAtom.cc
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
#define DEBUG_LEVEL_FULL



#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/symbolTable.h>
#include <cando/kinematics/rootAtomInfo.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/delayedBondedAtom.h>

namespace kinematics
{

    RefCountedAtomHandle DelayedBondedAtom::stubAtom3(AtomTree_sp tree) const
    {
	if ( this->_DelayType == delayForFollowingResidueBond1 )
	{
//	    chem::ConstitutionAtomIndex0N id = this->_DelayAtomId;
	    chem::AtomId bond0Id(this->_Id.moleculeId(),this->_Id.residueId(),this->_DelayAtomId);
	    RefCountedAtomHandle outBond0 = tree->lookup(bond0Id);
	    int idx = outBond0.get()->firstNonJumpChildIndex();
	    RefCountedAtomHandle nextResidueBond0 = outBond0.get()->child(idx);
	    chem::AtomId nextResidueBond0AtomId = nextResidueBond0.get()->id();
	    ASSERTF(bond0Id.residueId() != nextResidueBond0AtomId.residueId(),
		    BF("The atom that is supposed to be in the next residue has the same residueId[%d]")
		    % bond0Id.residueId() );
	    RootAtomInfo const* rootAtomInfo = nextResidueBond0.get()->rootAtomInfo();
	    ASSERTF(rootAtomInfo!=NULL,BF("The RootAtomInfo must never be NULL"));
	    chem::ConstitutionAtomIndex0N bond1Id = rootAtomInfo->_Bond1Id;
	    chem::AtomId bond1AtomId(nextResidueBond0AtomId.moleculeId(),
				     nextResidueBond0AtomId.residueId(),
				     bond1Id );
	    return tree->lookup(bond1AtomId);
	} else
	{
	    // Delay is for InternalResidueAtom
//	    chem::ConstitutionAtomIndex0N id = this->_DelayAtomId;
	    chem::AtomId bond0Id(this->_Id.moleculeId(),this->_Id.residueId(),this->_DelayAtomId);
	    RefCountedAtomHandle outBond0 = tree->lookup(bond0Id);
	    return outBond0;
	}
    }


    core::Symbol_sp DelayedBondedAtom::typeSymbol() const {_OF(); return _sym_delayed;};



};
