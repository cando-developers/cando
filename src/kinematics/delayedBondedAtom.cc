#define DEBUG_LEVEL_FULL



#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/kinematics/symbolTable.h>
#include <cando/kinematics/rootAtomInfo.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/delayedBondedAtom.h>

namespace kinematics
{

    RefCountedAtomHandle DelayedBondedAtom::stubAtom3(AtomTree_sp tree) const
    {_G();
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
