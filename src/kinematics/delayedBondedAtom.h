#ifndef	kinematics_delayedBondedAtom_H
#define kinematics_delayedBondedAtom_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondedAtom.h>



namespace kinematics
{

    typedef enum { delayForFollowingResidueBond1, delayForInternalResidueAtom } DelayType;


    class DelayedBondedAtom : public BondedAtom
    {
	friend class CheckpointOutPlugAtom_O;
	friend class CheckpointAtom_O;
	friend class Checkpoint_O;
    public:
	static const NodeType nodeType = delayedBondedAtom;
    protected:
	//! Store information that describes what atom this needs to build its children
	DelayType	_DelayType;
	/*! Store the index of the atom that we are delaying for if
	  it is in the same Residue or if we are delaying for the incoming Bond1 atom
	  in the next residue then store the index of the atom in our Residue that
	  connects to the next residue */
	chem::ConstitutionAtomIndex0N	_DelayAtomId;

    private:
    public:
    DelayedBondedAtom() : BondedAtom() {};
    DelayedBondedAtom(const chem::AtomId& atomId,const string& comment) : BondedAtom(atomId,comment) {};

	virtual core::Symbol_sp typeSymbol() const;

	/*! Return the stubAtom3 */
	virtual RefCountedAtomHandle stubAtom3(AtomTree_sp tree) const;

	/*! set the Delay type and the id of the atom we are delaying for
	  If type == delayForFollowingResidueBond1 then the atomId is
	  the atom within our residue that is bonded to the incoming Bond0 atom of
	  the next residue and we are delaying for that residues incoming Bond1 atom
	  If type == delayForInternalResidueAtom then atomId is the id
	  of the atom that we are delaying for */
	void setDelay(DelayType type, chem::ConstitutionAtomIndex0N atomId);
    };



};
#endif
