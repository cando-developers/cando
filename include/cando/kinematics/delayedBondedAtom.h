/*
    File: delayedBondedAtom.h
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
#ifndef	kinematics_delayedBondedJoint_H
#define kinematics_delayedBondedJoint_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondedAtom.h>



namespace kinematics
{

    typedef enum { delayForFollowingResidueBond1, delayForInternalResidueAtom } DelayType;

    FORWARD(DelayedBondedJoint);
    class DelayedBondedJoint_O : public BondedJoint_O
    {
	LISP_CLASS(kinematics,KinPkg,DelayedBondedJoint_O,"DelayedBondedAtom",BondedJoint_O);
    public:
	static const NodeType nodeType = delayedBondedAtom;
    public:
	//! Store information that describes what atom this needs to build its children
	DelayType	_DelayType;
	/*! Store the index of the atom that we are delaying for if
	  it is in the same Residue or if we are delaying for the incoming Bond1 atom
	  in the next residue then store the index of the atom in our Residue that
	  connects to the next residue */
	chem::ConstitutionAtomIndex0N	_DelayAtomId;
    public:
    DelayedBondedJoint_O() : BondedJoint_O() {};
    DelayedBondedJoint_O(const chem::AtomId& atomId,const string& comment) : BondedJoint_O(atomId,comment) {};

	virtual core::Symbol_sp typeSymbol() const;

	/*! Return the stubAtom3 */
	virtual Joint_sp stubAtom3(AtomTree_sp tree) const;

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
