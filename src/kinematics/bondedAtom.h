#ifndef	kinematics_bondedAtom_H
#define kinematics_bondedAtom_H

#include "core/foundation.h"
#include "kinFoundation.h"
#include "chem/atomId.h"
#include "atom.h"



namespace kinematics
{

    class BondedAtom : public Atom
    {
    public:
	static const NodeType nodeType = bondedAtom;
	static const int MaxChildren = 5;
    protected:
	int		_NumberOfChildren;
	RefCountedAtomHandle	_Children[MaxChildren];
	Real		_Phi;
	Real		_Theta;
	Real		_Distance;
	bool		_DofChangePropagatesToYoungerSiblings;
    private:
	/*! Bonded atoms can have different numbers of children wrt JumpAtoms */
	virtual int _maxNumberOfChildren() const { return MaxChildren;};
	/*! Return the current number of children */
	virtual int _numberOfChildren() const {return this->_NumberOfChildren;};
	/*! Return a reference to the indexed child */
	virtual RefCountedAtomHandle& _child(int idx) {return this->_Children[idx];};
	virtual RefCountedAtomHandle const& _child(int idx) const {return this->_Children[idx];};
	/*! Set a child */
	virtual void _setChild(int idx, const RefCountedAtomHandle& h) { this->_Children[idx] = h; };
	/*! Delete the child at the given index */
	virtual void _releaseChild(int idx);
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
	virtual void _insertChild(int idx, const RefCountedAtomHandle& c );
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
	virtual void _appendChild(const RefCountedAtomHandle& c) {this->_Children[this->_NumberOfChildren++] = c; };
	/*! Delete all of the children for the destructor */
	virtual void _releaseAllChildren();

    public:
    BondedAtom() : Atom(), _NumberOfChildren(0) {};
    BondedAtom(const chem::AtomId& atomId, const string& comment) : Atom(atomId,comment), _NumberOfChildren(0) {};

	virtual core::Symbol_sp typeSymbol() const;

	/*! Return the stubAtom1 */
	virtual RefCountedAtomHandle stubAtom1() const { return this->atomHandle();}

	/*! Return the stubAtom2 */
	virtual RefCountedAtomHandle stubAtom2() const { return this->parent();};

	/*! Return the stubAtom3 */
	virtual RefCountedAtomHandle stubAtom3(AtomTree_sp tree) const;

	/*! Update the internal coordinates */
	virtual void updateInternalCoords(Stub& stub,
					  bool const recursive,
					  AtomTree_sp at);


	bool keepDofFixed(DofType dof,AtomTree_sp at) const;


	string asString() const;


	/*! Update the external coordinates */
	virtual void updateXyzCoords(AtomTree_sp at);


	/*! Update the external coordinates using the input stub */
	virtual void updateXyzCoords(Stub& stub,AtomTree_sp at);



	/*! Get the value of the DOF */
	double dof(DofType const& dof) const;


    };



};
#endif
