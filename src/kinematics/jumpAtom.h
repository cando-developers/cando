#ifndef	kinematics_jumpAtom_H
#define kinematics_jumpAtom_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/atom.h>
#include <cando/kinematics/jump.h>
#include <cando/chem/atomId.h>

namespace kinematics
{

    class JumpAtom : public Atom
    {
    public:
	static const NodeType nodeType = jumpAtom;

    protected:
	Jump			_Jump;
	/*! JumpAtoms can have unlimited numbers of children */
	vector< RefCountedAtomHandle >	_Children;
    protected:
	/*! Bonded atoms can have different numbers of children wrt JumpAtoms */
	virtual int _maxNumberOfChildren() const { return INT_MAX;};
	/*! Return the current number of children */
	virtual int _numberOfChildren() const {return this->_Children.size();};
	/*! Get the value of the indexed child */
	virtual RefCountedAtomHandle& _child(int idx) {return this->_Children[idx];};
	virtual RefCountedAtomHandle const& _child(int idx) const {return this->_Children[idx];};
	/*! Set the value of a child */
	virtual void _setChild(int idx,const RefCountedAtomHandle& v) {this->_Children[idx] = v;};
	/*! Delete the child at the given index */
	virtual void _releaseChild(int idx); //  {this->_Children.erase(this->_Children.begin()+idx);};
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
	virtual void _insertChild(int idx, const RefCountedAtomHandle& c); //  {this->_Children.insert(this->_Children.begin()+idx,c);};
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
	virtual void _appendChild(const RefCountedAtomHandle& c); //  {this->_Children.push_back(c);};

	/*! Delete all children for dtor */
	virtual void _releaseAllChildren();

    public:
	/*! Empty ctor */
	JumpAtom() {};

    JumpAtom(const chem::AtomId& atomId, const string& comment) : Atom(atomId,comment) {};



	virtual core::Symbol_sp typeSymbol() const;


	/*! Update the internal coordinates */
	virtual void updateInternalCoords(Stub& stub,
					  bool const recursive,
	    AtomTree_sp at);


	/*! Yes, this is a JumpAtom */
	bool isJump() const { return true;};

	/*! Return the stubAtom1 */
	virtual RefCountedAtomHandle stubAtom1() const;

	/*! Return the stubAtom2 */
	virtual RefCountedAtomHandle stubAtom2() const;

	/*! Return the stubAtom3 */
	virtual RefCountedAtomHandle stubAtom3(AtomTree_sp at) const;

	bool keepDofFixed(DofType dof) const;


	/*! Update the external coordinates */
	virtual void updateXyzCoords(AtomTree_sp at);


	/*! Update the external coordinates using the input stub */
	virtual void updateXyzCoords(Stub& stub,AtomTree_sp at);


	/*! Return the DOF */
	double dof(DofType const& dof) const;
	
    };



};
#endif
