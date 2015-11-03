#ifndef	kinematics_atom_H
#define kinematics_atom_H

#include <clasp/core/foundation.h>
#include <cando/geom/vector3.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/dofType.h>
#include <cando/kinematics/stub.fwd.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/atomHandle.h>
#include <cando/kinematics/pool.h>


// Currently a maximum of 5 children hanging off of any atom
// means we can handle octahedral coordination but no more

namespace kinematics
{

 


#define	DEBUG_KIN_ATOM	1


    class RootAtomInfo;


/*! Atom stores a link to a parent atom and a list of children atoms that
  are reference counted.  The children atoms can be of different types and the order of the 
  children is enforced as JumpAtoms | BondedAtom atoms.
 */

    class Atom : public PoolMember
    {
	friend class Atom_O;
	friend class AtomTree_O;
	friend class BondedAtom;
    protected:
	WeakAtomHandle	_Me;
	//! Point to the parent atom (also used to contruct linked list of unused PoolMembers)
	WeakAtomHandle	_Parent;
	chem::AtomId	_Id;
	Vector3		_Position;
#if DEBUG_KIN_ATOM
	string		_Comment;
	Vector3		_TestPosition;
#endif
    private:
	/*! Track my position in my owner's list of Atoms with modified DOFs
	  -1 when my dofs have not changed since the last update_coords
	*/
	int		_DofRefoldIndex;
    protected:
	/*! Bonded atoms can have different numbers of children wrt JumpAtoms */
 	virtual int _maxNumberOfChildren() const = 0;
	/*! Return the current number of children */
	virtual int _numberOfChildren() const = 0;
	/*! Return a reference to the indexed child */
	virtual RefCountedAtomHandle& _child(int idx) = 0;
	/*! Return a reference to the indexed child */
	virtual RefCountedAtomHandle const& _child(int idx) const = 0;
	/*! Set a value of a child */
	virtual void _setChild(int idx,const RefCountedAtomHandle& atom) = 0;
	/*! Delete the child at the given index */
	virtual void _releaseChild(int idx) = 0;
	/*! Insert the child at the given index - this does the
	  work of opening up a space and putting the new value in */
	virtual void _insertChild(int idx, const RefCountedAtomHandle& atom) = 0;
	/*! Insert the child at the given index - this does the work
	  of opening up a space and putting the new value in */
	virtual void _appendChild(const RefCountedAtomHandle& atom) = 0;

	/*! Destructors need to delete all Children */
  	virtual void _releaseAllChildren() = 0;
    public:

	Atom() : _Parent(), _Id() {};
	Atom(const chem::AtomId& atomId, const string& comment) :
	    _Parent(),
	    _Id(atomId)
#if DEBUG_KIN_ATOM
	    , _Comment(comment)
#endif
	{};

	RefCountedAtomHandle atomHandle() const { return this->_Me;};
	/*! Set the AtomHandle */
	void __setAtomHandle(const WeakAtomHandle& handle) { this->_Me = handle;};


	/*! Return the lisp object */
	core::Lisp_sp lisp() const;

	chem::AtomId id() const { return this->_Id;};
	virtual core::Symbol_sp typeSymbol() const;
	virtual string asString() const;

	virtual void setNextUnusedMember(uint idx) { this->_Parent.setNextUnusedMember(idx);};
	virtual uint nextUnusedMember() const { return this->_Parent.getNextUnusedMember();};

	/*! Set the parent, the parent is a weak reference so we don't get reference cycles */

	void setParent(const RefCountedAtomHandle& parent) { this->_Parent = parent;};

	/*! Return a RefCountedAtomHandle for the parent */
	RefCountedAtomHandle parent() const { return this->_Parent; };

	/*! Insert the child before the (before) index. */
	void insertChild( int before, RefCountedAtomHandle& child );

	/*! Insert the child.
	  If the child is a JumpAtom then put it as the first child.
	  If it's a Bonded atom then put it before all the existing BondedAtoms.
	*/
	void insertChild(RefCountedAtomHandle& child);

	
	/*! Return the index of the first non-jump child */
	int firstNonJumpChildIndex() const;


	/*! Append the child. 
	  If the child is a JumpAtom then put it before all the other non-JumpAtoms.
	  If it is a Bonded atom then put it at the end.
	*/
	void appendChild(RefCountedAtomHandle& handle);

	/*! Erase the child from the Children using the index into the _Children array */
	void erase(int idx );


	RefCountedAtomHandle& child(int idx) { return this->_child(idx);};
	/*! Return a reference to the indexed child */
	RefCountedAtomHandle const& child(int idx) const { return this->_child(idx);};


	/*! Root nodes will return RootAtomInfo structures */
	virtual RootAtomInfo const* rootAtomInfo() const;

	/*! For debugging dump the AtomTree */
	void recursiveDumpChildrenIntoStringStream(const string& prefix,
						   stringstream& out);

	/*! Update the internal coordinates */
	virtual void updateInternalCoords(Stub& stub,
					  bool const recursive,
					  AtomTree_sp at 	) = 0;

	/*! Return true if this atom is a JumpAtom (or subclass) */
	virtual bool isJump() const { return false;};


	/*! Return the i(th) non-jump atom */
	RefCountedAtomHandle getNonJumpAtom(int idx) const;

	/*! Return true if the stub is defined */
	virtual bool stubDefined() const;


	/*! Return the stub for this atom */
	Stub getStub(AtomTree_sp at) const;

	/*! Return the position */
	Vector3 position() const { return this->_Position;};

	/*! Set the position */
	void position(const Vector3& pos) { this->_Position = pos;};

	/*! Return the input stub atom */
	inline RefCountedAtomHandle inputStubAtom0() const
	{_OF();
	    ASSERTF(this->parent().isDefined(),BF("Parent isn't defined"));
	    return this->parent();
	}

	/*! Return the input stub atom */
	inline RefCountedAtomHandle inputStubAtom1() const
	{_OF();
	    ASSERTF(this->parent().isDefined(),BF("Parent isn't defined"));
	    return this->parent().get()->stubAtom1();
	}

	/*! Return the input stub atom */
	inline RefCountedAtomHandle inputStubAtom2() const
	{_OF();
	    ASSERTF(this->parent().isDefined(),BF("Parent isn't defined"));
	    return this->parent().get()->stubAtom2();
	}


	/*! Return the previous child to this one */
	RefCountedAtomHandle previousChild(const RefCountedAtomHandle& child) const;

	/*! Return the previous sibling of this atom */
	RefCountedAtomHandle previousSibling() const
	{
	    if ( this->parent().isDefined() )
	    {
		return this->parent().get()->previousChild(this->atomHandle());
	    } else
	    {
		return RefCountedAtomHandle::undefined();
	    }
	}


	/*! Return the input stub atom3
	  It is either its parents stubAtom3 or its previous sibling */
	RefCountedAtomHandle inputStubAtom3(AtomTree_sp at) const;



	/*! Return the stubAtom1 */
	virtual RefCountedAtomHandle stubAtom1() const = 0;

	/*! Return the id of stubAtom */
	chem::AtomId stubAtom1Id() const { return this->stubAtom1().get()->id();};

	/*! Return the stubAtom2 */
	virtual RefCountedAtomHandle stubAtom2() const = 0;

	/*! Return the id of stubAtom */
	chem::AtomId stubAtom2Id() const { return this->stubAtom2().get()->id();};

	/*! Return the stubAtom3 */
	virtual RefCountedAtomHandle stubAtom3(AtomTree_sp tree) const = 0;

	/*! Return the id of stubAtom */
	chem::AtomId stubAtom3Id(AtomTree_sp at) const { return this->stubAtom3(at).get()->id();};

	/*! Return true if we want to keep the Dof fixed */
	bool keepDofFixed(DofType dof) const { return false;};

	/*! Walk the tree */
	void walkChildren(const AtomTreeWalkFunctor& functor);

	/*! Walk the tree but stay in one residue */
	void walkResidueTree(int residueId, const AtomTreeWalkFunctor& functor);


	/*! Update the external coordinates after calculating the input stub */
	virtual void updateXyzCoords(AtomTree_sp at);

	/*! Update the external coordinates using the input stub */
	virtual void updateXyzCoords(Stub& stub,AtomTree_sp at) {THROW_HARD_ERROR(BF("Subclass must implement"));};

	/*! Ensure proper function of the output-sensitive refold subroutine
	  derived classes must invoke this function during their updateXyzCoords subroutines
	*/
	void noteXyzUpToDate()
	{
	    this->_DofRefoldIndex = 0;
	}


	/*! Return the input stub */
	Stub getInputStub(AtomTree_sp at) const;


	/*! Return the value of the DOF */
	virtual double dof(DofType const& dof) const { THROW_HARD_ERROR(BF("SubClass must implement"));};

    };





};


#endif
