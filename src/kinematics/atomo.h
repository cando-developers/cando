#ifndef	_kinematics_AtomO_H
#define _kinematics_AtomO_H

#include <clasp/core/common.h>
#include <cando/candoBase/vector3.h>
#include <clasp/core/multipleValues.fwd.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/atomHandle.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

    FORWARD(AtomTree);

    /*! @class Wraps an AtomTree and an AtomHandle and provides
      safe, reference counted access to the underlying Atom.
    */
    class Atom_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(kinematics,KinPkg,Atom_O,"Atom");
    private:
	AtomTree_sp		_AtomTree;
	RefCountedAtomHandle	_Handle;
    public:
    Atom_O() : _Handle()
	{
	    // This should only happen when the nil object for this class is created
	    this->_AtomTree.reset();
	}
	Atom_O(const AtomTree_sp& tree, const RefCountedAtomHandle& handle);


	//! Copy ctor increments ref count
	Atom_O(const Atom_O& other);
	virtual ~Atom_O();

    public:
	static Atom_sp create(const RefCountedAtomHandle& handle)
	{
	    GC_ALLOCATE(Atom_O,atom);
	    atom->_AtomTree = handle.atomTree();
	    atom->_Handle = handle;
	    return atom;
	}
	static Atom_sp create(const AtomTree_sp& tree, uint handleIndex)
	{
	    GC_ALLOCATE(Atom_O,atom);
	    atom->_AtomTree = tree;
	    RefCountedAtomHandle handle(tree.get(),handleIndex);
	    try
	    {
		atom->_Handle = handle;
		handle.reset();
	    } catch (...)
	      {
		  handle.reset();
		  throw;
	      }
	    return atom;
	}


	RefCountedAtomHandle& handle() { return this->_Handle;};
	RefCountedAtomHandle const & handle() const { return this->_Handle;};


	Atom const& operator*() const
	{
	    return *(this->_Handle.get());
	}

	Atom* operator->()
	{
	    return this->_Handle.get();
	}


	Atom* get()
	{
	    if ( this->_Handle.notDefined()) return NULL;
	    return this->_Handle.get();
	}

	Atom const* get() const
	{
	    if ( this->_Handle.notDefined()) return NULL;
	    return this->_Handle.get();
	}


	Atom_O& operator=(const Atom_O& other);

	string asString() const;

	/*! Return the type of atom as a string */
	core::Symbol_sp typeSymbol() const;

	/*! Return the AtomId */
	chem::AtomId atomId() const;

	/*! Return the AtomId as MultipleValues */
	core::T_mv atomIdValues() const;

	/*! Return the comment for the Atom */
	string comment() const;

	/*! Append the atom as a child of this atom, this wraps appendChild */
	void appendChild(Atom_sp atom);


	/*! Insert the atom as a child of this atom, this wraps insertChild */
	void insertChild(Atom_sp atom);


	/*! Erase the atom */
	void eraseChild(Atom_sp atom);

	/*! Wrap setComment */
	Atom_sp setComment(const string& comment);

	/*! Return the index of the child of this Atom */
	int indexOfChild(Atom_sp atom);

	/*! Return a string representation */
	string __repr__() const;

	/*! Set the position of the Atom */
	void setPosition(const Vector3& pos);

	/*! Get the position of the Atom */
	Vector3 getPosition() const;

	/*! Get parent */
	Atom_sp parent() const;


	/*! Get stubAtom1 */
	Atom_sp stubAtom1() const;

	/*! Get stubAtom2 */
	Atom_sp stubAtom2() const;

	/*! Get stubAtom3 */
	Atom_sp stubAtom3(AtomTree_sp tree) const;


#if 0
	/*! Build a Molecule using the graph of monomers described by MonomerNode as a plan */
	void buildChainUsingchainNode(int moleculeId, ChainNode_sp sequenceNode );
#endif
    };







}; /* kinematics */

TRANSLATE(kinematics::Atom_O);

#endif /* _kinematics_Atom_H */


