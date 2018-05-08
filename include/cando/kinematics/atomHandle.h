/*
    File: atomHandle.h
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
#ifndef kinematics_atomHandle_H
#define kinematics_atomHandle_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/kinematics/kinematicsPackage.h>
//#include "pool.h"
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/bondId.fwd.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/bondId.h>


namespace kinematics
{

    class Atom;


    class RefCountedAtomHandle;

    class WeakAtomHandle
    {
    public:
	static const uint Undef = UndefinedUnsignedInt;
    protected:
	/*! Raw pointer to AtomTree - its dangerous but AtomHandles within an AtomTree should
	  never outlive their AtomTree */
	AtomTree_O*	_Tree;
	uint		_HolderIndex;
    public:
    WeakAtomHandle() : _Tree(NULL), _HolderIndex(Undef) {};
    WeakAtomHandle(AtomTree_O* tree, uint handleIndex) : _Tree(tree),_HolderIndex(handleIndex) {};

	/*! Return the AtomTree_sp for this WeakAtomHandle */
	AtomTree_sp atomTree() const;

	/*! Return true if the Handle is undefined */
	bool notDefined() const { return this->_HolderIndex == Undef;};


	/*! Return true if the Handle is defined (not undefined)*/
	bool isDefined() const { return this->_HolderIndex != Undef;};

	/*! Assignment operator */
	WeakAtomHandle& operator=(const WeakAtomHandle& other)
	{
	    this->_Tree = other._Tree;
	    this->_HolderIndex = other._HolderIndex;
	    return *this;
	}


	/*! Comparison operator */
	bool operator==(const WeakAtomHandle& h) const
	{
	    if ( this->_Tree != h._Tree ) return false;
	    return (this->_HolderIndex == h._HolderIndex);
	}

        								
	/*! Return a RefCountedAtomHandle for this WeakHandle */
	RefCountedAtomHandle lock() const;

	/*!When the HolderIndex is needed as an index for a linked list of
	  unused members, use this to set it */
	void setNextUnusedMember(uint idx) { this->_HolderIndex = idx;};
	/*! When the HolderIndex is needed as an index for a linked list of
	  unused members, use this to get it */
	uint getNextUnusedMember() const { return this->_HolderIndex;};

	/*! Return true if this handle points to a JumpAtom */
	bool isJump() const;

	/*! Reset this Handle to Undefined */
	void resetWeak()
	{
	    this->_Tree = NULL;
	    this->_HolderIndex = Undef;
	}

	/*! Return the AtomHolder for this Handle */
	AtomHolder* holder();


	/*! Return the AtomHolder for this Handle */
	AtomHolder const* holder() const;

	/*! Return the reference count */
	int refCount() const;

	/*! Return the AtomHolder index for this handle */
	uint holderIndex() const { return this->_HolderIndex;};


	/*! Return a pointer to the atom referenced by this handle */
	Atom* get();

	/*! Return a pointer to the atom referenced by this handle */
	Atom const* get() const;

	template <class U>
	    U* getAs() { return dynamic_cast<U*>(this->get());};
    };




    class RefCountedAtomHandle : public WeakAtomHandle
    {
    public:

    RefCountedAtomHandle() : WeakAtomHandle() {};

	RefCountedAtomHandle(AtomTree_O* tree, uint handleIndex );

	RefCountedAtomHandle(const WeakAtomHandle& handle);

	RefCountedAtomHandle(const RefCountedAtomHandle& handle);

	/*! Return a RefCountedAtomHandle for an undefined atom */
	static RefCountedAtomHandle undefined()
	{
	    RefCountedAtomHandle h(NULL,Undef);
	    return h;
	}


	/*! Comparison operator */
	bool operator==(const RefCountedAtomHandle& h) const
	{
	    if ( this->_Tree != h._Tree ) return false;
	    return (this->_HolderIndex == h._HolderIndex);
	}


	RefCountedAtomHandle& operator=(const RefCountedAtomHandle& other)
	    {
		// Ignore self assignments
		if ( this->_Tree == other._Tree && this->_HolderIndex == other._HolderIndex ) return *this;
		// If we already point to an Atom then decrement its ref count
		if ( this->isDefined() ) this->decrementRefCountAndReleaseIfZero();
		// Copy the contents
		this->_Tree = other._Tree;
		this->_HolderIndex = other._HolderIndex;
		// If we now point to an atom then increment its ref count
		if ( this->isDefined() ) this->incrementRefCount();
		return *this;
	    }


	~RefCountedAtomHandle()
	{
	    if (this->isDefined()) this->decrementRefCountAndReleaseIfZero();
	}


	/*! Reset the shared handle */
	void reset();


	/*! Increment the reference count */
	void incrementRefCount();

	/*! Decrement the reference count, release the handle if ref count reaches 0 */
	bool decrementRefCountAndReleaseIfZero();


    };


};


#endif
