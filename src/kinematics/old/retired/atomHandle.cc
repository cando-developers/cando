/*
    File: atomHandle.cc
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
#define	DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <cando/kinematics/pointTree.h>
#include <cando/kinematics/pointHandle.h>

#define	DEBUG_REFCOUNT	0

namespace kinematics
{

    PointTree_sp WeakAtomHandle::pointTree() const
    {
	if ( this->_Tree==NULL )
	{
	    THROW_HARD_ERROR(BF("Can not return PointTree for undefined WeakAtomHandle"));
	}
	return this->_Tree->sharedThis<PointTree_O>();
    }


    bool WeakAtomHandle::isJump() const
    {
	AtomHolder const* holder = this->holder();
	return ( holder->_Type == jumpPoint || holder->_Type == originJumpPoint );
    }

    RefCountedAtomHandle WeakAtomHandle::lock() const
    {
	RefCountedAtomHandle h(*this);
	return h;
    }

    int WeakAtomHandle::refCount() const
    {
	return this->holder()->refCount();
    };


    Atom* WeakAtomHandle::get()
    {_OF();
	if ( this->notDefined() )
	{
	    SIMPLE_ERROR(BF("This handle is Undefined and you tried to get its pointed to object"));
	}
	HARD_ASSERTF((int)this->_HolderIndex<this->_Tree->numberOfAtoms(),
		     BF("Illegal handle index[%d] must be less than %d")
		     % this->_HolderIndex % this->_Tree->numberOfAtoms() );
	ASSERT((int)this->_HolderIndex<(int)this->_Tree->_AtomHolders.size());
	AtomHolder& atomHolder = this->_Tree->_AtomHolders[this->_HolderIndex];
	Atom* result = NULL;
	switch (atomHolder._Type)
	{
	case jumpPoint:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_JumpPoints.size());
	    result = &(this->_Tree->_JumpPoints[atomHolder._NodeIndex]);
	    break;
	case bondedPoint:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_BondedPoints.size());
	    result = &(this->_Tree->_BondedPoints[atomHolder._NodeIndex]);
	    break;
	case rootBondedPoint:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_RootBondedPoints.size());
	    result = &(this->_Tree->_RootBondedPoints[atomHolder._NodeIndex]);
	    break;
	case originJumpPoint:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_OriginJumpPoints.size());
	    result = &(this->_Tree->_OriginJumpPoints[atomHolder._NodeIndex]);
	    break;
	case delayedBondedPoint:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_DelayedBondedPoints.size());
	    result = &(this->_Tree->_DelayedBondedPoints[atomHolder._NodeIndex]);
	    break;
	default:
	    SIMPLE_ERROR(BF("Could not dereference AtomHolder - unknown type"));
	}
	if ( result == NULL )
	{
	    SIMPLE_ERROR(BF("About to dereference AtomHandle but it will return NULL"));
	}
	return result;
    }

    Atom const* WeakAtomHandle::get() const
    {_OF();
	if ( this->notDefined() )
	{
	    SIMPLE_ERROR(BF("This handle is Undefined and you tried to get its pointed to object"));
	}
	HARD_ASSERTF((int)this->_HolderIndex<this->_Tree->numberOfAtoms(),
		     BF("Illegal handle index[%d] must be less than %d")
		     % this->_HolderIndex % this->_Tree->numberOfAtoms() );
	AtomHolder& atomHolder = this->_Tree->_AtomHolders[this->_HolderIndex];
	Atom* result = NULL;
	switch (atomHolder._Type)
	{
	case jumpPoint:
	    result = &(this->_Tree->_JumpPoints[atomHolder._NodeIndex]);
	    break;
	case bondedPoint:
	    result = &(this->_Tree->_BondedPoints[atomHolder._NodeIndex]);
	    break;
	case rootBondedPoint:
	    result = &(this->_Tree->_RootBondedPoints[atomHolder._NodeIndex]);
	    break;
	case originJumpPoint:
	    result = &(this->_Tree->_OriginJumpPoints[atomHolder._NodeIndex]);
	    break;
	case delayedBondedPoint:
	    result = &(this->_Tree->_DelayedBondedPoints[atomHolder._NodeIndex]);
	    break;
	default:
	    SIMPLE_ERROR(BF("Could not dereference AtomHolder - unknown type"));
	}
	if ( result == NULL )
	{
	    SIMPLE_ERROR(BF("About to dereference AtomHandle but it will return NULL"));
	}
	return result;
    }



    AtomHolder* WeakAtomHandle::holder()
    {_OF();
	ASSERTF(!this->notDefined(), BF("Undefined WeakAtomHandle"));
	return &(this->_Tree->_AtomHolders[this->_HolderIndex]);
    }

    const AtomHolder* WeakAtomHandle::holder() const
    {_OF();
	ASSERTF(!this->notDefined(), BF("Undefined WeakAtomHandle"));
	return &(this->_Tree->_AtomHolders[this->_HolderIndex]);
    }





    RefCountedAtomHandle::RefCountedAtomHandle(PointTree_O* tree, uint handleIndex ) : WeakAtomHandle(tree,handleIndex)
    {
	if ( !this->notDefined() )
	{
	    this->incrementRefCount();
	}
    };







    RefCountedAtomHandle::RefCountedAtomHandle(const WeakAtomHandle& handle) : WeakAtomHandle(handle)
	{
	    if ( !this->notDefined() )
	    {
		this->incrementRefCount();
	    }
	};



    RefCountedAtomHandle::RefCountedAtomHandle(const RefCountedAtomHandle& handle)
    {
	this->_Tree = handle._Tree;
	this->_HolderIndex = handle._HolderIndex;
	if ( this->boundp() )
	{
	    this->incrementRefCount();
	}
    }






    void RefCountedAtomHandle::reset()
    {
	if ( !this->notDefined() )
	{
	    this->decrementRefCountAndReleaseIfZero();
	}
	this->resetWeak();
    }


#if 0
	/*! Set the shared handle */
    void RefCountedAtomHandle::assign(const WeakAtomHandle& handle)
	{
	    if ( this->_HolderIndex == handle.holderIndex() ) return;
	    if ( !this->notDefined() )
	    {
		this->decrementRefCountAndReleaseIfZero();
	    }
	    this->assignWeak(handle);
	    if ( (int)handle.holderIndex() < this->_Tree->numberOfAtoms() )
	    {
		this->incrementRefCount();
	    }
	}
#endif


    void RefCountedAtomHandle::incrementRefCount()
    {_OF();
#if DEBUG_REFCOUNT
	printf("\n%s:%d - About to increment refCount->[%d] for handle[%d] nodeType[%s]",
	       __FILE__, __LINE__, this->holder()->refCount(), this->_HolderIndex,
	       this->holder()->typeAsString().c_str());
	LOG(BF("REFCOUNT-- About to increment refCount->[%d] for handle[%d] nodeType[%s]")
	    % this->holder()->refCount() % this->_HolderIndex % this->holder()->typeAsString() );
#endif
	LOG(BF("Incrementing RefCount for handle[%d]") % this->_HolderIndex );
	ASSERTF(!this->notDefined(),BF("Tried to increment RefCount of undefined handle"));
	ASSERTF((int)this->_HolderIndex<this->_Tree->numberOfAtoms(),BF("Illegal Handle for incrementRefCount"));
	ASSERTF(this->holder()->_Type!=unused,
		BF("incrementRefCount>>You tried to access an unused node[%d]") % this->_HolderIndex );
	this->_Tree->_AtomHolders[this->_HolderIndex].refAdd();
	LOG(BF("REFCOUNT++ incremented refCount->[%d] for handle[%d]")
	    % this->holder()->refCount() % this->_HolderIndex );
#if DEBUG_REFCOUNT
	printf("\n%s:%d - incremented refCount[%d] for handle[%d]",
	       __FILE__, __LINE__, this->holder()->refCount(), this->_HolderIndex);
#endif
    }

    bool RefCountedAtomHandle::decrementRefCountAndReleaseIfZero()
    {_OF();
#if DEBUG_REFCOUNT
	printf("\n%s:%d - About to decrement refCount->[%d] for handle[%d] nodeType[%s]",
	       __FILE__, __LINE__, this->holder()->refCount(), this->_HolderIndex,
	       this->holder()->typeAsString().c_str());
	LOG(BF("REFCOUNT-- About to decrement refCount->[%d] for handle[%d] nodeType[%s]")
	    % this->holder()->refCount() % this->_HolderIndex % this->holder()->typeAsString() );
#endif
	ASSERTF(!this->notDefined(),BF("Tried to increment RefCount of undefined handle"));
	ASSERTF((int)this->_HolderIndex<this->_Tree->numberOfAtoms(),BF("Illegal Handle for incrementRefCount"));
	ASSERTF(this->holder()->_Type!=unused,
		BF("decrementRefCountAndReleaseIfZero>>You tried to access an unused node[%d]") % this->_HolderIndex );
	bool isZero = this->_Tree->_AtomHolders[this->_HolderIndex].refDec();
	LOG(BF("REFCOUNT-- decremented refCount->[%d] for handle[%d]")
	    % this->holder()->refCount() % this->_HolderIndex );
#if DEBUG_REFCOUNT
	printf("\n%s:%d - decremented refCount->[%d] for handle[%d]",
	       __FILE__, __LINE__, this->holder()->refCount(), this->_HolderIndex );
	if ( isZero )
	{
	    printf("\n%s:%d - releasing handle[%d]", __FILE__, __LINE__, this->_HolderIndex);
	}
#endif
	if ( isZero ) 
	{
	    this->_Tree->releaseHandle(*this);
	    this->resetWeak();
	}
	return isZero;
    }




};
