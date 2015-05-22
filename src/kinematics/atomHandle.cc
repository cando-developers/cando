#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "atomTree.h"
#include "atomHandle.h"

#define	DEBUG_REFCOUNT	0

namespace kinematics
{
    core::Lisp_sp WeakAtomHandle::lisp() const
    {
	if ( this->_Tree == NULL )
	{
	    THROW_HARD_ERROR(BF("Tried to get lisp object for undefined WeakAtomHandle"));
	}
	return this->_Tree->lisp();
    }


    AtomTree_sp WeakAtomHandle::atomTree() const
    {
	if ( this->_Tree==NULL )
	{
	    THROW_HARD_ERROR(BF("Can not return AtomTree for undefined WeakAtomHandle"));
	}
	return this->_Tree->sharedThis<AtomTree_O>();
    }


    bool WeakAtomHandle::isJump() const
    {
	AtomHolder const* holder = this->holder();
	return ( holder->_Type == jumpAtom || holder->_Type == originJumpAtom );
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
	case jumpAtom:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_JumpAtoms.size());
	    result = &(this->_Tree->_JumpAtoms[atomHolder._NodeIndex]);
	    break;
	case bondedAtom:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_BondedAtoms.size());
	    result = &(this->_Tree->_BondedAtoms[atomHolder._NodeIndex]);
	    break;
	case rootBondedAtom:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_RootBondedAtoms.size());
	    result = &(this->_Tree->_RootBondedAtoms[atomHolder._NodeIndex]);
	    break;
	case originJumpAtom:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_OriginJumpAtoms.size());
	    result = &(this->_Tree->_OriginJumpAtoms[atomHolder._NodeIndex]);
	    break;
	case delayedBondedAtom:
	    ASSERT((int)atomHolder._NodeIndex<(int)this->_Tree->_DelayedBondedAtoms.size());
	    result = &(this->_Tree->_DelayedBondedAtoms[atomHolder._NodeIndex]);
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
	case jumpAtom:
	    result = &(this->_Tree->_JumpAtoms[atomHolder._NodeIndex]);
	    break;
	case bondedAtom:
	    result = &(this->_Tree->_BondedAtoms[atomHolder._NodeIndex]);
	    break;
	case rootBondedAtom:
	    result = &(this->_Tree->_RootBondedAtoms[atomHolder._NodeIndex]);
	    break;
	case originJumpAtom:
	    result = &(this->_Tree->_OriginJumpAtoms[atomHolder._NodeIndex]);
	    break;
	case delayedBondedAtom:
	    result = &(this->_Tree->_DelayedBondedAtoms[atomHolder._NodeIndex]);
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





    RefCountedAtomHandle::RefCountedAtomHandle(AtomTree_O* tree, uint handleIndex ) : WeakAtomHandle(tree,handleIndex)
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
	if ( this->isDefined() )
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
