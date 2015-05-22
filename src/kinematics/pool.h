#ifndef	kinematics_pool_H
#define kinematics_pool_H

#include "core/foundation.h"
#include "core/object.h"
#include "core/lisp.h"

namespace kinematics
{


/*! Virtual base class for any object that can be put in a Pool */

    class PoolMember
    {
    public:

	virtual uint nextUnusedMember() const = 0;
	virtual void setNextUnusedMember(uint idx) = 0;
    };


    template <class T>
	class Pool
    {
    private:
	core::Lisp_sp	_Lisp;
	string		_Comment;
	uint		_FirstUnusedMember;
	vector<T>	_Pool;
	T		_Empty;
    public:
	int size() const { return this->_Pool.size(); };
	core::Lisp_sp lisp() const {return _Lisp;};

	/*! Index into the pool */
	T& operator[](uint i) { return this->_Pool[i];};

	/*! Index into the pool */
	T const &  operator[](uint i) const { return this->_Pool[i];};

	/*! Return the index to the linked list of unused members */
	uint firstUnusedMember() const { return this->_FirstUnusedMember;};

	/*! increase the size of the pool and move new entries into an empty list */
	void resize(int newSize)
	{_OF();
	    int currentSize = this->_Pool.size();
	    ASSERTF(currentSize<newSize,BF("The requested size[%d] of the Pool must be larger than the old size[%d]") % newSize % currentSize );
	    this->_Pool.resize(newSize,this->_Empty);
	    // Add the new entries to the empty list
	    for ( int i=newSize-1; i>=currentSize; i-- )
	    {
		this->_Pool[i].setNextUnusedMember(this->_FirstUnusedMember);
		this->_FirstUnusedMember = i;
	    };
	};

	/*! Allocate a node, return its index and if there are none in the empty list then increase the size of 
	  the pool */
	uint allocate(const T& entry)
	{_OF();
	    if ( this->_FirstUnusedMember == UndefinedUnsignedInt )
	    {
		this->resize(this->_Pool.size()*2);
	    }
	    uint memberIndex = this->_FirstUnusedMember;
	    PoolMember& alloc = this->_Pool[memberIndex];
	    this->_FirstUnusedMember = alloc.nextUnusedMember();
	    this->_Pool[memberIndex] = entry;
	    LOG(BF("ALLOC pool[%s] allocating idx[%d]")
		% this->_Comment % memberIndex );
	    return memberIndex;
	}


	/*! Allocate a node, return its index and if there are none in the empty list then increase the size of 
	  the pool */
	void release(uint idx)
	{_OF();
	    LOG(BF("RELEASE pool[%s] releasing idx[%d]")
		% this->_Comment
		% idx );
	    PoolMember& releaseMember = this->_Pool[idx];
	    releaseMember.setNextUnusedMember(this->_FirstUnusedMember);
	    this->_FirstUnusedMember = idx;
	}

	string dump() const
	{
	    stringstream ss;
	    ss << "Pool[" << this->_Comment << "]";
	    ss << "  _FirstUnusedMember = " << this->_FirstUnusedMember << endl;
	    for ( int zz=0; zz<(int)this->_Pool.size(); zz++ )
	    {
		ss << "    Element["<<zz<<"] = " << this->_Pool[zz].nextUnusedMember() << endl;
	    }
	    return ss.str();
	}
	// ctor
    Pool(core::Lisp_sp lisp, const string& comment) : _Lisp(lisp),_Comment(comment), _FirstUnusedMember(UndefinedUnsignedInt) {};

    };



};

#endif // kinematics_pool_H
