/*
    File: pool.h
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
#ifndef	kinematics_pool_H
#define kinematics_pool_H

#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>

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
	    ss << "  _FirstUnusedMember = " << this->_FirstUnusedMember << std::endl;
	    for ( int zz=0; zz<(int)this->_Pool.size(); zz++ )
	    {
              ss << "    Element["<<zz<<"] = " << this->_Pool[zz].nextUnusedMember() << std::endl;
	    }
	    return ss.str();
	}
	// ctor
    Pool(core::Lisp_sp lisp, const string& comment) : _Lisp(lisp),_Comment(comment), _FirstUnusedMember(UndefinedUnsignedInt) {};

    };



};

#endif // kinematics_pool_H
