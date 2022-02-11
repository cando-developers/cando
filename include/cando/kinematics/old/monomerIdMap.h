/*
    File: monomerIdMap.h
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
#ifndef	_MonomerIdMap_H
#define _MonomerIdMap_H

#include <clasp/core/common.h>
#include <cando/kinematics/monomerId.h>

namespace kinematics
{

template <class T>
class MonomerIdMap
{
public:
    typedef	vector<T>		MonomerMap;
    typedef	vector<MonomerMap>	SequenceMap;
public:
	// Functions here
    /*! Resize the number of sequences */
    void resize(uint nseqs)
    {
	HARD_ASSERTF(nseqs>0,BF("Illegal number of sequences[%d]") % nseqs );
	this->_SequenceMap.resize(nseqs);
    };
    /*! Resize the number of monomers in a sequence */
    void resize(uint seq, uint nmon)
    {
	HARD_ASSERTF(seq<this->_SequenceMap.size(),BF("Illegal sequence[%d]") % seq );
	HARD_ASSERTF(nmon>=0,BF("Illegal number of monomers[%d]") % nmon );
	this->_SequenceMap[seq].resize(nmon);
    }

    /*! Resize the number of atoms in a monomer of a sequence */
    void resize(uint seq, uint mon, uint natoms)
    {
	HARD_ASSERTF(seq<this->_SequenceMap.size(),BF("Illegal sequence[%d]") % seq );
	HARD_ASSERTF(mon<this->_SequenceMap[seq].size(),BF("Illegal monomer[%d]") % mon );
	HARD_ASSERTF(natoms>=0,BF("Illegal number of atoms[%d]") % natoms );
	this->_SequenceMap[seq][mon].resize(natoms);
    }



public: // Access
    T const& operator[](const MonomerId& monomerId) const
    {
	HARD_ASSERTF((monomerId.chainId()!=MonomerId::UndefId) && (monomerId.monomerId()!=MonomerId::UndefId),
		     BF("MonomerId is undefined"));
	HARD_ASSERTF(monomerId.chainId()<(int)this->_SequenceMap.size(),
		     BF("Sequence index[%d] is out of range - number of sequences[%d]")
		     % monomerId.chainId() % this->_SequenceMap.size() );
	HARD_ASSERTF(monomerId.monomerId()<(int)this->_SequenceMap[monomerId.chainId()].size(),
		     BF("Monomer index[%d] is out of range - number of monomers[%d] in sequence[%d]")
		     % monomerId.monomerId()
		     % this->_SequenceMap[monomerId.chainId()].size()
		     % monomerId.chainId() );
	return this->_SequenceMap[monomerId.chainId()][monomerId.monomerId()];
    }
    
    T& operator[](const MonomerId& monomerId)
    {
	HARD_ASSERTF((monomerId.chainId()!=MonomerId::UndefId) && (monomerId.monomerId()!=MonomerId::UndefId),
		     BF("MonomerId is undefined"));
	HARD_ASSERTF(monomerId.chainId()<this->_SequenceMap.size(),
		     BF("Sequence index[%d] is out of range - number of sequences[%d]")
		     % monomerId.chainId() % this->_SequenceMap.size() );
	HARD_ASSERTF(monomerId.monomerId()<this->_SequenceMap[monomerId.chainId()].size(),
		     BF("Monomer index[%d] is out of range - number of monomers[%d] in sequence[%d]")
		     % monomerId.monomerId()
		     % this->_SequenceMap[monomerId.chainId()].size()
		     % monomerId.chainId() );
	return this->_SequenceMap[monomerId.chainId()][monomerId.monomerId()];
    }
    

private:
	// instance variables here
    T		_DefaultValue;
    SequenceMap	_SequenceMap;
};

}; /* kinematics */


#endif /* _MonomerIdMap_H */
