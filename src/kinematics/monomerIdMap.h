#ifndef	_MonomerIdMap_H
#define _MonomerIdMap_H

#include "core/common.h"
#include "monomerId.h"

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


