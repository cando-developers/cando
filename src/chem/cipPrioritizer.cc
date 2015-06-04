#define	DEBUG_LEVEL_NONE

#include <cando/chem/cipPrioritizer.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/molecule.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/loop.h>
#include <clasp/core/sort.h>
#include <cando/chem/bond.h>
#include <cando/chem/virtualAtom.h>
#include <clasp/core/wrappers.h>



namespace chem {





uint atomicNumber(Atom_sp a)
{
    return a->getAtomicNumber();
}


class	OrderByBase
{
protected:
	CipPrioritizer_sp	prior;
public:
    void setCipPrioritizer(CipPrioritizer_sp p) { this->prior = p;};
    virtual bool operator() (Atom_sp x, Atom_sp y) = 0;
};


class	OrderByP : public OrderByBase
{
public:
    bool operator()(Atom_sp x, Atom_sp y )
    {
	if ( prior->getP(x) <= prior->getP(y) ) return true;
	return false;
    }
};


class	OrderByS : public OrderByBase
{
public:
    int cmpByS(Atom_sp x, Atom_sp y )
{_G();
        ASSERTNOTNULL(x);
        ASSERTNOTNULL(y);
	vector<int>&	xv = prior->getS(x);
	vector<int>&	yv = prior->getS(y);
	vector<int>::iterator	xi = xv.begin();
	vector<int>::iterator	yi = yv.begin();
	bool done=false;
	int res;
	done =  (( xi == xv.end() ) || ( yi == yv.end() ));
	while ( !done )
	{
	    if ( *xi < *yi )
	    {
		res = -1;
		goto done;
	    }
	    if ( *xi > *yi ) 
	    {
		res = 1;
		goto done;
	    }
	    xi++;
	    yi++;
	    done |=  (( xi == xv.end() ) || ( yi == yv.end() ));
	}
	if ( xi==xv.end() && yi==yv.end() ) 
	{
	    res = 0;
	    goto done;
	}
	if ( xi==xv.end() )
	{
	    	// x is shorter than y so x<=y
	    res = -1;
	    goto done;
	}
	res = 1;
    done: 
        return res;
    }


    bool operator()(Atom_sp x, Atom_sp y )
    {_G();
        int order = this->cmpByS(x,y);
	if ( order < 0 ) return true;
	if ( order == 0 ) return true;
	return false;
    }
};



void	CipPrioritizer_O::initialize()
{
    this->Base::initialize();
}

int	CipPrioritizer_O::getP(Atom_sp a)
{_OF();
    ASSERT_lessThan(a->getRelativePriority(),this->_p.size());
    return this->_p[a->getRelativePriority()];
}

vector<int>& CipPrioritizer_O::getS(Atom_sp a)
{_OF();
    ASSERTNOTNULL(a);
    if ( ! ( a->getRelativePriority() < this->_s.size() ) )
    {
	LOG(BF("Bad priority for atom") );
	LOG(BF("Bad priority for atom: %s") % a->description().c_str()  );
	LOG(BF("   priority value = %d") % a->getRelativePriority()  );
    }
    ASSERT_lessThan(a->getRelativePriority(),this->_s.size());
    return this->_s[a->getRelativePriority()];
}



/*!
 * partition the atoms between cBegin and cEnd-1 into classes based on orderFn
 */
    vector<int>	partition(gctools::Vec0<Atom_sp> atoms,int indexOffset, OrderByS& orderer)
{
    vector<int> part;
    part.push_back(indexOffset);
    for ( uint i=1; i<atoms.size(); i++ )
    {
	if ( orderer.cmpByS(atoms[i-1],atoms[i])!=0 )
	{
	    part.push_back(i+indexOffset);
	}
    }
    return part;
}
    

void CipPrioritizer_O::assignPriorities(Matter_sp matter)
{_G();
    CipPrioritizer_sp prior;
    prior = CipPrioritizer_O::create();
    prior->assignCahnIngoldPrelogPriorityToAtomsRelativePriority(matter);
}


//
// ! Assign the CIP priority to the RelativePriority of each Atom in the molecule
void CipPrioritizer_O::assignCahnIngoldPrelogPriorityToAtomsRelativePriority(Matter_sp molOrAgg)
{_G();
    gctools::Vec0<Atom_sp>	mAtoms;
    vector<int> newC;
    vector<int> C;
    // 
    //  Determination of topological RS chirality
    //  see: http://www.chemcomp.com/journal/chiral.htm
    // The objective of the algorithm is to assign a priority, p(i) to every 
    // atom i in the molecule so that if p(i)<p(j) then atom i has strictly 
    // lower CIP priority than atom j with equality only occurring when the 
    // atoms are indistinguishable (topologically). In essence, the algorithm 
    // maintains an ordered list of equivalence classes of atoms. Each atom 
    // in an equivalence class is assigned the priority of the -1 of the 
    // class in the sorted list. The algorithm repeatedly splits classes 
    // (maintaining the ordering) until no changes occur. The final priority 
    // assignments are then output as the CIP priorities.
    // 
    // The algorithm proceeds as follows:
    // 
    // 1. 	For each atom i, set p(i) equal to the atomic number of atom i. 
    // Set the initial partition to be an ordered list of classes (C1,...,Ck) 
    // such that for each atom i in class Cr and each atom j in class Cs we 
    // have (i) p(i)<p(j) iff r<s; and (ii) p(i)=p(j) iff r=s.
    //
    // Use the RelativePriority value of each atom to store its position in the 
    // p and s arrays
    //
    // 
    { _BLOCK_TRACE("Stage 1 - setting initial priorities");
	this->_p.clear();
	Loop l;
	l.loopTopGoal(molOrAgg,ATOMS);
	while ( l.advanceLoopAndProcess() )
	{
	    Atom_sp a = l.getAtom();
	    if ( a.isA<VirtualAtom_O>() ) continue;
	    a->setRelativePriority(mAtoms.size());
	    mAtoms.push_back(a);
	    this->_p.push_back(atomicNumber(a));
	}

	OrderByP byP;
	byP.setCipPrioritizer(this->sharedThis<CipPrioritizer_O>());
	if ( mAtoms.begin() != mAtoms.end() )
	{
	    core::sort::quickSort(mAtoms.begin(),mAtoms.end(),byP);
	}

	C.clear();
	int pcur = -1;
	for ( uint ai=0; ai<mAtoms.size(); ai++ )
	{
	    Atom_sp a = mAtoms[ai];
	    uint idx = a->getRelativePriority();
	    if (pcur != this->_p[idx])
	    {
		C.push_back(ai);
		pcur = this->_p[idx];
	    }
	}
#ifdef	DEBUG_ON
	int cidx = 0;
	LOG(BF("Assigning priorities STAGE1") );
	for ( gctools::Vec0<Atom_sp>::iterator mit=mAtoms.begin(); 
		    mit!=mAtoms.end(); mit++ )
	{
	    stringstream ss;
	    ss << "Atom: " << (*mit)->getName();
	    ss << " priority: " << this->getP((*mit));
	    ASSERT_lessThan(cidx,(int)(C.size()+1));
	    if ( cidx < (int)(C.size()) && mit-mAtoms.begin() == C[cidx] )
	    {
		cidx++;
	    }
	    ss << " C-class: " << cidx-1;
	    LOG(BF("%s") % ss.str().c_str()  );
	}
#endif
    }
    // C.append(len(mAtoms))

    // print "Classes go from C[i] to C[i+1] unless last i then end, C=",C


    // 2.	For each atom i set s(i)=abc...z to be an ordered list of 
    //  neighboring p(j) numbers in decreasing order accounting for bond 
    //  multiplicities with repeated values (i.e., if atom i is double 
    //  bonded to an atom with priority p, then put p in the list twice).

    {_BLOCK_TRACE("Stage 2 to 6");
		    // previously this was at the top of step 3
	bool didPartition = true;
	while ( didPartition )
	{
	    this->_s.clear();
	    this->_s.resize(mAtoms.size());
	    { _BLOCK_TRACE("Stage 2");
		for ( gctools::Vec0<Atom_sp>::iterator mi =mAtoms.begin();
				 mi!=mAtoms.end(); mi++ )
		{
		    Atom_sp myatom = *mi;
		    LOG(BF("About to fill mys") );
		    vector<int>	mys;
		    for ( gctools::Vec0<Bond_sp>::iterator bi=myatom->bonds_begin();
				 bi!=myatom->bonds_end(); bi++ )
		    {
			int bondOrder = 1;
			if ( (*bi)->getOrder() == singleBond ) bondOrder = 1;
			if ( (*bi)->getOrder() == doubleBond ) bondOrder = 2;
			if ( (*bi)->getOrder() == tripleBond ) bondOrder = 3;
			for ( int oi=0; oi<bondOrder; oi++ )
			{
			    mys.push_back(this->_p[(*bi)->getOtherAtom(myatom)->getRelativePriority()]);
			}
		    }
		    LOG(BF("About to sort %d mys objects") % mys.size()  );
		    core::sort::quickSort(mys.begin(),mys.end());
		    LOG(BF("Done sort") );
		    core::sort::reverse(mys.begin(),mys.end());
		    LOG(BF("Done reverse") );
		    // print "atom(%s) mys = %s"%(myatom.getName(),str(mys))
		    this->_s[myatom->getRelativePriority()] = mys;
		}
	    }
#ifdef	DEBUG_ON
	    { _BLOCK_TRACEF(BF("STAGE2 results - there are %d classes")% C.size() );
    //	vector< vector<int> >::iterator sit;
                gctools::Vec0<Atom_sp>::iterator ait;
	    LOG(BF("Reverse sorted neighbor priorities for each atom") );
	    for ( ait=mAtoms.begin(); ait!=mAtoms.end(); ait++ )
	    {
		vector<int>::iterator zit;
		stringstream ss;
		for ( zit=this->_s[(*ait)->getRelativePriority()].begin();
			    zit!=this->_s[(*ait)->getRelativePriority()].end();
			    zit++ )
		{
		    ss << " " << *zit;
		}
		LOG(BF("  reverse sorted priorities of neighbors of %s - %s") % (*ait)->getName().c_str() % ss.str().c_str()  );
	    }
	}
#endif
	// 3.	For each class Cr, partition the atoms in the class into 
	//  ordered subclasses (S1,...,Sk) such that for each atom i in subclass 
	//  Sr and each atom j in subclass Ss we have (i) s(i)<s(j) iff r<s; 
	//  and (ii) s(i)=s(j) iff r=s. (Note, the s(i) strings are compared 
	//  lexicographically.)


	// didPartition = True
	// while ( didPartition ):
	    newC.clear();
	    {_BLOCK_TRACE("Stage 3");
		didPartition = false;
		uint classIndex;
		for ( classIndex = 0; classIndex < C.size(); classIndex++ )
		{ _BLOCK_TRACEF(BF(" classIndex = %u/%u") % classIndex%C.size() );
		    uint classBegin = C[classIndex];
		    uint classEnd;
		    if ( classIndex == C.size()-1 )
		    {
			classEnd = mAtoms.size();
		    } else
		    {
			classEnd = C[classIndex+1];
		    }
		    LOG(BF("Looking at class (%u)-(%u)") % classBegin % classEnd  );
                    gctools::Vec0<Atom_sp> S;
		    {_BLOCK_TRACE("Filling S");
			for ( uint ci = classBegin; ci<classEnd; ci++ )
			{
			    Atom_sp a = mAtoms[ci];
			    S.push_back(a);
			    // print "Before sort S.atom(%s) s=%s"%(a.getName(),str(s[a.getRelativePriority()]))
			}
		    }
		    OrderByS byS;
		    {_BLOCK_TRACE("Sorting mAtoms and S");
			LOG(BF("Setting up prioritizer") );
			byS.setCipPrioritizer(this->sharedThis<CipPrioritizer_O>());
			LOG(BF("about to sort number of elements = %d") % S.size()  );
#ifdef DEBUG_ON
			LOG(BF("Contents of S") );
			for ( gctools::Vec0<Atom_sp>::iterator ssi=S.begin();ssi!=S.end();ssi++)
			{
			    LOG(BF("    %s") % (*ssi)->description().c_str() );
			}
#endif
			core::sort::quickSort(S.begin(),S.end(),byS);
			LOG(BF("done sort") );
		    }
//		    int i = 0;
		    LOG(BF("About to dump sort results") );
#ifdef	DEBUG_ON
                    gctools::Vec0<Atom_sp>::iterator iiS;
for ( iiS=S.begin(); iiS!=S.end(); iiS++)
{
    stringstream sz;
    sz << "After sort S[" << iiS-S.begin() << "].atom("<<(*iiS)->getName();
    sz << " neighbors: ";
    vector<int>::iterator zit;
    for ( zit=this->_s[(*iiS)->getRelativePriority()].begin();
		zit!=this->_s[(*iiS)->getRelativePriority()].end();
		zit++ )
    {
	sz << " " << *zit;
    }

    LOG(BF("%s") % sz.str().c_str()  );
}
#endif
		//    for a in S:
		// 	print "After sort S[%d].atom(%s) s=%s"%(classBegin+i,a.getName(),str(s[a.getRelativePriority()]))
		// 	i += 1
		    LOG(BF("Replacing mAtoms from (%u)-(%u) with S") % classBegin % classEnd );
		    uint si, mi;
		    for ( si=0, mi=classBegin; mi!=classEnd; mi++, si++ ) mAtoms[mi] = S[si];
		    LOG(BF("Partitioning based on s") );
		    vector<int> partitionS = partition(S,classBegin,byS);
#if DEBUG_ON
		    vector<int>::iterator pSi;
		    stringstream Ss;
		    for ( pSi=partitionS.begin(); pSi!=partitionS.end(); pSi++)
		    {
			Ss << " " << *pSi;
		    }
		    LOG(BF(" partitionS = %s") % Ss.str().c_str()  );
#endif
		    if ( partitionS.size()>1 ) 
		    {
			LOG(BF(" Did partition of S") );
			didPartition = true;
		    }
		    for ( vector<int>::iterator psi=partitionS.begin(); psi!=partitionS.end(); psi++ )
		    {
			newC.push_back(*psi);
		    }
		}
	    }
	// 4.	If every class was partitioned into only one subclass then 
	//  terminate with p(i) as the priority of atom i.
	    { _BLOCK_TRACE("Stage 4");
		if ( !didPartition ) break;
	    }

	    // print "Repartitioning"
	// 5.	Form a new partition of all the atoms by concatenating 
	//  	all of the computed subclasses of all of the classes (in the 
	//  	same sequence as the original classes).

	    { _BLOCK_TRACEF(BF("Stage 5 - there are now %d classes") % newC.size());
		C = newC;
	    }
	    // print "Regenerated C = ", C

	// 6.	For each class Cr and for each atom i in Cr set 
	// 	p(i) to r and go to Step 3.
	// 
	    // print "Resetting p's"
	    { _BLOCK_TRACE("Stage 6" );
		uint classBegin, classEnd;
		for ( uint classIndex= 0; classIndex < C.size(); classIndex++ )
		{
		    classBegin = C[classIndex];
		    if ( classIndex == C.size()-1 )
		    {
			classEnd = mAtoms.size();
		    } else
		    {
			classEnd = C[classIndex+1];
		    }
		    for ( uint cci =classBegin; cci!=classEnd; cci++ )
		    {
			Atom_sp a = mAtoms[cci];
			this->_p[a->getRelativePriority()] = classIndex;
		    }
		}
	    }
	}
    }
    // The partitioning steps can be effected with sorting and since 
    // all other steps require linear time, we have that each iteration 
    // of the algorithm requires O(nlogn) time (assuming bounded degree 
    // of all n atoms). At most n iterations are required giving a total 
    // running time of O(n2logn).
    // 
    // Once CIP priorities have been assigned to every atom in a 
    //  molecule it is a simple matter to order the neighbors of 
    //  each atom and compute the appropriate signed volume tests 
    //  on the atomic coordinates to make topological R, S, cis, 
    //  or trans assignments.
    // 
    // Once the basic chiral assignments have been made, the 
    // classes can be further partitioned taking the initial 
    // assignment into account. This process will create new 
    // chiral centers based on the chirality of the branches.
    //  now write the CIP relative priority into RelativePriority 
    { _BLOCK_TRACE("Writing priorities to atoms");
	for ( gctools::Vec0<Atom_sp>::iterator ai=mAtoms.begin(); ai!=mAtoms.end(); ai++ )
	{
	    uint relPriority = this->_p[(*ai)->getRelativePriority()];
	    (*ai)->setRelativePriority(relPriority);
	    LOG(BF("Assigned to atom: %s priority: %d") % (*ai)->getName().c_str() % relPriority  );
	}
    }
}


bool orderByPriority( Atom_sp p1, Atom_sp p2 )
{
    return p2->getRelativePriority() < p1->getRelativePriority();
}

bool orderByName( Atom_sp p1, Atom_sp p2 )
{
  return core::SymbolComparer::order(p1->getName(),p2->getName())<=0;
}



void CipPrioritizer_O::defineStereochemicalConfigurationsForAllAtoms(Matter_sp molOrAgg)
{_G();
    this->assignCahnIngoldPrelogPriorityToAtomsRelativePriority(molOrAgg);
    Loop l;
    l.loopTopGoal(molOrAgg,ATOMS);
    while ( l.advanceLoopAndProcess() )
    {
	Atom_sp a = l.getAtom();
	if ( a.isA<VirtualAtom_O>() ) continue;
	if ( a->numberOfBonds() == 4 )
	{
	    //  get the names of the atoms in order of priority
            gctools::Vec0<Atom_sp> priority;
	    for ( gctools::Vec0<Bond_sp>::iterator bi = a->bonds_begin(); bi!=a->bonds_end(); bi++ )
	    {
		priority.push_back((*bi)->getOtherAtom(a));
	    }
	    core::sort::quickSort(priority.begin(), priority.end(), orderByPriority);

	    //  now figure out if they are four different priorities
	    uint prevPriority = UndefinedUnsignedInt;
	    int diff = 0;
	    for ( gctools::Vec0<Atom_sp>::iterator pi=priority.begin(); 
	    		pi!=priority.end(); pi++ )
	    {
		if ( (*pi)->getRelativePriority() != prevPriority )
		{
		    diff += 1;
		}
		prevPriority = (*pi)->getRelativePriority();
	    }
	    if ( diff == 4 )
	    {
		//  they are all different, this is a chiral atom
		// If the stereochemistry type has already been assigned then don't overwrite it.
		if ( a->getStereochemistryType() != chiralCenter )
		{
		    a->setStereochemistryType(chiralCenter);
		    a->setConfiguration(undefinedConfiguration);
		}
	    } else 
	    {
		//  there are less than 4 different priorities that means
		//  we treat this as a proChiral center, methyls will also order
		//  their atoms
		a->setStereochemistryType(prochiralCenter);
		a->setConfiguration(undefinedConfiguration);
	    }
	} else
	{
	    a->setConfiguration(undefinedConfiguration);
	    a->setStereochemistryType(undefinedCenter);
	}
    }
}




    void CipPrioritizer_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<CipPrioritizer_O>("makeCipPrioritizer")
	.def("assignCahnIngoldPrelogPriorityToAtomsRelativePriority",
	     &CipPrioritizer_O::assignCahnIngoldPrelogPriorityToAtomsRelativePriority)
    .def("defineStereochemicalConfigurationsForAllAtoms",
    		&CipPrioritizer_O::defineStereochemicalConfigurationsForAllAtoms )
    ;

}
    void CipPrioritizer_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,CipPrioritizer,"","",_lisp)
    ;
#endif //]

}


    EXPOSE_CLASS(chem,CipPrioritizer_O);
};



