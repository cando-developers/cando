/*
    File: cipPrioritizer.cc
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
#define	DEBUG_LEVEL_FULL

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
#include <clasp/core/hashTableEq.h>
#include <cando/chem/bond.h>
#include <cando/chem/virtualAtom.h>
#include <clasp/core/wrappers.h>



namespace chem {





uint atomicNumber(Atom_sp a)
{
  return a->getAtomicNumber();
}

struct	OrderByBase
{
  OrderByBase() {};
  CipPrioritizer_sp	prior;
  void setCipPrioritizer(CipPrioritizer_sp p) { this->prior = p;};
  virtual bool operator() (AtomPriority& x, AtomPriority& y) = 0;
};

struct OrderByPriority : public OrderByBase {
public:
  bool operator()(AtomPriority& x, AtomPriority& y )
  {
    return y._relativePriority < x._relativePriority;
  }
  OrderByPriority() : OrderByBase() {};
};

struct	OrderByP : public OrderByBase
{
public:
  bool operator()(AtomPriority& x, AtomPriority& y )
  {
    if ( prior->getP(x) <= prior->getP(y) ) return true;
    return false;
  }
  OrderByP() : OrderByBase() {};
};


struct	OrderByS : public OrderByBase
{
  OrderByS() {};
  int cmpByS(AtomPriority& x, AtomPriority& y )
  {
    ASSERTNOTNULL(x);
    ASSERTNOTNULL(y);
    gctools::Vec0<int>&	xv = prior->getS(x);
    gctools::Vec0<int>&	yv = prior->getS(y);
    gctools::Vec0<int>::iterator	xi = xv.begin();
    gctools::Vec0<int>::iterator	yi = yv.begin();
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


  bool operator()(AtomPriority& x, AtomPriority& y )
  {
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

int	CipPrioritizer_O::getP(AtomPriority& a)
{
  ASSERT_lessThan(a._relativePriority,this->_p.size());
  return this->_p[a._relativePriority];
}

gctools::Vec0<int>& CipPrioritizer_O::getS(AtomPriority& a)
{
  ASSERTNOTNULL(a);
  if ( ! ( a._relativePriority < this->_s.size() ) )
  {
    LOG("Bad priority for atom" );
    LOG("Bad priority for atom: {}" , a._Atom->description().c_str()  );
    LOG("   priority value = {}" , a._relativePriority  );
  }
  ASSERT_lessThan(a._relativePriority,this->_s.size());
  return this->_s[a._relativePriority];
}



/*!
 * partition the atoms between cBegin and cEnd-1 into classes based on orderFn
 */
vector<int>	partition(gctools::Vec0<AtomPriority> atoms,int indexOffset, OrderByS& orderer)
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


CL_LISPIFY_NAME(chem:assign-priorities-hash-table);
DOCGROUP(cando);
CL_DEFUN core::HashTable_sp CipPrioritizer_O::assignPrioritiesHashTable(Matter_sp matter)
{
  core::HashTable_sp cip = core::HashTableEq_O::create_default();
  CipPrioritizer_sp prior;
  if (gc::IsA<Molecule_sp>(matter)) {
    prior = CipPrioritizer_O::create();
    prior->assignCahnIngoldPrelogPriorityToAtomsRelativePriority(matter,cip);
  } else if (gc::IsA<Aggregate_sp>(matter)) {
    Loop l;
    l.loopTopGoal(matter,MOLECULES);
    while ( l.advanceLoopAndProcess() ) {
      Molecule_sp mol = l.getMolecule();
      prior = CipPrioritizer_O::create();
      prior->assignCahnIngoldPrelogPriorityToAtomsRelativePriority(mol,cip);
    }
  } else {
    TYPE_ERROR(matter,core::Cons_O::createList(_sym_Aggregate_O,_sym_Molecule_O));
  }
  return cip;
}


//
// ! Assign the CIP priority to the RelativePriority of each Atom in the molecule
CL_DOCSTRING(R"dx(
Assign the topological RS chirality
see: http://www.chemcomp.com/journal/chiral.htm  (dead link now - check Wayback Machine)

An Efficient Algorithm For The Determination Of Topological RS Chirality
P. Labute
Chemical Computing Group Inc.
1010 Sherbrooke Street W, Suite 910; Montreal, Quebec; Canada H3A 2R7.
November 20, 1996

Determination of chiral centers and labeling by the RS system as proposed by Cohen, Ingold and Prelog is important in computational chemistry for a number of reasons. Aside the usual uses of nomenclature, specification of chiral constraints during energy minimization and detection of symmetries, proper labeling of chirality is extremely useful for 2D renderings of conformations. A number of algorithms have been proposed. Many are incorrect or do not correspond to the CIP system. Others are computationally expensive.

We present an efficient algorithm for the assignment of CIP priorities to every atom in a molecule. This information can then be used to assign R, S, cis or trans labels.

)dx")
CL_LISPIFY_NAME("assignCahnIngoldPrelogPriorityToAtomsRelativePriority");
CL_DEFMETHOD void CipPrioritizer_O::assignCahnIngoldPrelogPriorityToAtomsRelativePriority(Matter_sp molOrAgg, core::HashTable_sp cip)
{
  gctools::Vec0<AtomPriority>	mAtoms;
  vector<int> newC;
  vector<int> C;
    // 
    //  Determination of topological RS chirality
    //  see: http://www.chemcomp.com/journal/chiral.htm (link dead see Wayback machine) paper included below
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
    // C.append(len(mAtoms))

    // print "Classes go from C[i] to C[i+1] unless last i then end, C=",C

  { 
    this->_p.clear();
    Loop l;
    l.loopTopGoal(molOrAgg,ATOMS);
    while ( l.advanceLoopAndProcess() )
    {
      Atom_sp a = l.getAtom();
      if ( a.isA<VirtualAtom_O>() ) continue;
      size_t relativePriority = mAtoms.size();
      cip->setf_gethash(a,core::clasp_make_fixnum(relativePriority));
      AtomPriority ap(a,relativePriority);
      mAtoms.push_back(ap);
      this->_p.push_back(atomicNumber(a)); // p(i) = atomic_number of atom
      LOG("Atom {} has atomicNumber() {}\n", _rep_(a), atomicNumber(a));
    }
    OrderByP byP;
    byP.setCipPrioritizer(this->sharedThis<CipPrioritizer_O>());
#ifdef	DEBUG_ON
    LOG("Assigning initial priorities STAGE1" );
    for ( auto mit=mAtoms.begin(); mit!=mAtoms.end(); mit++ ) {
      LOG("\nAtom: {} _relativePriority {} p-> {}\n", _rep_((*mit)._Atom), (*mit)._relativePriority, this->getP(*mit));
    }
    for ( auto cidx=0; cidx<C.size(); cidx++ ) {
      LOG("\nC[{}] = {}\n", cidx, C[cidx]);
    }
#endif
    if ( mAtoms.size()>1 ) {
      sort::quickSortVec0(mAtoms, 0, mAtoms.size(), byP);
    }

    C.clear();
    int pcur = -1;
    for ( uint ai=0; ai<mAtoms.size(); ai++ ) {
      Atom_sp a = mAtoms[ai]._Atom;
      uint idx = mAtoms[ai]._relativePriority;
      if (pcur != this->_p[idx])
      {
        C.push_back(ai);
        pcur = this->_p[idx];
      }
    }
#ifdef	DEBUG_ON
    LOG("Assigning priorities STAGE1" );
    for ( auto mit=mAtoms.begin(); mit!=mAtoms.end(); mit++ ) {
      LOG("\nAtom: {} _relativePriority {} p-> {}\n", _rep_((*mit)._Atom), (*mit)._relativePriority, this->getP(*mit));
    }
    for ( auto cidx=0; cidx<C.size(); cidx++ ) {
      LOG("\nC[{}] = {}\n", cidx, C[cidx]);
    }
#endif
  }

    // 2.	For each atom i set s(i)=abc...z to be an ordered list of 
    //  neighboring p(j) numbers in decreasing order accounting for bond 
    //  multiplicities with repeated values (i.e., if atom i is double 
    //  bonded to an atom with priority p, then put p in the list twice).

  {
		    // previously this was at the top of step 3
    bool didPartition = true;
    while ( didPartition )
    {
      this->_s.clear();
      this->_s.resize(mAtoms.size());
      { 
        for ( gctools::Vec0<AtomPriority>::iterator mi =mAtoms.begin();
              mi!=mAtoms.end(); mi++ )
        {
          Atom_sp myatom = (*mi)._Atom;
          LOG("About to fill mys" );
          gctools::Vec0<int>	mys;
          for ( gctools::Vec0<Bond_sp>::iterator bi=myatom->bonds_begin();
                bi!=myatom->bonds_end(); bi++ )
          {
            int bondOrder = 1;
            if ( Bond_O::singleBondP((*bi)->getRawOrder()) ) bondOrder = 1;
            if ( (*bi)->getRawOrder() == doubleBond ) bondOrder = 2;
            if ( (*bi)->getRawOrder() == tripleBond ) bondOrder = 3;
            for ( int oi=0; oi<bondOrder; oi++ )
            {
              mys.push_back(this->_p[(*bi)->getOtherAtom(myatom)->getRelativePriority(cip)]);
            }
          }
          LOG("About to sort {} mys objects" , mys.size()  );
          if (mys.size()>1) {
            int* begin = &mys[0];
            int* end = &mys[mys.size()];
            sort::quickSortVec0(mys, 0, mys.size());
            LOG("Done sort" );
            sort::reverse(begin,end);
            LOG("Done reverse" );
          }
          this->_s[myatom->getRelativePriority(cip)] = mys;
        }
      }
#ifdef	DEBUG_ON
      { 
        LOG("Reverse sorted neighbor priorities for each atom" );
        for ( auto ait=mAtoms.begin(); ait!=mAtoms.end(); ait++ )
        {
          stringstream ss;
          for ( auto zit=this->_s[(*ait)._relativePriority].begin();
                zit!=this->_s[(*ait)._relativePriority].end();
                zit++ )
          {
            ss << " " << *zit;
          }
          LOG("  reverse sorted priorities of neighbors of {} - {}" , _rep_((*ait)._Atom), ss.str().c_str()  );
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
      {
        didPartition = false;
        uint classIndex;
        for ( classIndex = 0; classIndex < C.size(); classIndex++ )
        { 
          uint classBegin = C[classIndex];
          uint classEnd;
          if ( classIndex == C.size()-1 )
          {
            classEnd = mAtoms.size();
          } else
          {
            classEnd = C[classIndex+1];
          }
          LOG("Looking at class ({})-({})" , classBegin , classEnd  );
          gctools::Vec0<AtomPriority> S;
          {
            for ( uint ci = classBegin; ci<classEnd; ci++ )
            {
              S.push_back(mAtoms[ci]);
            }
          }
          OrderByS byS;
          {
            LOG("Setting up prioritizer" );
            byS.setCipPrioritizer(this->sharedThis<CipPrioritizer_O>());
            LOG("about to sort number of elements = {}" , S.size()  );
#ifdef DEBUG_ON
            LOG("Contents of S" );
            for ( auto ssi=S.begin();ssi!=S.end();ssi++)
            {
              LOG("    {} relativePriority {}" , _rep_((*ssi)._Atom), (*ssi)._relativePriority );
            }
#endif
            if (S.size()>1) sort::quickSortVec0(S, 0, S.size(), byS);
            LOG("done sort" );
          }
//		    int i = 0;
          LOG("About to dump sort results" );
#ifdef	DEBUG_ON
          for ( auto iiS=S.begin(); iiS!=S.end(); iiS++)
          {
            stringstream sz;
            sz << "After sort S[" << iiS-S.begin() << "].atom("<<_rep_((*iiS)._Atom);
            sz << " neighbors: ";
            for ( auto zit=this->_s[(*iiS)._relativePriority].begin();
                  zit!=this->_s[(*iiS)._relativePriority].end();
                  zit++ )
            {
              sz << " " << *zit;
            }

            LOG("{}" , sz.str().c_str()  );
          }
#endif
          LOG("Replacing mAtoms from ({})-({}) with S" , classBegin , classEnd );
          uint si, mi;
          for ( si=0, mi=classBegin; mi!=classEnd; mi++, si++ ) {
            mAtoms[mi]._Atom = S[si]._Atom;
            mAtoms[mi]._relativePriority = S[si]._relativePriority;
          }
          LOG("Partitioning based on s" );
          vector<int> partitionS = partition(S,classBegin,byS);
#if DEBUG_ON
          vector<int>::iterator pSi;
          stringstream Ss;
          for ( pSi=partitionS.begin(); pSi!=partitionS.end(); pSi++)
          {
            Ss << " " << *pSi;
          }
          LOG(" partitionS = {}" , Ss.str().c_str()  );
#endif
          if ( partitionS.size()>1 ) 
          {
            LOG(" Did partition of S" );
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
      { 
        if ( !didPartition ) break;
      }

	    // print "Repartitioning"
	// 5.	Form a new partition of all the atoms by concatenating 
	//  	all of the computed subclasses of all of the classes (in the 
	//  	same sequence as the original classes).

      { 
        C = newC;
      }
	    // print "Regenerated C = ", C

	// 6.	For each class Cr and for each atom i in Cr set 
	// 	p(i) to r and go to Step 3.
	// 
	    // print "Resetting p's"
      { 
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
            AtomPriority& a = mAtoms[cci];
            this->_p[a._relativePriority] = classIndex;
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
  { 
    for ( gctools::Vec0<AtomPriority>::iterator ai=mAtoms.begin(); ai!=mAtoms.end(); ai++ )
    {
      uint relPriority = this->_p[(*ai)._relativePriority];
      cip->setf_gethash((*ai)._Atom,core::clasp_make_fixnum(relPriority));
      LOG("Assigned to atom: {} priority: {}" , _rep_((*ai)._Atom->getName()) , relPriority  );
    }
  }
}


bool orderByName( Atom_sp p1, Atom_sp p2 )
{
  return core::SymbolComparer::order(p1->getName(),p2->getName())<=0;
}

CL_DOCSTRING("Return (values atom-to-cip-relative-priorities-map atom-to-stereochemistry-type)")
CL_DOCSTRING_LONG(R"dx(For each atom in the matter calculate the relative CIP priority and the stereochemistry type.
Return these two values in a pair of hash-tables each keyed on the atom)dx")
CL_LISPIFY_NAME("calculateStereochemicalTypeForAllAtoms");
CL_DEFMETHOD core::HashTable_mv CipPrioritizer_O::calculateStereochemistryTypeForAllAtoms( Matter_sp molOrAgg )
{
  core::HashTable_sp cip = core::HashTableEq_O::create_default();
  core::HashTable_sp stereochemistryType = core::HashTableEq_O::create_default();
  this->assignCahnIngoldPrelogPriorityToAtomsRelativePriority( molOrAgg, cip );
  Loop l;
  l.loopTopGoal( molOrAgg, ATOMS );
  while ( l.advanceLoopAndProcess() ) {
    Atom_sp a = l.getAtom();
    if ( a.isA<VirtualAtom_O>() ) continue;
    if ( a->numberOfBonds() == 4 ) {
	    //  get the names of the atoms in order of priority
      gctools::Vec0<AtomPriority> priority;
      for ( gctools::Vec0<Bond_sp>::iterator bi = a->bonds_begin(); bi!=a->bonds_end(); bi++ ) {
        AtomPriority ap( (*bi)->getOtherAtom(a), a->getRelativePriority(cip) );
        priority.push_back(ap);
      }
      OrderByPriority orderByPriority;
      sort::quickSortVec0( priority, 0, priority.size(), orderByPriority );
	    //  now figure out if they are four different priorities
      uint prevPriority = UndefinedUnsignedInt;
      int diff = 0;
      for ( gctools::Vec0<AtomPriority>::iterator pi=priority.begin(); pi!=priority.end(); pi++ ) {
        if ( (*pi)._relativePriority != prevPriority ) {
          diff += 1;
        }
        prevPriority = (*pi)._relativePriority;
      }
      if ( diff == 4 ) {
		//  they are all different, this is a chiral atom
		// If the stereochemistry type has already been assigned then don't overwrite it.
        if ( a->getStereochemistryType() != chiralCenter ) {
          stereochemistryType->setf_gethash( a, chemkw::_sym_chiral );
        }
      } else {
		//  there are less than 4 different priorities that means
		//  we treat this as a proChiral center, methyls will also order
		//  their atoms
        stereochemistryType->setf_gethash( a, chemkw::_sym_prochiral );
      }
    } else {
      stereochemistryType->setf_gethash( a, chemkw::_sym_undefinedCenter );
    }
  }
  return Values(cip, stereochemistryType );
}


void doCalculateStereochemistry( bool useStructure, Molecule_sp matter, core::HashTable_sp cips, core::HashTable_sp stereochemistryTypes, core::HashTable_sp configurations ) {
  CipPrioritizer_sp prior;
  prior = CipPrioritizer_O::create();
  core::HashTable_mv cip = prior->calculateStereochemistryTypeForAllAtoms(matter);
  core::MultipleValues &values = core::lisp_multipleValues();
  core::T_sp ret2 = values.second( cip.number_of_values());
  if (ret2.nilp()) {
    SIMPLE_ERROR("NIL was returned rather than a hash-table for stereochemistryType");
  }
  core::HashTable_sp stereochemistryType = gc::As<core::HashTable_sp>(ret2);
  Loop la;
  la.loopTopGoal(matter,ATOMS);
  while (la.advanceLoopAndProcess() ) {
    Atom_sp atm = la.getAtom();
    core::T_sp st = stereochemistryType->gethash(atm);
    if ( st.notnilp() ) {
      ConfigurationEnum ce(undefinedConfiguration);
      if (useStructure) ce = atm->calculateStereochemicalConfiguration(cip);
      else stereochemistryTypes->setf_gethash( atm, st );
      configurations->setf_gethash( atm, translate::to_object<ConfigurationEnum>::convert(ce) );
      cips->setf_gethash( atm, cip->gethash(atm) );
    }
  }
}


CL_DOCSTRING(R"dx(Calculate the stereochemistry for each atom in the aggregate or molecule)dx");
CL_DOCSTRING_LONG(R"dx(If use-structure is T then use the three-dimensional structure to calculate the configuration of each atom.
Otherwise pull the configuration out of the atoms _Configuration slot.
Return (values atom-to-stereochemistry-type-hash-table atom-to-configuration-hash-table atom-to-relative-cip-hash-table) )dx")
CL_LAMBDA(matter &key use-structure)
DOCGROUP(cando);
CL_DEFUN core::HashTable_mv chem__calculateStereochemistry( Matter_sp matter, bool useStructure ) {
  core::HashTableEq_sp stereochemistryTypes = core::HashTableEq_O::create_default();
  core::HashTableEq_sp configurations = core::HashTableEq_O::create_default();
  core::HashTableEq_sp cips = core::HashTableEq_O::create_default();
  if (gc::IsA<Molecule_sp>(matter)) {
    doCalculateStereochemistry( useStructure, gc::As_unsafe<Molecule_sp>( matter ), cips, stereochemistryTypes, configurations );
  } else if (gc::IsA<Aggregate_sp>(matter)) {
    Loop l;
    l.loopTopGoal(matter,MOLECULES);
    while ( l.advanceLoopAndProcess() ) {
      Molecule_sp mol = l.getMolecule();
      doCalculateStereochemistry( useStructure, mol, cips, stereochemistryTypes, configurations );
    }
  } else {
    TYPE_ERROR(matter,core::Cons_O::createList(_sym_Aggregate_O,_sym_Molecule_O));
  }
  return Values( stereochemistryTypes, configurations, cips );
}










};

#if 0
An Efficient Algorithm For The Determination Of Topological RS Chirality


P. Labute
Chemical Computing Group Inc.
1010 Sherbrooke Street W, Suite 910; Montreal, Quebec; Canada H3A 2R7.


November 20, 1996


Abstract. In this work, we describe an efficient method to detect and label chiral centers as per the RS system.
INTRODUCTION


Determination of chiral centers and labeling by the RS system as proposed by Cohen, Ingold and Prelog is important in computational chemistry for a number of reasons. Aside the usual uses of nomenclature, specification of chiral constraints during energy minimization and detection of symmetries, proper labeling of chirality is extremely useful for 2D renderings of conformations. A number of algorithms have been proposed. Many are incorrect or do not correspond to the CIP system. Others are computationally expensive.

We present an efficient algorithm for the assignment of CIP priorities to every atom in a molecule. This information can then be used to assign R, S, cis or trans labels.

METHOD


The objective of the algorithm is to assign a priority, p(i) to every atom i in the molecule so that if p(i)<p(j) then atom i has strictly lower CIP priority than atom j with equality only occurring when the atoms are indistinguishable (topologically). In essence, the algorithm maintains an ordered list of equivalence classes of atoms. Each atom in an equivalence class is assigned the priority of the -1 of the class in the sorted list. The algorithm repeatedly splits classes (maintaining the ordering) until no changes occur. The final priority assignments are then output as the CIP priorities.

The algorithm proceeds as follows:

For each atom i, set p(i) equal to the atomic number of atom i. Set the initial partition to be an ordered list of classes (C1,…,Ck) such that for each atom i in class Cr and each atom j in class Cs we have (i) p(i)<p(j) iff r<s; and (ii) p(i)=p(j) iff r=s.
For each atom i set s(i)=abc…z to be an ordered list of neighboring p(j) numbers in decreasing order accounting for bond multiplicities with repeated values (i.e., if atom i is double bonded to an atom with priority p, then put p in the list twice).
For each class Cr, partition the atoms in the class into ordered subclasses (S1,…,Sk) such that for each atom i in subclass Sr and each atom j in subclass Ss we have (i) s(i)<s(j) iff r<s; and (ii) s(i)=s(j) iff r=s. (Note, the s(i) strings are compared lexicographically.)
If every class was partitioned into only one subclass then terminate with p(i) as the priority of atom i.
Form a new partition of all the atoms by concatenating all of the computed subclasses of all of the classes (in the same sequence as the original classes).
For each class Cr and for each atom i in Cr set p(i) to r and go to Step 3.
The partitioning steps can be effected with sorting and since all other steps require linear time, we have that each iteration of the algorithm requires O(nlogn) time (assuming bounded degree of all n atoms). At most n iterations are required giving a total running time of O(n2logn).

Once CIP priorities have been assigned to every atom in a molecule it is a simple matter to order the neighbors of each atom and compute the appropriate signed volume tests on the atomic coordinates to make topological R, S, cis, or trans assignments.

Once the basic chiral assignments have been made, the classes can be further partitioned taking the initial assignment into account. This process will create new chiral centers based on the chirality of the branches.

DISCUSSION


The RS system uses bond orders which, from a computational standpoint, is unfortunate since (i) a particular resonance structure must be chosen; (ii) the "phantom" atoms, although easy to deal with, are an unnecessary complication.

As an alternative we propose the following system. The algorithm to compute assignments is the same as the one presented except for the initial assignment of priorities. With CIP priorities, the initial string is made up of the neighboring atomic numbers taking multiplicities into account. The new system need only replace the initial priority assignment of atomic number with a code taking into account further properties.

Initially each atom is assigned a code of the following form:

(a128+s)8+b
where a is the atomic number of the atom, s is the isotope number, b is the hybridization of the atom coded as sp=6, sp2=5, d3sp3=4, d2sp3=3, dsp3=2, sp3=1, and the ground state being coded as 0.

This system has the advantage that the hybridization state is invariant under resonance, is easy to detect, assignable in the absence of hydrogens, and will still be close to the original CIP priorities.

#endif
