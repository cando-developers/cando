/*
    File: topology.cc
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


//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/topology.h>
#include <cando/chem/coupling.h>
#include <cando/chem/monomer.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <clasp/core/environment.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <clasp/core/wrappers.h>

namespace chem
{

void TopologyAtomInfo_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,fftype),this->_fftype);
  node->field(INTERN_(kw,charge),this->_AtomicCharge);
  node->field(INTERN_(kw,index),this->_ConstitutionAtomIndex);
  this->Base::fields(node);
}


void Topology_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,name),this->_Name);
  node->field(INTERN_(kw,constitution),this->_Constitution);
  node->field(INTERN_(kw,notrainers),this->_SuppressTrainers);
  node->field(INTERN_(kw,atominfo),this->_AtomInfo);
  node->field(INTERN_(kw,atomprops),this->_StereoisomerAtomProperties);
  node->field(INTERN_(kw,properties),this->_Properties);
  node->field(INTERN_(kw,defaultoutplug),this->_DefaultOutPlugName);
  node->field(INTERN_(kw,plugs),this->_Plugs);
  this->Base::fields(node);
}


Topology_mv Topology_O::makeTopologyFromResidue(chem::Residue_sp residue, core::Symbol_sp topologyName,  core::T_sp tconstitution )
{
  Constitution_sp constitution;
  if ( constitution.nilp() ) {
    SIMPLE_ERROR(BF("The constitution was nil - you need to pass a constitution"));
  }
  constitution = gc::As<Constitution_sp>(tconstitution);
  ConstitutionAtoms_sp ca = constitution->getConstitutionAtoms();
  Topology_sp topology = Topology_O::make(topologyName, constitution, _Nil<core::T_O>() );
  core::VectorTNs_sp vec = core::VectorTNs_O::make(residue->numberOfAtoms(),_Nil<core::T_O>());
  int atomIndex = 0;
  for ( auto ai = residue->begin_atoms(); ai!=residue->end_atoms(); ++ai, ++atomIndex) {
    Atom_sp atom = *ai;
    ConstitutionAtomIndex0N constitutionIndex = constitution->getConstitutionAtoms()->index(atom->getName());
    TopologyAtomInfo_sp tai = TopologyAtomInfo_O::make(atom->getType(),
                                                       atom->getCharge(),
                                                       constitutionIndex);
    vec->rowMajorAset(atomIndex,tai);
  }
  topology->setf_atomInfo(vec);
  return Values(topology,constitution);
}

Topology_sp Topology_O::make(core::Symbol_sp name, Constitution_sp constitution, core::List_sp plugs)
{
  Topology_sp me = gctools::GC<Topology_O>::allocate(name,constitution);
    //
    // If "plugs" option is outgoing plugs only then we should
    // change the name to "outPlugs"
    //
  me->_Plugs.clear();
  core::List_sp curPlug;
  for ( auto curPlug : plugs ) {
    auto p = core::oCar(curPlug).as<Plug_O>();
    me->addPlug(p->getName(),p);
  }
  return me;
};


CL_DEFMETHOD Residue_sp Topology_O::build_residue() const
{
  Residue_sp res = Residue_O::make(this->getName());
  gctools::Vec0<Atom_sp> atoms;
  atoms.resize(this->_AtomInfo->arrayTotalSize());
  res->resizeContents(this->_AtomInfo->arrayTotalSize());
  ConstitutionAtoms_sp constitutionAtoms = this->_Constitution->getConstitutionAtoms();
  size_t idx = 0;
  for ( size_t idx=0, idxEnd(this->_AtomInfo->arrayTotalSize()); idx<idxEnd; ++idx ) {
    TopologyAtomInfo_sp ai = gc::As<TopologyAtomInfo_sp>(this->_AtomInfo->rowMajorAref(idx));
    Atom_sp atom = Atom_O::create();
    ConstitutionAtom_sp ca = (*constitutionAtoms)[ai->constitutionAtomIndex()];
    atom->setName(ca->_AtomName);
    atom->setElement(ca->_Element);
    atom->setType(ai->_fftype);
    atom->setCharge(ai->_AtomicCharge);
    atom->turnOnFlags(needsBuild);
    atoms[ai->_ConstitutionAtomIndex] = atom;
    res->putMatter(idx,atom); // order atoms as in Topology
  }
  for ( size_t i=0, iEnd(constitutionAtoms->numberOfAtoms()); i<iEnd; ++i ) {
    Atom_sp fromAtom = atoms[i];
    ConstitutionAtom_sp ca = (*constitutionAtoms)[i];
    for ( auto bi=ca->_Bonds.begin(); bi!=ca->_Bonds.end(); ++bi )
    {
      Atom_sp toAtom = atoms[(*bi)->_ToAtomIndex];
      if ( fromAtom->atLowerUniqueAtomOrderThan(toAtom) ) {
        BondOrder order = (*bi)->_BondOrder;
        fromAtom->bondTo(toAtom,order);
      }
    }
  }
  return res;
}



#if 0
    void __checkExtractFragment(ExtractFragmentBase_sp extractFrag,adapt::StringSet_sp missingAtomNames, map<string,OverlappingFragments>& overlaps)
    {

	Fragment_sp frag = extractFrag->getFragment();
	adapt::StringSet_sp atomNames = frag->getAtomNames();
	for ( adapt::StringSet_O::iterator ni = atomNames->begin(); ni!=atomNames->end(); ni++ )
	{
	    if ( missingAtomNames->contains(*ni) )
	    {
		missingAtomNames->remove(*ni);
	    }
	    if ( overlaps.count(*ni)==0 )
	    {
		OverlappingFragments over;
		over._TimesSeen = 1;
		over._Fragments.insert(frag->getName());
		overlaps[*ni] = over;
	    } else
	    {
		OverlappingFragments& over = overlaps[*ni];
		over._TimesSeen = over._TimesSeen + 1;
		over._Fragments.insert(frag->getName());
	    }
	}
    }

    void Topology_O::throwIfExtractFragmentsAreNotExclusive(ConstitutionAtoms_sp residue)
    {
	uint numAtoms;
	map<string,OverlappingFragments>	overlaps;
	adapt::StringSet_sp missingAtomNames = residue->atomNamesAsStringSet();

	__checkExtractFragment(this->_ExtractCoreFragment,missingAtomNames,overlaps);
	for ( gctools::Vec0<ExtractFragment_sp>::iterator it=this->_ExtractFragments.begin();
	      it!=this->_ExtractFragments.end(); it++ )
	{
	    __checkExtractFragment((*it),missingAtomNames,overlaps);
	}
	if ( missingAtomNames->size() != 0 )
	{
	    stringstream ss;
	    ss << "The Topology " << this->getName() << " ExtractFragments do not include all atoms" << std::endl;
	    ss << " of the Constitution.  The following atom names are missing: " << missingAtomNames->asString();
	    SIMPLE_ERROR(BF("%s")%ss.str());
	}
	stringstream so;
	bool sawOverlaps = false;
	for( map<string,OverlappingFragments>::iterator oi = overlaps.begin(); oi!=overlaps.end(); oi++ )
	{
	    if ( oi->second._TimesSeen > 1 )
	    {
		sawOverlaps = true;
		so << "Multiple references for atom(" <<oi->first<<") in fragments: ";
		for ( set<string>::iterator si=oi->second._Fragments.begin(); si!=oi->second._Fragments.end();si++ )
		{
		    so << *si << " ";
		}
	    }
	}
	if ( sawOverlaps )
	{
	    stringstream se;
	    se << "In definition of Topology(" << this->getName() << ") there were overlapping ExtractFragment definitions"<<std::endl;
	    se << so.str();
	    SIMPLE_ERROR(BF(se.str()));
	}
    }
#endif

struct OverlappingFragments {
  uint	_TimesSeen;
  set<string>	_Fragments;
};

CL_LISPIFY_NAME(connect_residues);
CL_DEFUN void connect_residues(Topology_sp prev_topology,
                               Residue_sp prev_residue,
                               core::Symbol_sp out_plug_name,
                               Topology_sp next_topology,
                               Residue_sp next_residue,
                               core::Symbol_sp in_plug_name)
{
//  printf("%s:%d  prev_topology@%p next_topology@%p\n", __FILE__, __LINE__, prev_topology.raw_(), next_topology.raw_());
//  printf("%s:%d  prev_residue@%p next_residue@%p\n", __FILE__, __LINE__, prev_residue.raw_(), next_residue.raw_());
//  printf("%s:%d  out_plug_name = %s    in_plug_name = %s\n", __FILE__, __LINE__, _rep_(out_plug_name).c_str(), _rep_(in_plug_name).c_str() );
  Plug_sp out_plug = prev_topology->plugNamed(out_plug_name);
  Plug_sp in_plug = next_topology->plugNamed(in_plug_name);
  core::Symbol_sp out_plug_atom_name = out_plug->getB0();
  core::Symbol_sp in_plug_atom_name = in_plug->getB0();
  if ( out_plug_atom_name.nilp() ) {
    SIMPLE_ERROR(BF("Could not find b0 in plug %s for topology %s") % _rep_(out_plug) % _rep_(prev_topology));
  }
  if ( in_plug_atom_name.nilp() ) {
    SIMPLE_ERROR(BF("Could not find b0 in plug %s for topology %s") % _rep_(in_plug) % _rep_(next_topology));
  }
//  printf("%s:%d  out_plug_atom_name = %s  in_plug_atom_name = %s\n", __FILE__, __LINE__, _rep_(out_plug_atom_name).c_str(), _rep_(in_plug_atom_name).c_str());
  Atom_sp out_atom = prev_residue->atomWithName(out_plug_atom_name);
  Atom_sp in_atom = next_residue->atomWithName(in_plug_atom_name);
//  printf("%s:%d  out_atom = %s  in_atom = %s\n", __FILE__, __LINE__, _rep_(out_atom).c_str(), _rep_(in_atom).c_str());
  BondOrder bo = in_plug->getBondOrder0();
  in_atom->bondTo(out_atom, bo);
  if (in_plug->getB1().notnilp()) {
    Atom_sp out_atom = prev_residue->atomWithName(out_plug->getB1());
    Atom_sp in_atom = next_residue->atomWithName(in_plug->getB1());
    BondOrder bo = in_plug->getBondOrder1();
    in_atom->bondTo(out_atom, bo);
  }
}


CL_LISPIFY_NAME("getMonomerContext");
CL_DEFMETHOD     MonomerContext_sp Topology_O::getMonomerContext(CandoDatabase_sp bdb)
    {
	MonomerContext_sp context = MonomerContext_O::create();
	EntityNameSet_sp selfSet = EntityNameSet_O::create();
	Constitution_sp constitution = this->getConstitution();
        (void)constitution;
	context->setFocus(selfSet);
	IMPLEMENT_ME();
    }


CL_LISPIFY_NAME("properties");
CL_DEFMETHOD     core::HashTableEq_sp Topology_O::properties() const
    {_OF();
	return this->_Properties;
    }


    bool Topology_O::hasInPlug()
    {
      for ( auto i : this->_Plugs ) if ( i.second->getIsIn() ) return true;
      return false;
    }




CL_LISPIFY_NAME("getInPlug");
CL_DEFMETHOD     Topology_O::plugType Topology_O::getInPlug()
    {
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    if ( i->second->getIsIn() ) 
	    {
		return i->second;
	    }
	}
	return _Nil<plugOType>();
    }


CL_LISPIFY_NAME("plugsAsList");
CL_DEFMETHOD     core::List_sp Topology_O::plugsAsList()
    {_OF();
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
	core::Cons_sp cur = first;
	LOG(BF("The number of plugs = %d") % this->_Plugs.size()  );
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    LOG(BF("Adding plug: %s") % i->second->getName());
	    core::Cons_sp one = core::Cons_O::create(i->second,_Nil<core::T_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	LOG(BF("Returning the plugs"));
	return first->cdr();
    }


CL_LISPIFY_NAME("plugsWithMatesAsList");
CL_DEFMETHOD     core::List_sp Topology_O::plugsWithMatesAsList()
    {
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
	core::Cons_sp cur = first;
	LOG(BF("The number of plugs = %d") % this->_Plugs.size()  );
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    // skip origin plugs
	    if ( !i->second.isA<PlugWithMates_O>() ) continue;
	    LOG(BF("Adding plug: %s") % _rep_(i->second) );
	    core::Cons_sp one = core::Cons_O::create(i->second,_Nil<core::T_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	LOG(BF("Returning the plugs") );
	return first->cdr();
    }

CL_LISPIFY_NAME("outPlugsAsList");
CL_DEFMETHOD     core::List_sp Topology_O::outPlugsAsList()
    {
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
	core::Cons_sp cur = first;
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    Plug_sp plug = i->second;
	    if ( !plug->getIsIn() )
	    {
		core::Cons_sp one = core::Cons_O::create(i->second,_Nil<core::T_O>());
		cur->setCdr(one);
		cur = one;
	    }
	}
	return first->cdr();
    }

    string	Topology_O::description() const
    {
	stringstream	ss;
	ss << "#<TOPOLOGY";
	Topology_O* me = const_cast<Topology_O*>(this);
	ss << " Name(" << me->getName() << ")";
	ss << " Constitution("<< me->getConstitution()->getName() << ")";
	ss << " Plugs:";
	for ( Plugs::const_iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
          ss << (BF("%s@%p ") % _rep_(i->second->getName()) % i->second->getName().get() ).str();
	}
	ss << ">";
	return ss.str();
    }


    bool Topology_O::hasMatchingPlugsWithMates(adapt::SymbolSet_sp plugSet)
    {
	adapt::SymbolSet_sp myPlugSet = adapt::SymbolSet_O::create();
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
          LOG(BF("Looking at plug[%s]") % _rep_(i->second->getName()) );
	    if (i->second.isA<PlugWithMates_O>() )
	    {
		myPlugSet->insert(i->first);
	    }
	}
	LOG(BF("Topology plugs are: %s") % myPlugSet->asString() );
	bool match = plugSet->equal(myPlugSet);
	LOG(BF("Do they match the plugs passed as arguments[%s] --> %d") % plugSet->asString() % match );
	return match;
    }



    bool	Topology_O::matchesMonomerEnvironment( Monomer_sp mon )
    {
	LOG(BF("Checking if monomer[%s] matches the topology environment") % mon->description() );
	uint numPlugsWithMates = 0;
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
          LOG(BF("Looking at plug[%s]") % _rep_(i->second->getName()) );
	    if (!i->second.isA<PlugWithMates_O>() )
	    {
		LOG(BF("It's not a PlugWithMates"));
		continue;
	    }
	    if ( !mon->hasCouplingWithPlugName(i->second->getName()) ) 
	    {
              LOG(BF("The monomer doesn't have a coupling with plug name[%s]") % _rep_(i->second->getName()) );
		return false;
	    }
	    numPlugsWithMates++;
	}
	if ( numPlugsWithMates != mon->numberOfCouplings() ) 
	{
	    LOG(BF("There is a mismatch with the number of plugs in the topology[%d] and the number of couplings for the monomer[%d]") % numPlugsWithMates % mon->numberOfCouplings() );
	    return false;
	}
	LOG(BF("They match"));
	return true;
    }



    RingClosingPlug_sp Topology_O::provideMissingRingClosingPlug( Monomer_sp mon )
    {
	RingClosingPlug_sp missingRingClosingPlug = _Nil<RingClosingPlug_O>();
	uint numPlugsWithMates = 0;
	for ( Plugs::iterator i=this->_Plugs.begin();
	      i!= this->_Plugs.end(); i++)
	{
	    if (!i->second.isA<PlugWithMates_O>() ) continue;
	    if ( i->second.isA<RingClosingPlug_O>() )
	    {
		missingRingClosingPlug = i->second.as<RingClosingPlug_O>();
		continue;
	    }
	    if ( !mon->hasCouplingWithPlugName(i->second->getName()) ) return _Nil<RingClosingPlug_O>();
	    numPlugsWithMates++;
	}
	if ( numPlugsWithMates != mon->numberOfCouplings() ) return _Nil<RingClosingPlug_O>();
	return missingRingClosingPlug;
    }




CL_LISPIFY_NAME("matchesContext");
CL_DEFMETHOD     bool	Topology_O::matchesContext(MonomerContext_sp cm)
    {
	uint numPlugsWithMates = 0;
	for ( Plugs::iterator i=this->_Plugs.begin();
	      i!= this->_Plugs.end(); i++)
	{
	    if (!i->second.isA<PlugWithMates_O>() ) continue;
	    if ( !cm->hasNeighborWithCouplingName(i->second->getName()) ) return false;
	    numPlugsWithMates++;
	}
	if ( numPlugsWithMates != cm->numberOfNeighbors() ) return false;
	return true;
    }

CL_LISPIFY_NAME("hasPlugNamed");
CL_DEFMETHOD     bool	Topology_O::hasPlugNamed(core::Symbol_sp name)
    {
	bool res = this->_Plugs.contains(name);
	LOG(BF("Result = %d") % res );
	return res;
    };

CL_LISPIFY_NAME("plugNamed");
CL_DEFMETHOD     Plug_sp Topology_O::plugNamed(core::Symbol_sp name)
    {
      bool res = this->_Plugs.contains(name);
      if ( !res ) {
        SIMPLE_ERROR(BF("Could not find plug with name %s") % _rep_(name));
      }
      return this->_Plugs.get(name);
    };





    void	Topology_O::initialize()
    {
	this->Base::initialize();
	this->_Name = _Nil<core::Symbol_O>();
	this->_Plugs.clear();
	this->_SuppressTrainers = false;
	this->_Properties = core::HashTableEq_O::create_default();
    }





CL_LISPIFY_NAME("lookupOrCreateStereoisomerAtoms");
CL_DEFMETHOD     StereoisomerAtoms_sp Topology_O::lookupOrCreateStereoisomerAtoms(core::Symbol_sp stereoisomerName)
    {_OF();
	adapt::SymbolMap<StereoisomerAtoms_O>::iterator it= this->_StereoisomerAtomProperties.find(stereoisomerName);
	StereoisomerAtoms_sp result;
	if ( it==this->_StereoisomerAtomProperties.end() )
	{
	    Constitution_sp constitution = this->getConstitution();
	    ASSERTF(constitution->hasStereoisomerWithName(stereoisomerName),BF("Could not find stereoisomer named[%s] in constitution[%s]") % _rep_(stereoisomerName) % constitution->getName());
	    ConstitutionAtoms_sp constitutionAtoms = constitution->getConstitutionAtoms();
	    result = StereoisomerAtoms_O::create(constitutionAtoms);
	    this->_StereoisomerAtomProperties.set(stereoisomerName,result);
	} else
	{
	    result = this->_StereoisomerAtomProperties.get(stereoisomerName);
	}
	return result;
    }

//
// Destructor
//








};
