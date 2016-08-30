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
#define	DEBUG_LEVEL_FULL


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
#include <clasp/core/binder.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <clasp/core/wrappers.h>

namespace chem
{

#ifdef XML_ARCHIVE
    void	Topology_O::archiveBase(core::ArchiveP node)
    {
	if ( node->saving() ) this->getConstitution();
	node->archiveWeakPointer("constitution",this->_WeakConstitution);
#if PRODUCTION_CODE   // FIXME use "name" only and remove the test for "_key"
	node->attribute("name",this->_Name);
#else
	if ( node->loading() )
	{
	    if ( node->hasAttribute("name") )
	    {
		node->attribute("name",this->_Name);
	    } else
	    {
		node->attribute("_key",this->_Name);
	    }
	} else
	{
	    node->attribute("name",this->_Name);
	}
#endif
	node->attributeIfNotDefault("ResidueNetCharge",this->_ResidueNetCharge,0);
	node->attributeIfNotDefault("suppressTrainers",this->_SuppressTrainers,false);
	node->archiveSymbolMap( "topologyPlugs", this->_Plugs );
	node->archiveSymbolMapIfDefined("stereoisomerAtomProperties",this->_StereoisomerAtomProperties);
	LOG(BF("About to get core for topology (%s)") % _rep_(this->_Name)  );
	node->attribute("properties",this->_Properties);
    }
#endif

Topology_sp Topology_O::makeTopologyFromResidue(chem::Residue_sp residue, core::T_sp constitution )
{
  IMPLEMENT_MEF(BF("Implement makeTopologyFromResidue"));
}

Topology_sp make_topology(core::Symbol_sp name, Constitution_sp constitution, int netCharge, core::List_sp plugs)
{
  Topology_sp me = gctools::GC<Topology_O>::allocate(name,constitution,netCharge);
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


#if 0

   core::T_sp Topology_O::__init__(core::Function_sp exec, core::List_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->_Name = translate::from_object<core::Symbol_O>::convert(env->lookup(ChemPkg,"name"));
	this->_ResidueNetCharge = translate::from_object<int>::convert(env->lookup(ChemPkg,"netCharge"));
	core::HashTableEq_sp properties = env->lookup(ChemPkg,"properties").as<core::HashTableEq_O>();
	//
	// If "plugs" option is outgoing plugs only then we should
	// change the name to "outPlugs"
	//
	core::List_sp curPlug = translate::from_object<core::List_O>::convert(env->lookup(ChemPkg,"plugs"));
	this->_Plugs.clear();
	for ( ; curPlug.notnilp(); curPlug = curPlug->cdr() )
	{
	    plugType p = curPlug->car<plugOType>();
	    this->addPlug(p->getName(),p);
	}
	ConstitutionAtoms_sp residue = translate::from_object<ConstitutionAtoms_O>::convert(env->lookup(ChemPkg,"constitutionAtoms"));
#if 0
	// Here handle the AtomTreeBuilder
	this->_AtomTreeTemplate = translate::from_object<kinematics::AtomTemplate_O>::convert(env->lookup(Topology_O::Package(),"atomTreeTemplate"));

	// Define the order of rotamer dihedrals
	this->_ChiList = env->lookup(_lisp->internWithPackageName(ChemPkg,"chiList")).as<kinematics::ChiList_O>();
	ASSERTF(this->_AtomTreeTemplate.notnilp(),BF("You must provide the atomTreeTemplate"));
#endif
	this->_Properties->addAllBindings(properties);
	return _Nil<core::T_O>();
    }

#endif

    struct OverlappingFragments {
	uint	_TimesSeen;
	set<string>	_Fragments;
    };

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
	DEPRECIATED();
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    if ( i->second->getIsIn() ) return true;
	}
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


CL_LISPIFY_NAME("plugsAsCons");
CL_DEFMETHOD     core::List_sp Topology_O::plugsAsCons()
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


CL_LISPIFY_NAME("plugsWithMatesAsCons");
CL_DEFMETHOD     core::List_sp Topology_O::plugsWithMatesAsCons()
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

CL_LISPIFY_NAME("outPlugsAsCons");
CL_DEFMETHOD     core::List_sp Topology_O::outPlugsAsCons()
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
	ss << this->Base::description();
	Topology_O* me = const_cast<Topology_O*>(this);
	ss << " Name(" << me->getName() << ")";
	ss << " Constitution("<< me->getConstitution()->getName() << ")";
	ss << " Plugs:";
	for ( Plugs::const_iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
          ss << (BF("%s@%p ") % _rep_(i->second->getName()) % i->second->getName().get() ).str();
	}
	ss << "]";
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



CL_LISPIFY_NAME("hasFlag");
CL_DEFMETHOD     bool Topology_O::hasFlag(core::Symbol_sp f) const
    {_OF();
	return this->_Flags.contains(f);
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
    {_OF();
#ifdef DEBUG_ON
      LOG(BF("Looking for plug name[%s@%p]") % _rep_(name) % name.get() );
	for (Plugs::const_iterator ki=this->_Plugs.begin(); ki!=this->_Plugs.end(); ki++ )
	{
          LOG(BF(" available Symbol key: %s@%p") % _rep_(ki->first) % ki->first.get() );
	}
#endif
	bool res = this->_Plugs.contains(name);
	LOG(BF("Result = %d") % res );
	return res;
    };

CL_LISPIFY_NAME("plugNamed");
CL_DEFMETHOD     Plug_sp Topology_O::plugNamed(core::Symbol_sp name)
    {_OF();
	LOG(BF("Looking for plug name[%s] available keys[%s]")
	    % _rep_(name)
	    % adapt::StringSet_O::create_fromKeysOfSymbolMap(this->_Plugs)->asString());
	bool res = this->_Plugs.contains(name);
	LOG(BF("Result = %d") % res );
	if ( !res ) return _Nil<Plug_O>();
	return this->_Plugs.get(name);
    };





    void	Topology_O::initialize()
    {
	this->Base::initialize();
	this->_Name = _Nil<core::Symbol_O>();
	this->_ResidueNetCharge = 0;
	this->_Plugs.clear();
	this->_TemporaryObject = _Nil<core::T_O>();
	this->_SuppressTrainers = false;
	this->_Properties = core::HashTableEq_O::create_default();
    }



CL_LISPIFY_NAME("setTemporaryObject");
CL_DEFMETHOD     void Topology_O::setTemporaryObject(core::T_sp o)
    {
	ASSERTNOTNULL(o);
	this->_TemporaryObject = o;
    }

CL_LISPIFY_NAME("getTemporaryObject");
CL_DEFMETHOD     core::T_sp Topology_O::getTemporaryObject()
    {
	ASSERTNOTNULL(this->_TemporaryObject);
	return this->_TemporaryObject;
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
