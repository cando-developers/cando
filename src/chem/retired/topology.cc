#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/common.h>
#include <clasp/core/stringSet.h>
#include <clasp/core/symbolSet.h>
#include <cando/chem/topology.h>
#include <cando/chem/coupling.h>
#include <cando/chem/monomer.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <clasp/core/environment.h>
#include <clasp/core/binder.h>
//#include "kinematics/atomTemplate.h"
//#include "kinematics/chiDihedrals.h"
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <clasp/core/wrappers.h>



namespace chem
{

#ifdef XML_ARCHIVE
    void	Topology_O::archiveBase(core::ArchiveP node)
    {_G();
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
	LOG(BF("About to get core for topology (%s)") % this->_Name->__repr__()  );
//	node->attributeIfNotNil("atomTreeTemplate",this->_AtomTreeTemplate);
//	node->attributeIfNotNil("chiList",this->_ChiList);
	node->attribute("properties",this->_Properties);
    }
#endif

#if INIT_TO_FACTORIES

#define ARGS_Topology_O_make "(name net_charge properties plugs residue)"
#define DECL_Topology_O_make ""
#define DOCS_Topology_O_make "make Topology"
    Topology_sp Topology_O::make(core::Symbol_sp name, int netCharge, core::HashTableEq_sp properties, core::Cons_sp curPlugs, ConstitutionAtoms_sp residue ) //, kinematics::AtomTemplate_sp atomTreeTemplate, kinematics::ChiList_sp chiList)
  {_G();
      GC_ALLOCATE(Topology_O, me );
    me->_Name = name;
    me->_ResidueNetCharge = netCharge;
    //
    // If "plugs" option is outgoing plugs only then we should
    // change the name to "outPlugs"
    //
    me->_Plugs.clear();
    core::Cons_sp curPlug;
    for ( ; curPlug.notnilp(); curPlug = curPlug->cdr() )
      {
	plugType p = curPlug->car<plugOType>();
	me->addPlug(p->getName(),p);
      }
#if 0
    // Here handle the AtomTreeBuilder
    me->_AtomTreeTemplate = atomTreeTemplate;
    
    // Define the order of rotamer dihedrals
    me->_ChiList = chiList;
    ASSERTF(me->_AtomTreeTemplate.notnilp(),BF("You must provide the atomTreeTemplate"));
#endif
    properties->mapHash( [&me] (core::T_sp key, core::T_sp val) {
            me->_Properties->setf_gethash(key.as<core::Symbol_O>(),val);
        } );
    return me;
  };
  
#else

   core::T_sp Topology_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
	this->_Name = translate::from_object<core::Symbol_O>::convert(env->lookup(ChemPkg,"name"));
	this->_ResidueNetCharge = translate::from_object<int>::convert(env->lookup(ChemPkg,"netCharge"));
	core::HashTableEq_sp properties = env->lookup(ChemPkg,"properties").as<core::HashTableEq_O>();
	//
	// If "plugs" option is outgoing plugs only then we should
	// change the name to "outPlugs"
	//
	core::Cons_sp curPlug = translate::from_object<core::Cons_O>::convert(env->lookup(ChemPkg,"plugs"));
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
    void __checkExtractFragment(ExtractFragmentBase_sp extractFrag,core::StringSet_sp missingAtomNames, map<string,OverlappingFragments>& overlaps)
    {_G();

	Fragment_sp frag = extractFrag->getFragment();
	core::StringSet_sp atomNames = frag->getAtomNames();
	for ( core::StringSet_O::iterator ni = atomNames->begin(); ni!=atomNames->end(); ni++ )
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
    {_G();
	uint numAtoms;
	map<string,OverlappingFragments>	overlaps;
	core::StringSet_sp missingAtomNames = residue->atomNamesAsStringSet();

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


    MonomerContext_sp Topology_O::getMonomerContext(CandoDatabase_sp bdb)
    {_G();
	MonomerContext_sp context = MonomerContext_O::create();
	EntityNameSet_sp selfSet = EntityNameSet_O::create();
	Constitution_sp constitution = this->getConstitution();
	context->setFocus(selfSet);
	IMPLEMENT_ME();
    }


    core::HashTableEq_sp Topology_O::properties() const
    {_OF();
	return this->_Properties;
    }


    bool Topology_O::hasInPlug()
    {_G();
	DEPRECIATED();
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    if ( i->second->getIsIn() ) return true;
	}
	return false;
    }




    Topology_O::plugType Topology_O::getInPlug()
    {_G();
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    if ( i->second->getIsIn() ) 
	    {
		return i->second;
	    }
	}
	return _Nil<plugOType>();
    }


    core::Cons_sp Topology_O::plugsAsCons()
    {_OF();
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
	core::Cons_sp cur = first;
	LOG(BF("The number of plugs = %d") % this->_Plugs.size()  );
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    LOG(BF("Adding plug: %s") % i->second->getName());
	    core::Cons_sp one = core::Cons_O::create(i->second,_Nil<core::Cons_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	LOG(BF("Returning the plugs"));
	return first->cdr();
    }


    core::Cons_sp Topology_O::plugsWithMatesAsCons()
    {_G();
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
	core::Cons_sp cur = first;
	LOG(BF("The number of plugs = %d") % this->_Plugs.size()  );
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    // skip origin plugs
	    if ( !i->second->isAssignableTo<PlugWithMates_O>() ) continue;
	    LOG(BF("Adding plug: %s") % i->second->__repr__() );
	    core::Cons_sp one = core::Cons_O::create(i->second,_Nil<core::Cons_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	LOG(BF("Returning the plugs") );
	return first->cdr();
    }

    core::Cons_sp Topology_O::outPlugsAsCons()
    {_G();
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
	core::Cons_sp cur = first;
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    Plug_sp plug = i->second;
	    if ( !plug->getIsIn() )
	    {
		core::Cons_sp one = core::Cons_O::create(i->second,_Nil<core::Cons_O>());
		cur->setCdr(one);
		cur = one;
	    }
	}
	return first->cdr();
    }



    Constitution_sp	Topology_O::getConstitution()
    {_OF();
	IMPLEMENT_MEF(BF("Handle new way of dealing with owners"));
#if 0
	ASSERTNOTNULL(this->_WeakConstitution);
	if ( this->_WeakConstitution.lock().nilp() )
	{
	    // Crawl up the owners chain and look
	    // for a Constitution
	    this->_WeakConstitution = this->ownerWithClass<Constitution_O>();
	}
	return this->_WeakConstitution.lock();
#endif
    } 

    string	Topology_O::description() const
    {
	stringstream	ss;
	ss << this->core::T_O::description();
	Topology_O* me = const_cast<Topology_O*>(this);
	ss << " Name(" << me->getName() << ")";
	ss << " Constitution("<< me->getConstitution()->getName() << ")";
	ss << " Plugs:";
	for ( Plugs::const_iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    ss << (BF("%s@%p ") % i->second->getName()->__repr__() % i->second->getName().get() ).str();
	}
	ss << "]";
	return ss.str();
    }


    bool Topology_O::hasMatchingPlugsWithMates(core::SymbolSet_sp plugSet)
    {_G();
	core::SymbolSet_sp myPlugSet = core::SymbolSet_O::create();
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    LOG(BF("Looking at plug[%s]") % i->second->getName()->__repr__() );
	    if (i->second->isAssignableTo<PlugWithMates_O>() )
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
    {_G();
	LOG(BF("Checking if monomer[%s] matches the topology environment") % mon->description() );
	uint numPlugsWithMates = 0;
	for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
	{
	    LOG(BF("Looking at plug[%s]") % i->second->getName()->__repr__() );
	    if (!i->second->isAssignableTo<PlugWithMates_O>() )
	    {
		LOG(BF("It's not a PlugWithMates"));
		continue;
	    }
	    if ( !mon->hasCouplingWithPlugName(i->second->getName()) ) 
	    {
		LOG(BF("The monomer doesn't have a coupling with plug name[%s]") % i->second->getName()->__repr__() );
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
    {_G();
	RingClosingPlug_sp missingRingClosingPlug = _Nil<RingClosingPlug_O>();
	uint numPlugsWithMates = 0;
	for ( Plugs::iterator i=this->_Plugs.begin();
	      i!= this->_Plugs.end(); i++)
	{
	    if (!i->second->isAssignableTo<PlugWithMates_O>() ) continue;
	    if ( i->second->isAssignableTo<RingClosingPlug_O>() )
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



    bool Topology_O::hasFlag(core::Symbol_sp f) const
    {_OF();
	return this->_Flags.contains(f);
    }


    bool	Topology_O::matchesContext(MonomerContext_sp cm)
    {
	uint numPlugsWithMates = 0;
	for ( Plugs::iterator i=this->_Plugs.begin();
	      i!= this->_Plugs.end(); i++)
	{
	    if (!i->second->isAssignableTo<PlugWithMates_O>() ) continue;
	    if ( !cm->hasNeighborWithCouplingName(i->second->getName()) ) return false;
	    numPlugsWithMates++;
	}
	if ( numPlugsWithMates != cm->numberOfNeighbors() ) return false;
	return true;
    }

    bool	Topology_O::hasPlugNamed(core::Symbol_sp name)
    {_OF();
#ifdef DEBUG_ON
	LOG(BF("Looking for plug name[%s@%p]") % name->__repr__() % name.get() );
	for (Plugs::const_iterator ki=this->_Plugs.begin(); ki!=this->_Plugs.end(); ki++ )
	{
	    LOG(BF(" available Symbol key: %s@%p") % ki->first->__repr__() % ki->first.get() );
	}
#endif
	bool res = this->_Plugs.contains(name);
	LOG(BF("Result = %d") % res );
	return res;
    };

    Plug_sp Topology_O::plugNamed(core::Symbol_sp name)
    {_OF();
	LOG(BF("Looking for plug name[%s] available keys[%s]")
	    % name->__repr__()
	    % core::StringSet_O::create_fromKeysOfSymbolMap(this->_Plugs,_lisp)->asString());
	bool res = this->_Plugs.contains(name);
	LOG(BF("Result = %d") % res );
	if ( !res ) return _Nil<Plug_O>();
	return this->_Plugs.get(name);
    };





    void	Topology_O::initialize()
    {
	this->Base::initialize();
	this->_WeakConstitution = _Nil<Constitution_O>();
	this->_Name = _Nil<core::Symbol_O>();
	this->_ResidueNetCharge = 0;
	this->_Plugs.clear();
	this->_TemporaryObject = _Nil<core::T_O>();
	this->_SuppressTrainers = false;
	this->_Properties = core::HashTableEq_O::create_default();
    }



    void Topology_O::setTemporaryObject(core::T_sp o)
    {_G();
	ASSERTNOTNULL(o);
	this->_TemporaryObject = o;
    }

    core::T_sp Topology_O::getTemporaryObject()
    {_G();
	ASSERTNOTNULL(this->_TemporaryObject);
	return this->_TemporaryObject;
    }


    StereoisomerAtoms_sp Topology_O::lookupOrCreateStereoisomerAtoms(core::Symbol_sp stereoisomerName)
    {_OF();
	core::SymbolMap<StereoisomerAtoms_O>::iterator it= this->_StereoisomerAtomProperties.find(stereoisomerName);
	StereoisomerAtoms_sp result;
	if ( it==this->_StereoisomerAtomProperties.end() )
	{
	    Constitution_sp constitution = this->getConstitution();
	    ASSERTF(constitution->hasStereoisomerWithName(stereoisomerName),BF("Could not find stereoisomer named[%s] in constitution[%s]") % stereoisomerName->__repr__() % constitution->getName());
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




    void Topology_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Topology_O>()
////	    .def_raw("core:__init__",&Topology_O::__init__,"(self &key name netCharge plugs constitutionAtoms properties suppressTrainers)")
	    .def("suppressTrainers",&Topology_O::suppressTrainers)
	    .def("topology-suppress-trainers",&Topology_O::suppressTrainers)
	    .def("getResidueNetCharge",&Topology_O::getResidueNetCharge)
	    .def("getMonomerContext",&Topology_O::getMonomerContext)
	    .def("getConstitution",&Topology_O::getConstitution)
	    .def("numberOfPlugs",&Topology_O::numberOfPlugs)
	    .def("properties",&Topology_O::properties)
	    .def("plugsAsCons",&Topology_O::plugsAsCons)
	    .def("plugsWithMatesAsCons",&Topology_O::plugsWithMatesAsCons)
	    .def("outPlugsAsCons",&Topology_O::outPlugsAsCons)
	    .def("getInPlug",&Topology_O::getInPlug)
	    .def("getName",&Topology_O::getName,"","","",false)
	    .def("setTemporaryObject",&Topology_O::setTemporaryObject)
	    .def("getTemporaryObject",&Topology_O::getTemporaryObject)
	    .def("hasPlugNamed",&Topology_O::hasPlugNamed)
	    .def("plugNamed",&Topology_O::plugNamed)
	    .def("hasFlag", &Topology_O::hasFlag)
	    .def("matchesContext", &Topology_O::matchesContext)
	    .def("lookupOrCreateStereoisomerAtoms",&Topology_O::lookupOrCreateStereoisomerAtoms)
	    ;
    }

    void Topology_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Topology,"","",_lisp)
	    .def("suppressTrainers",&Topology_O::suppressTrainers)
	    .def("topology-suppress-trainers",&Topology_O::suppressTrainers)
	    .def("getResidueNetCharge",&Topology_O::getResidueNetCharge)
	    .def("getMonomerContext",&Topology_O::getMonomerContext)
	    .def("getConstitution",&Topology_O::getConstitution)
	    .def("numberOfPlugs",&Topology_O::numberOfPlugs)
	    .def("plugsAsCons",&Topology_O::plugsAsCons)
	    .def("plugsWithMatesAsCons",&Topology_O::plugsWithMatesAsCons)
	    .def("outPlugsAsCons",&Topology_O::outPlugsAsCons)
	    .def("getInPlug",&Topology_O::getInPlug)
//	    .def("getName",&Topology_O::getName,"","","",false)
	    .def("setTemporaryObject",&Topology_O::setTemporaryObject)
	    .def("getTemporaryObject",&Topology_O::getTemporaryObject)
	    .def("hasPlugNamed",&Topology_O::hasPlugNamed)
	    .def("plugNamed",&Topology_O::plugNamed)
	    .def("hasFlag", &Topology_O::hasFlag)
	    .def("matchesContext", &Topology_O::matchesContext)
	    .def("lookupOrCreateStereoisomerAtoms",&Topology_O::lookupOrCreateStereoisomerAtoms)
	    ;
#endif
    }



    EXPOSE_CLASS(chem,Topology_O);
};


