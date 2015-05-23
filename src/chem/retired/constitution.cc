#define	DEBUG_LEVEL_FULL
//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/common.h>
#include <clasp/core/symbolSet.h>
#include <clasp/core/symbol.h>
#include <clasp/core/stringList.h>
#include <clasp/core/environment.h>
#include <cando/chem/symbolTable.h>
#include <cando/chem/constitution.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/loop.h>
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/monomer.h>
#include <cando/chem/plug.h>
#include <cando/chem/complexRestraints.h>
#include <cando/chem/topology.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/constitutionAtoms.h>
#include <clasp/core/wrappers.h>



namespace chem {


#if INIT_TO_FACTORIES

#define ARGS_Constitution_O_make "(name comment meta_constitution constitution_atoms stereo_information plugs topologies)"
#define DECL_Constitution_O_make ""
#define DOCS_Constitution_O_make "make Constitution"
    Constitution_sp Constitution_O::make(core::Symbol_sp name, const string& comment, core::Symbol_sp metaConstitution, ConstitutionAtoms_sp constitutionAtoms, StereoInformation_sp stereoInformation, core::Cons_sp plugs, core::Cons_sp topologies)
  {_G();
      GC_ALLOCATE(Constitution_O, me );
      me->_Name = name;
      me->_Comment = comment;
      me->_MetaConstitution = metaConstitution;
      me->_ConstitutionAtoms = constitutionAtoms;
      IMPLEMENT_MEF(BF("Handle setOwner"));
//      me->_ConstitutionAtoms->setOwner(me);
      me->_StereoInformation = stereoInformation;
      IMPLEMENT_MEF(BF("Handle setOwner"));
//      me->_StereoInformation->setOwner(me);
      me->_StereoInformation->validate();
      {_BLOCK_TRACE("Adding plugs to Constitution");
	  IMPLEMENT_MEF(BF("Handle setOwnerOfAllEntries"));
//	  plugs->setOwnerOfAllEntries(me);
	  me->_PlugsByName.clear();
	  for ( ; plugs.notnilp(); plugs = plugs->cdr() )
	  {
	      Plug_sp p = plugs->car<Plug_O>();
	      _BLOCK_TRACEF(BF("Adding plug[%s]") % p->getName() );
	      if ( me->_PlugsByName.count(p->getName())>0 )
	      {
		  SIMPLE_ERROR(BF("There is already a plug named: %s") %p->getName()->__repr__() );
	      }
	      core::Symbol_sp key = p->getName();
	      ASSERTF(!key->isKeywordSymbol(),BF("Don't use keyword symbols for plug names"));
	      me->_PlugsByName.set(p->getName(),p);
	  }
      }
      {_BLOCK_TRACE("Adding topologies to Constitution");
	  IMPLEMENT_MEF(BF("Handle setOwnerOfAllEntries"));
//	  topologies->setOwnerOfAllEntries(me);
	  me->_Topologies.clear();
	  for ( ; topologies.notnilp(); topologies = topologies->cdr() )
	  {
	      Topology_sp t = topologies->car<Topology_O>();
	      _BLOCK_TRACEF(BF("Adding topology[%s]@%p") % t->getName() % (void*)(t.get()) );
	      if ( me->_Topologies.count(t->getName())>0 )
	      {
		  SIMPLE_ERROR(BF("There is already a topology named: "+t->getName()->__repr__() ));
	      }
	      me->_Topologies.set(t->getName(),t);
	  }
      }
      return me;
  };

#else

    core::T_sp Constitution_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{_G();
    this->_Name = translate::from_object<core::Symbol_O>::convert(env->lookup(ChemPkg,"name"));
    this->_Comment = translate::from_object<string>::convert(env->lookup(ChemPkg,"comment"));
    this->_MetaConstitution = translate::from_object<core::Symbol_O>::convert(env->lookup(ChemPkg,"metaConstitution"));
    this->_ConstitutionAtoms = translate::from_object<ConstitutionAtoms_O>::convert(env->lookup(ChemPkg,"constitutionAtoms"));
    this->_ConstitutionAtoms->setOwner(this->sharedThis<Constitution_O>());
    this->_StereoInformation = translate::from_object<StereoInformation_O>::convert(env->lookup(ChemPkg,"stereoInformation"));
    this->_StereoInformation->setOwner(this->sharedThis<Constitution_O>());
    this->_StereoInformation->validate();
    {_BLOCK_TRACE("Adding plugs to Constitution");
	core::Cons_sp c = translate::from_object<core::Cons_O>::convert(env->lookup(ChemPkg,"plugs"));
	c->setOwnerOfAllEntries(this->sharedThis<Constitution_O>());
	this->_PlugsByName.clear();
	for ( ; c.notnilp(); c = c->cdr() )
	{
	    Plug_sp p = c->car<Plug_O>();
	    _BLOCK_TRACEF(BF("Adding plug[%s]") % p->getName() );
	    if ( this->_PlugsByName.count(p->getName())>0 )
	    {
		SIMPLE_ERROR(BF("There is already a plug named: "+p->getName()->__repr__() ));
	    }
	    core::Symbol_sp key = p->getName();
	    ASSERTF(!key->isKeywordSymbol(),BF("Don't use keyword symbols for plug names"));
	    this->_PlugsByName.set(p->getName(),p);
	}
    }
    {_BLOCK_TRACE("Adding topologies to Constitution");
	core::Cons_sp c = translate::from_object<core::Cons_O>::convert(env->lookup(ChemPkg,"topologies"));
	c->setOwnerOfAllEntries(this->sharedThis<Constitution_O>());
	this->_Topologies.clear();
	for ( ; c.notnilp(); c = c->cdr() )
	{
	    Topology_sp t = c->car<Topology_O>();
	    _BLOCK_TRACEF(BF("Adding topology[%s]@%p") % t->getName() % (void*)(t.get()) );
	    if ( this->_Topologies.count(t->getName())>0 )
	    {
		SIMPLE_ERROR(BF("There is already a topology named: "+t->getName()->__repr__() ));
	    }
	    this->_Topologies.set(t->getName(),t);
	}
    }
    return _Nil<core::T_O>();
}

#endif

    void Constitution_O::addStereoisomersToCandoDatabase(CandoDatabase_sp db)
    {_OF();
	for ( stereoisomerIterator si=this->begin_Stereoisomers(); 
	      si!=this->end_Stereoisomers(); si++ )
	{
	    db->addEntity(*si);
	}
    }




/*!
	Return a copy of the residue that this constitution defines
*/
void	Constitution_O::makeResidueConsistentWithStereoisomerNamed(Residue_sp res,
								   core::Symbol_sp stereoisomerName)
{_G();
    core::Symbol_sp	fullName, pdbName;
    Stereoisomer_sp				si;
    string					atomName;
    Constitution_sp				residueConstitution, stereoisomerConstitution;
    Atom_sp					aa;
    CandoDatabase_sp			bdb;
    residueConstitution = res->getConstitution();
    stereoisomerConstitution = getCandoDatabase()->constitutionForNameOrPdb(stereoisomerName);

    if (!residueConstitution->isSameAs(stereoisomerConstitution)) {
	LOG(BF("residueConstitution = %s")%residueConstitution->__repr__());
	SIMPLE_ERROR(BF("The residue constitution(%s) is different from the "
				   "stereoisomer(%s) constitution(%s) that you want to change it too")
		     % residueConstitution->getName() % stereoisomerName % stereoisomerConstitution->getName() );
    }
    bdb = getCandoDatabase();
    fullName = bdb->getMonomerNameForNameOrPdb(stereoisomerName);
    pdbName = bdb->getPdbNameForNameOrPdb(stereoisomerName);
    res->setName(fullName);
    res->setPdbName(pdbName);
    //
    //
    // Set chiral restraints
    //
    //
    si = this->_StereoInformation->getStereoisomer(stereoisomerName);
    gctools::Vec0<StereoConfiguration_sp>::iterator	sci;
    for (sci=si->_Configurations_begin();sci!=si->_Configurations_end();sci++){
	aa = res->atomWithName((*sci)->getAtomName());
	LOG(BF("Setting the configuration of atom(%s) to(%s)") % aa->description().c_str() % (*sci)->getConfiguration().c_str()  ); //
	if ( (*sci)->getConfiguration() == chemkw::_sym_S ) {
	    aa->setConfiguration( S_Configuration );
	} else if ( (*sci)->getConfiguration() == chemkw::_sym_R ) {
	    aa->setConfiguration( R_Configuration );
	}
    }
    //
    //
    // Set dihedral restraints for E/Z pi bonds
    //
    //
    gctools::Vec0<ComplexRestraint_sp>::iterator	tpi;
    for ( tpi=this->_StereoInformation->begin_ComplexRestraints(); 
	  tpi!=this->_StereoInformation->end_ComplexRestraints(); tpi++ )
    {
	(*tpi)->fillRestraints(res);
    }
}

    bool	Constitution_O::recognizesStereoisomerName(core::Symbol_sp nm)
{ 
    return this->_StereoInformation->recognizesNameOrPdb(nm);
};


RepresentativeList_sp	Constitution_O::expandedRepresentativeList() const
{_G();
    gctools::Vec0<Stereoisomer_sp>::const_iterator	si;
    RepresentativeList_sp allRepresentatives  = RepresentativeList_O::create();
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++)
    {
	RepresentativeList_sp oneList = (*si)->expandedRepresentativeList();
	allRepresentatives->vectorPushExtend(oneList);
    }
    return allRepresentatives;
}



    RingClosingPlug_sp Constitution_O::getMissingRingClosingPlug(Monomer_sp mon, Monomer_sp mate)
{_G();
    RingClosingPlug_sp plug, missing;
    missing = _Nil<RingClosingPlug_O>();
    TopologyMap::iterator	ti;
    gctools::SmallOrderedSet<Topology_sp>	candidateTopologies;
    for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ ) 
    {
	plug = (ti->second)->provideMissingRingClosingPlug(mon);
	if ( plug.notnilp() ) 
	{
	    if ( plug->recognizesRingClosingMate(mate->getName()) )
	    {
		candidateTopologies.insert(ti->second);
		missing = plug;
	    }
	}
    }
    if ( candidateTopologies.size() > 1 )
    {
	stringstream ss;
	ss << "In Constitution(" << this->constitutionName() << " MetaConstitution(" << this->getMetaConstitutionName() << ")" << std::endl;
	ss << "there are more than one topologies that match the current monomer environment: " << mon->description() << std::endl;
	ss << "that are missing a ring closing plug"<<std::endl;
	SIMPLE_ERROR(BF("%s") % ss.str() );
    }
    if ( candidateTopologies.size() == 0 )
    {
	SIMPLE_ERROR(BF("There are no Topologies with missing ring closing plugs"));
    }
    return missing;
}

    core::Symbol_sp Constitution_O::pdbFromNameOrPdb(core::Symbol_sp nm)
    {
	return this->_StereoInformation->pdbFromNameOrPdb(nm);
    };


core::StringList_sp Constitution_O::getMonomerNamesAsStringList()
{
	return this->_StereoInformation->getMonomerNamesAsStringList();
};
core::SymbolSet_sp	Constitution_O::getMonomerNamesAsSymbolSet()
{
    return this->_StereoInformation->getMonomerNamesAsSymbolSet();
};

core::StringList_sp Constitution_O::getPdbNamesAsStringList() {
	return this->_StereoInformation->getPdbNamesAsStringList();
};




    bool Constitution_O::hasStereoisomerWithName(core::Symbol_sp stereoisomerName)
{_G();
    gctools::Vec0<Stereoisomer_sp>::iterator	si;
    LOG(BF("Looking for Stereoisomer with name(%s)") % stereoisomerName->__repr__() );
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++)
    {
	LOG(BF("Looking at stereoisomer name(%s)") % (*si)->getName()->__repr__() );
	if ( (*si)->getName() == stereoisomerName )
	{
	    LOG(BF("Found it!!!!") ); //
	    return true;
	}
    }
    LOG(BF("Could not find it") ); //
    return false;
}

    Stereoisomer_sp Constitution_O::getStereoisomerWithName(core::Symbol_sp stereoisomerName) const
{_G();
    const_stereoisomerIterator	si;
    LOG(BF("Looking for Stereoisomer with name(%s)") % stereoisomerName->__repr__() );
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++)
    {
	LOG(BF("Looking at stereoisomer name(%s)") % (*si)->getName()->__repr__() );
	if ( (*si)->getName() == stereoisomerName )
	{
	    LOG(BF("Found it!!!!") ); //
	    return *si;
	}
    }
    LOG(BF("Could not find it") ); //
    return _Nil<Stereoisomer_O>();
}



core::Cons_sp Constitution_O::stereoisomersAsCons() { return this->_StereoInformation->stereoisomersAsCons(); };
core::Cons_sp Constitution_O::topologiesAsCons() {
    core::Cons_sp result = _Nil<core::Cons_O>();
    for ( TopologyMap::iterator it = this->_Topologies.begin(); it!=this->_Topologies.end(); ++it ) {
        result = core::Cons_O::create(it->second,result);
    }
    return result;
}

core::Cons_sp Constitution_O::plugsAsCons() {
    core::Cons_sp result = _Nil<core::Cons_O>();
    for ( PlugMap::iterator it = this->_PlugsByName.begin(); it!=this->_PlugsByName.end(); ++it ) {
        result = core::Cons_O::create(it->second,result);
    }
    return result;
}
core::Cons_sp Constitution_O::plugsWithMatesAsCons()
{_G();
    core::Cons_sp first = _Nil<core::Cons_O>();
    PlugMap::iterator mi;
    for ( mi=this->_PlugsByName.begin(); mi!=this->_PlugsByName.end(); mi++ )
    {
	if ( mi->second->hasMates() )
	{
	    first = core::Cons_O::create(mi->second,first);
	}
    }
    return first;
}

Constitution_O::stereoisomerIterator Constitution_O::begin_Stereoisomers() 
{
    return this->_StereoInformation->begin_Stereoisomers();
};
Constitution_O::stereoisomerIterator Constitution_O::end_Stereoisomers() 
{
    return this->_StereoInformation->end_Stereoisomers();
};

Constitution_O::const_stereoisomerIterator Constitution_O::begin_Stereoisomers() const
{
    return this->_StereoInformation->begin_Stereoisomers();
};
Constitution_O::const_stereoisomerIterator Constitution_O::end_Stereoisomers() const
{
    return this->_StereoInformation->end_Stereoisomers();
};


/*!
	Return a copy of the residue that this constitution defines
*/
    Residue_sp	Constitution_O::createResidueForStereoisomerName(core::Symbol_sp stereoisomerName)
{_G();
    Residue_sp res = this->_ConstitutionAtoms->makeResidue();
    IMPLEMENT_MEF(BF("Handle owners"));
//    res->setOwner(this->sharedThis<Constitution_O>());
    this->makeResidueConsistentWithStereoisomerNamed(res,stereoisomerName);
    return res;
}


    core::Symbol_sp	Constitution_O::nameFromNameOrPdb(core::Symbol_sp nm)
{ _G();
    ASSERTNOTNULL(this->_StereoInformation);
    return this->_StereoInformation->nameFromNameOrPdb(nm);
};



string	Constitution_O::description() const
{
stringstream	ss;
    ss << this->className() << "[";
    ss << "Name(";
    if ( this->_Name.not_pointerp() ) ss << "--UNDEFINED--";
    else ss << this->_Name->__repr__();
    ss <<") ";
    ss << "]";
    return ss.str();
}

core::SymbolSet_sp Constitution_O::getPlugNames()
{_G();
    core::SymbolSet_sp ss = core::SymbolSet_O::create();
    PlugMap::iterator pi;
    for ( pi=this->_PlugsByName.begin(); pi!=this->_PlugsByName.end(); pi++ )
    {
	ss->insert(pi->first);
    }
    return ss;
}


    Topology_sp	Constitution_O::simplestTopologyWithPlugNamed(core::Symbol_sp name)
{_G();
    TopologyMap::iterator	ti;
Topology_sp			tres;
int				connects,temp;
tres = _Nil<Topology_O>();
    connects = 9999999;
    for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ )
    {
	LOG(BF("Looking for plug[%s]@%p in topology[%s]") % name->__repr__() % name.get() % ti->second->description());
	if ( ti->second->hasPlugNamed(name) )
	{
	    LOG(BF("Found the plug[%s]") % name->__repr__() );
	    temp = ti->second->numberOfPlugs();
	    if ( temp < connects )
	    {
		LOG(BF("Found one with only %d plugs") % temp );
		connects = temp;
		tres = ti->second;
	    } else
	    {
		LOG(BF("Did not find the plug[%s]") % name->__repr__() );
	    }
	}
    }
#ifdef	DEBUG_ON
    if ( tres.notnilp() )
    {
        LOG(BF("simplest topology is named: %s") % tres->getName()->__repr__()  );
    } else
    {
        LOG(BF("there is no topology that has the plug[%s]") % name->__repr__() );
    }
#endif
    return tres;
}


    Topology_sp	Constitution_O::topologyWithName(core::Symbol_sp name) const
    {_OF();
	Topology_sp			tres;
	ASSERTF(this->_Topologies.contains(name),BF("There is no topology with name[%s] in constitution[%s]") % name->__repr__() % this->getName()->__repr__() );
	return this->_Topologies.get(name);
    }






Topology_sp	Constitution_O::getTopologyForMonomerEnvironment(Monomer_sp mon )
{_G();
    TopologyMap::iterator	ti;
Topology_sp			tres;
tres = _Nil<Topology_O>();
    for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ ) 
    {
	LOG(BF("Checking topology(%s)")% ti->first->__repr__() );
	if ( ti->second->matchesMonomerEnvironment(mon) ) 
	{
	    return ti->second;
	}
    }
    LOG(BF("No topology matched the environment"));
    return tres;
}


Topology_sp	Constitution_O::getTopologyForContext(MonomerContext_sp cont )
{_G();
    TopologyMap::iterator	ti;
Topology_sp			tres;
tres = _Nil<Topology_O>();
    for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ ) {
	if ( ti->second->matchesContext(cont) ) {
	    return ti->second;
	}
    }
    return tres;
}





RepresentedEntityNameSet_sp	Constitution_O::asGroup()
{_G();
RepresentedEntityNameSet_sp				group;
gctools::Vec0<Stereoisomer_sp>::iterator	si;
    group = RepresentedEntityNameSet_O::create(); //,getCandoDatabase());

//    group->setName(this->getName());
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++){
	group->addMonomerName((*si)->getName());
    }
    return group;
}






#ifdef	COUPLING_MAP //[
//
//	defineCouplingIn
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingIn( string in )
{
core::VectorStrings	out;

    this->defineCouplingInOutVector(in,out);
}


//
//	defineCouplingInOut1
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingInOut1( string in, string out1 )
{
core::VectorStrings out;

    out.push_back(out1);
    this->defineCouplingInOutVector(in,out);
}



//
//	defineCouplingInOut2
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingInOut2( string in,
						string out1,
						string out2 )
{
core::VectorStrings out;

    out.push_back(out1);
    out.push_back(out2);
    this->defineCouplingInOutVector(in,out);
}



//
//	defineCouplingInOut3
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingInOut3( string in,
						string out1,
						string out2,
						string out3 )
{
core::VectorStrings out;

    out.push_back(out1);
    out.push_back(out2);
    out.push_back(out3);
    this->defineCouplingInOutVector(in,out);
}




//
//	defineCouplingInOutVector
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingInOutVector( string in, core::VectorStrings out )
{
Topology_sp	link;

    link = new_Topology_sp();
    link->setIn(in);
    link->setOut(out);
    this->links.push_back(link);
}

//
//	defineCouplingOut1
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingOut1( 	string out1 )
{
core::VectorStrings out;

    out.push_back(out1);
    this->defineCouplingOutVector(out);
}
//
//	defineCouplingOut2
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingOut2( 	string out1,
						string out2 )
{
core::VectorStrings out;

    out.push_back(out1);
    out.push_back(out2);
    this->defineCouplingOutVector(out);
}
//
//	defineCouplingOut3
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingOut3( 	string out1,
						string out2,
						string out3 )
{
core::VectorStrings out;

    out.push_back(out1);
    out.push_back(out2);
    out.push_back(out3);
    this->defineCouplingOutVector(out);
}





//
//	defineCouplingOutVector
//
//	Define a valid constellation of links for this Constitution
//
void	Constitution_O::defineCouplingOutVector( core::VectorStrings out )
{
Topology_sp	link;

    link = new_Topology_sp();
    link->setIn("");
    link->setOut(out);
    this->links.push_back(link);
}
//
//	isTopologyValid
//
//	Return true if this Constitution recognizes this coupling map
//
bool	Constitution_O::isTopologyValid( Topology_sp cm )
{_G();
    gctools::Vec0<Topology_sp>::iterator	i;

#ifdef	DEBUG_ON
    LOG(BF("==========  Comparing coupling map VVVVVVV") ); //
    cm->dump();
#endif

    for ( i=this->links.begin(); i!= this->links.end(); i++ ) {
#ifdef	DEBUG_ON
    LOG(BF("----------  To constitution coupling map") ); //
    (*i)->dump();
#endif

	if ( (*i)->matchesTopology(cm) ) {
	    return true;
            LOG(BF("      YES!!!") ); //
        }
    }
    LOG(BF("      Nope!!!") ); //
    return false;
}

#else //] [


#endif //]



void	Constitution_O::initialize()
{_OF();
    this->Base::initialize();
    this->_Comment = "";
    this->_Name = _Nil<core::Symbol_O>();
    this->_ConstitutionAtoms = ConstitutionAtoms_O::create();
    this->_StereoInformation = StereoInformation_O::create();
}

//
// Destructor
//




#ifdef XML_ARCHIVE
void	Constitution_O::archiveBase(core::ArchiveP node)
{_G();
    this->Base::archiveBase(node);
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
    node->attributeIfNotDefault<string>("comment",this->_Comment,"");
    node->attributeIfNotNil("metaConstitution",this->_MetaConstitution);
    node->attribute("constitutionAtoms",this->_ConstitutionAtoms);
    node->attribute("stereoInformation", this->_StereoInformation );
    node->archiveSymbolMap( "constitutionPlugs", this->_PlugsByName );
    node->archiveSymbolMap( "topologies", this->_Topologies );
//    node->archiveMap( "fragments", this->_Fragments );
//    node->archiveMap( "frames", this->_Frames );
}
#endif








void Constitution_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<Constitution_O>()
//	.def_raw("core:__init__",&Constitution_O::__init__,"(self &key name comment metaConstitution constitutionAtoms stereoInformation plugs topologies)")
//	.def("setCandoDatabase",&Constitution_O::setCandoDatabase)
//	.def("getCandoDatabase",&Constitution_O::getCandoDatabase)
	.def("stereoisomersAsCons", &Constitution_O::stereoisomersAsCons)
	.def("topologiesAsCons", &Constitution_O::topologiesAsCons)
	.def("plugsWithMatesAsCons", &Constitution_O::plugsWithMatesAsCons)
	.def("setConstitutionName", &Constitution_O::setConstitutionName)
	.def("constitutionName", &Constitution_O::constitutionName)
	.def("setMetaConstitutionName", &Constitution_O::setMetaConstitutionName)
	.def("getMetaConstitutionName", &Constitution_O::getMetaConstitutionName)
	.def("asGroup",&Constitution_O::asGroup)
	.def("getPlugNamed",&Constitution_O::getPlugNamed)
//	.def("getMonomerNameAsStringList", &Constitution_O::getMonomerNamesAsStringList) // this fails in to_object conversion
	.def("getMonomerNameAsStringSet", &Constitution_O::getMonomerNamesAsSymbolSet)
	.def("getPdbNamesAsStringList", &Constitution_O::getPdbNamesAsStringList)
	.def("setComment", &Constitution_O::setComment)
	.def("getComment", &Constitution_O::getComment)
	.def("hasPlugNamed",&Constitution_O::hasPlugNamed)
	.def("simplestTopologyWithPlugNamed",&Constitution_O::simplestTopologyWithPlugNamed)
	.def("getTopologyForContext",&Constitution_O::getTopologyForContext)
	.def("getConstitutionAtoms", &Constitution_O::getConstitutionAtoms)
	.def("createResidueForStereoisomerName", &Constitution_O::createResidueForStereoisomerName)
	.def("topologyWithName",&Constitution_O::topologyWithName)
	;
}

    void Constitution_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Constitution,"","",_lisp)
//	.def("setCandoDatabase",&Constitution_O::setCandoDatabase)
//	.def("getCandoDatabase",&Constitution_O::getCandoDatabase)
	.def("stereoisomersAsCons", &Constitution_O::stereoisomersAsCons)
	.def("topologiesAsCons", &Constitution_O::topologiesAsCons)
	.def("plugsWithMatesAsCons", &Constitution_O::plugsWithMatesAsCons)
	.def("setConstitutionName", &Constitution_O::setConstitutionName)
	.def("constitutionName", &Constitution_O::constitutionName)
	.def("setMetaConstitutionName", &Constitution_O::setMetaConstitutionName)
	.def("getMetaConstitutionName", &Constitution_O::getMetaConstitutionName)
	.def("asGroup",&Constitution_O::asGroup)
	.def("getPlugNamed",&Constitution_O::getPlugNamed)
//	.def("getMonomerNameAsStringList", &Constitution_O::getMonomerNamesAsStringList)  // This fails in to_object conversion
	.def("getMonomerNameAsStringSet", &Constitution_O::getMonomerNamesAsSymbolSet)
	.def("getPdbNamesAsStringList", &Constitution_O::getPdbNamesAsStringList)
	.def("setComment", &Constitution_O::setComment)
	.def("getComment", &Constitution_O::getComment)
	.def("hasPlugNamed",&Constitution_O::hasPlugNamed)
	.def("simplestTopologyWithPlugNamed",&Constitution_O::simplestTopologyWithPlugNamed)
	.def("getTopologyForContext",&Constitution_O::getTopologyForContext)
	.def("getConstitutionAtoms", &Constitution_O::getConstitutionAtoms)
	.def("createResidueForStereoisomerName", &Constitution_O::createResidueForStereoisomerName)
	.def("topologyWithName",&Constitution_O::topologyWithName)
	;
#endif
}




EXPOSE_CLASS(chem,Constitution_O);

}; // namespace chem

