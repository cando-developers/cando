       
//
// (C) 2004 Christian E. Schafmeister
//

#define	TURN_DEBUG_OFF

#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <clasp/core/stringSet.h>
#include <clasp/core/environment.h>
#include <clasp/core/symbolSet.h>
#include <clasp/core/stringList.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/constitution.h>
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/restraint.h>
#include <cando/chem/complexRestraints.h>
#include <cando/chem/entity.h>
#include <clasp/core/wrappers.h>




namespace chem {




/*! Return +1 for "R" and -1 for "S"
*/
int	StereoConfiguration_O::getMoeConfiguration()
{
    if ( this->_Configuration == "R" ) {
	return 1;
    }
	// If S
    return -1;
}


    SYMBOL_EXPORT_SC_(ChemKwPkg,R);
    SYMBOL_EXPORT_SC_(ChemKwPkg,S);

#ifdef XML_ARCHIVE
    void	StereoConfiguration_O::archive(core::ArchiveP node)
{_G();
    node->attribute("atom",this->_AtomName);
    node->attribute("config",this->_Configuration);
}
#endif


/*
__BEGIN_DOC(classes.StereoConfiguration.!class.StereoConfiguration)
\requiredKeyed{atomName:}{Text::atom}
\requiredKeyed{config:}{Text::configuration}

Provide the atom name and the stereo-configuration \sa{configuration} of "R" or "S".
__END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_StereoConfiguration_O_make "(atom_name config)"
#define DECL_StereoConfiguration_O_make ""
#define DOCS_StereoConfiguration_O_make "make StereoConfiguration"
StereoConfiguration_sp StereoConfiguration_O::make(string atomName, string config)
  {_G();
      GC_ALLOCATE(StereoConfiguration_O, me );
    me->_AtomName = atomName;
    me->_Configuration = config;
    return me;
  };

#else
    core::T_sp StereoConfiguration_O::__init__(core::Function_sp exec, core::Cons_sp args,
					      core::Environment_sp env, core::Lisp_sp lisp)
{_G();
    this->_AtomName = env->lookup(Pkg(),"atomName").as<core::Str_O>()->get();
    this->_Configuration = translate::from_object<string>::convert(env->lookup(Pkg(),"config"));
    return _Nil<core::T_O>();
}

#endif

/*! Create a list of StereoConfiguration(s)
 * Match atom names to configurations
 */
core::Cons_sp StereoConfiguration_O::create_multiple(core::Cons_sp atomNames, core::Cons_sp configurations, core::Lisp_sp lisp )
{_G();
    ASSERT_eq(atomNames->length(), configurations->length() );
    core::Cons_sp list = _Nil<core::Cons_O>();
    core::Cons_sp curName = atomNames;
    core::Cons_sp curConfig = configurations;
    while ( curName.notnilp() )
    {
	GC_ALLOCATE(StereoConfiguration_O, one );
	one->setAtomName(curName->car<core::Str_O>()->get());
	one->setConfiguration(curConfig->car<core::Str_O>()->get());
	list = core::Cons_O::create(one,list,lisp);
	curName = curName->cdr();
	curConfig = curConfig->cdr();
    }
    return list;
}



    
    
#define ARGS_af_StereoConfiguration_create_multiple "(atomNames configurations)"
#define DECL_af_StereoConfiguration_create_multiple ""
#define DOCS_af_StereoConfiguration_create_multiple "StereoConfiguration_create_multiple"
    core::T_sp af_StereoConfiguration_create_multiple(core::Cons_sp atomNames, core::Cons_sp configurations)
    {_G();
	return StereoConfiguration_O::create_multiple(atomNames,configurations,_lisp);
    }

core::Cons_sp StereoConfiguration_O::stereochemicalPermutations(uint numberOfCenters,core::Lisp_sp lisp)
{_G();
    core::Cons_sp list = _Nil<core::Cons_O>();
    uint pow = 1 << numberOfCenters;
    for ( uint i=0; i<pow; i++ )
    {
	uint map = i;
	core::Cons_sp one = _Nil<core::Cons_O>();
	for ( uint ix=0; ix<numberOfCenters; ix++ )
	{
	    if ( map & 1 )
	    {
		one = core::Cons_O::create(core::Str_O::create("R"),one,lisp);
	    } else
	    {
		one = core::Cons_O::create(core::Str_O::create("S"),one,lisp);
	    }
	    map = map >> 1;
	}
	list = core::Cons_O::create(one,list,lisp);
    }
    return list;
}




    
    
#define ARGS_af_stereochemicalPermutations "(numcenters)"
#define DECL_af_stereochemicalPermutations ""
#define DOCS_af_stereochemicalPermutations "stereochemicalPermutations"
    core::T_sp af_stereochemicalPermutations(core::Fixnum_sp numcenters)
    {_G();
	uint numberOfCenters = numcenters->get();
	core::Cons_sp results = StereoConfiguration_O::stereochemicalPermutations(numberOfCenters,_lisp);
	return results;
    }





void	Stereoisomer_O::initialize()
{_OF();
    this->Base::initialize();
    this->_Name = _Nil<core::Symbol_O>();
    this->_Pdb = _Nil<core::Symbol_O>();
    this->_WeakConstitution = _Nil<Constitution_O>();
    this->_Enantiomer = _Nil<core::Symbol_O>();
}


#ifdef XML_ARCHIVE
    void	Stereoisomer_O::archive(core::ArchiveP node)
{
    if ( node->saving() ) this->getConstitution();
    this->Base::archiveBase(node);
    node->archiveWeakPointer("constitution",this->_WeakConstitution);
    node->attribute("name",this->_Name );
    node->attribute("pdb", this->_Pdb );
    node->attributeIfNotNil("enantiomer", this->_Enantiomer);
    node->archiveVector0OfObjectsSubClassOf<StereoConfiguration_O>( this->_Configurations );
}
#endif

core::SymbolSet_sp	Stereoisomer_O::expandedNameSet()
{_G();
core::SymbolSet_sp	ss;
    ss = core::SymbolSet_O::create();
    ss->insert(this->getName());
    return ss;
}

	/*! Build and return an ObjectList that contains one RepresentedEntityNameSet
	 * that has this Stereoisomers name as its representative and
	 * a ContainedNames list with only this Stereoisomers name
	 */
RepresentativeList_sp Stereoisomer_O::expandedRepresentativeList() const
{_G();
    RepresentativeList_sp ss = RepresentativeList_O::create();
    RepresentedEntityNameSet_sp gr = RepresentedEntityNameSet_O::create();
    gr->setName(this->getName());
    gr->setRepresentative(this->getName());
    gr->addEntityName(this->getName());
    ss = RepresentativeList_O::create();
    ss->append(gr);
    return ss;
}


string Stereoisomer_O::__repr__() const
{
    stringstream ss;
    ss << "Stereoisomer-name(" << this->_Name->__repr__() << ")-pdbName(" << this->_Pdb->__repr__() << ")";
    return ss.str();
}

string Stereoisomer_O::getConfigurationForCenter( const string& centerName )
{_G();
    for ( gctools::Vec0<StereoConfiguration_sp>::iterator it=this->_Configurations.begin();
          it!=this->_Configurations.end(); it++ )
    {
	if ( (*it)->getAtomName() == centerName )
	{
	    return (*it)->getConfiguration();
	}
    }
    return "unknownConfig";
}

Constitution_sp Stereoisomer_O::getConstitution()
{_G();
    ASSERTNOTNULL(this->_WeakConstitution);
    		// If the constitution is undefined then
		// crawl up the owners to find one
    if ( this->_WeakConstitution.lock().nilp() )
    {
	IMPLEMENT_MEF(BF("Handle ownerWithClass"));
//	this->_WeakConstitution = this->ownerWithClass<Constitution_O>();
    }
    return this->_WeakConstitution.lock();
}


/*
__BEGIN_DOC(classes.Stereoisomer.!class.Stereoisomer)
\scriptCmdRet{Stereoisomer}{name:(Text::name) pdb:(Text::pdb) configs:(List::)}{StereoIsomer::}
__END_DOC
*/

#if INIT_TO_FACTORIES

#define ARGS_Stereoisomer_O_make "(name pdb configs)"
#define DECL_Stereoisomer_O_make ""
#define DOCS_Stereoisomer_O_make "make ChemDraw"
Stereoisomer_sp Stereoisomer_O::make(core::Symbol_sp name, core::Symbol_sp pdb, core::Cons_sp configs)
  {_G();
      GC_ALLOCATE(Stereoisomer_O, me );
    me->_Name = name;
    me->_Pdb = pdb;
    IMPLEMENT_MEF(BF("Handle setOwnerOfAllEntries"));
//    configs->setOwnerOfAllEntries(me);
    me->_Configurations.clear();
    for ( ; configs.notnilp(); configs = configs->cdr() )
      {
	me->_Configurations.push_back(configs->car<StereoConfiguration_O>() );
      }
    return me;
  };

#else

    core::T_sp Stereoisomer_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{_G();
    this->_Name = translate::from_object<core::Symbol_O>::convert(env->lookup(Pkg(),"name"));
    this->_Pdb = translate::from_object<core::Symbol_O>::convert(env->lookup(Pkg(),"pdb"));
    core::Cons_sp configs = translate::from_object<core::Cons_O>::convert(env->lookup(Pkg(),"configs"));
    configs->setOwnerOfAllEntries(this->sharedThis<Stereoisomer_O>());
    this->_Configurations.clear();
    for ( ; configs.notnilp(); configs = configs->cdr() )
    {
	this->_Configurations.push_back(configs->car<StereoConfiguration_O>() );
    }
    return _Nil<core::T_O>();
}

#endif

/*
__BEGIN_DOC(classes.Stereoisomer.!class.MultiStereoisomers)
\scriptCmdRet{MultiStereoisomers}{nameTemplate:(Text::template) centers:(List::) configs:(List::)}{Cons::stereoisomers}
__END_DOC
*/




#define ARGS_af_multiStereoisomers "(name-template centers configs)"
#define DECL_af_multiStereoisomers ""
#define DOCS_af_multiStereoisomers "multiStereoisomers"
core::T_sp af_multiStereoisomers(const string& nameTemplate, core::Cons_sp centers, core::Cons_sp configs)
{_G();
    IMPLEMENT_MEF(BF("Implement multiStereoisomers"));
};





void	StereoInformation_O::addStereoisomer(Stereoisomer_sp s)
{_OF();
    if ( this->_NameOrPdbToStereoisomer.count(s->getName())>0 ) {
	SIMPLE_ERROR(BF("addStereoisomer monomer name (%s) has already been used") % s->getName());
    }
    if ( this->_NameOrPdbToStereoisomer.count(s->getPdb())>0 ) {
	SIMPLE_ERROR(BF("addStereoisomer pdb name (%s) has already been used") % s->getPdb() );
    }
    this->_NameOrPdbToStereoisomer.set(s->getName(),s);
    this->_NameOrPdbToStereoisomer.set(s->getPdb(),s);
    this->_Stereoisomers.append(s);
}


void StereoInformation_O::validate()
{_OF();
    for ( gctools::Vec0<Stereoisomer_sp>::iterator it=this->_Stereoisomers.begin(); it!=this->_Stereoisomers.end(); it++ )
    {
	if ( (*it)->getName().nilp() ) 
	{
	    SIMPLE_ERROR(BF("StereoInformation has stereoisomer with blank name"));
	}
	if ( (*it)->getPdb().nilp() ) 
	{
	    SIMPLE_ERROR(BF("StereoInformation has stereoisomer with blank pdb name"));
	}
    }
}



#if 0
/*
__BEGIN_DOC(classes.StereoInformation.!instance.addProChiralCenter)
\scriptMethod{self}{addProChiralCenter}{ProChiralCenter::center}

Add the \sa{center} to the StereoInformation object.
__END_DOC
*/
void	StereoInformation_O::addProChiralCenter(RPProChiralCenter s)
{
    this->_ProChiralCenters.push_back(s);
}
#endif


    Stereoisomer_sp	StereoInformation_O::getStereoisomer(core::Symbol_sp  nameOrPdb )
    {_OF();
    return this->_NameOrPdbToStereoisomer.get(nameOrPdb);
}

#ifdef XML_ARCHIVE
void	StereoInformation_O::archive(core::ArchiveP node)
{_G();
    node->archiveVector0OfObjectsSubClassOf<Stereoisomer_O>(this->_Stereoisomers);
//    node->archiveVector0OfObjectsSubClassOf<O_ProChiralCenter>( this->_ProChiralCenters);
    node->archiveVector0( "complexRestraints", this->_ComplexRestraints);
    node->needsFinalization();
}
#endif

#ifdef XML_ARCHIVE
bool	StereoInformation_O::loadFinalize(core::ArchiveP node)
{_G();
    gctools::Vec0<Stereoisomer_sp>::iterator	si;
    for ( si=this->_Stereoisomers.begin(); si!=this->_Stereoisomers.end(); si++ )
    {
    	LOG(BF("Adding stereoisomer with name: %s") % (*si)->getName().c_str()  );
        this->_NameOrPdbToStereoisomer.set((*si)->getName(), (*si));
        this->_NameOrPdbToStereoisomer.set((*si)->getPdb(), (*si));
    }
    return true;
}
#endif


/*
__BEGIN_DOC(classes.StereoInformation.!class.StereoInformation)
\requiredKeyed{stereoisomers:}{List::stereoisomers}
\optionalKeyed{proChiralCenters:}{List::}
\optionalKeyed{constrainedPiBonds:}{List::}
__END_DOC
*/


#if INIT_TO_FACTORIES

#define ARGS_StereoInformation_O_make "(stereoisomers restraints)"
#define DECL_StereoInformation_O_make ""
#define DOCS_StereoInformation_O_make "make StereoInformation"
StereoInformation_sp StereoInformation_O::make(core::Cons_sp stereoisomers, core::Cons_sp restraints)
  {_G();
      GC_ALLOCATE(StereoInformation_O, me );
    me->_Stereoisomers.fillFromCons(stereoisomers);
    IMPLEMENT_MEF(BF("Handle setOwnerOfAllEntries"));
//    stereoisomers->setOwnerOfAllEntries(me);
    me->_NameOrPdbToStereoisomer.clear();
    for ( gctools::Vec0<Stereoisomer_sp>::iterator it = me->_Stereoisomers.begin();
	  it != me->_Stereoisomers.end(); it++ )
      {
	if ( me->_NameOrPdbToStereoisomer.count((*it)->getName()) >0 )
	  {
	      SIMPLE_ERROR(BF("Stereoisomer name: "+(*it)->getName()->__repr__()+" has already been defined"));
	  }
	if ( me->_NameOrPdbToStereoisomer.count((*it)->getPdb()) >0 )
	  {
	      SIMPLE_ERROR(BF("Stereoisomer pdb: "+(*it)->getPdb()->__repr__()+" has already been defined"));
	  }
	me->_NameOrPdbToStereoisomer.set((*it)->getName(), (*it));
	me->_NameOrPdbToStereoisomer.set((*it)->getPdb(), (*it));
      }
    return me;
  };

#else

    core::T_sp 	StereoInformation_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{_OF();
    this->Base::__init__(exec,args,env,lisp);
    core::Cons_sp stereoisomers = env->lookup(Pkg(),"stereoisomers").as<core::Cons_O>();
    core::Cons_sp restraints = env->lookup(Pkg(),"restraints").as<core::Cons_O>();
    this->_Stereoisomers.fillFromCons(stereoisomers);
    stereoisomers->setOwnerOfAllEntries(this->sharedThis<StereoInformation_O>());
    this->_NameOrPdbToStereoisomer.clear();
    for ( gctools::Vec0<Stereoisomer_sp>::iterator it = this->_Stereoisomers.begin();
    		it != this->_Stereoisomers.end(); it++ )
    {
	if ( this->_NameOrPdbToStereoisomer.count((*it)->getName()) >0 )
	{
	    SIMPLE_ERROR(BF("Stereoisomer name: "+(*it)->getName()->__repr__()+" has already been defined"));
	}
	if ( this->_NameOrPdbToStereoisomer.count((*it)->getPdb()) >0 )
	{
	    SIMPLE_ERROR(BF("Stereoisomer pdb: "+(*it)->getPdb()->__repr__()+" has already been defined"));
	}
	this->_NameOrPdbToStereoisomer.set((*it)->getName(), (*it));
	this->_NameOrPdbToStereoisomer.set((*it)->getPdb(), (*it));
	return _Nil<core::T_O>();
    }

//    core::Cons_sp proChiralCenters = dict->getAndRemoveOrDefault("proChiralCenters",_Nil<core::Cons_O>()).as<core::Cons_O>();
//    proChiralCenters->setOwnerOfAllEntries(this->sharedThis<StereoInformation_O>());
//    this->_ProChiralCenters.fillFromCons(proChiralCenters);
    restraints->setOwnerOfAllEntries(this->sharedThis<StereoInformation_O>());
    this->_ComplexRestraints.fillFromCons(restraints);
    return _Nil<core::T_O>();
}
#endif



//
//	getMonomerNames
//
//!	Return a core::StringList_sp of monomer names for this mold
core::StringList_sp StereoInformation_O::getMonomerNamesAsStringList()
{
core::StringList_sp				names;
gctools::Vec0<Stereoisomer_sp>::iterator	mnpi;
    names = core::StringList_O::create();
    for ( mnpi=this->_Stereoisomers.begin();
		mnpi!=this->_Stereoisomers.end(); mnpi++ ) {
	names->append((*mnpi)->getName()->symbolNameAsString());
    }
    return names;
}




//!	Return a core::StringSet_sp of monomer names for this Constitution
core::SymbolSet_sp StereoInformation_O::getMonomerNamesAsSymbolSet()
{
core::SymbolSet_sp				names;
gctools::Vec0<Stereoisomer_sp>::iterator	mnpi;
    names = core::SymbolSet_O::create();
    for ( mnpi=this->_Stereoisomers.begin();
		mnpi!=this->_Stereoisomers.end(); mnpi++ ) {
	names->insert((*mnpi)->getName());
    }
    return names;
}



//
//	getPdbNames
//
//!	Return a core::VectorStrings of monomer pdb names for this mold
core::StringList_sp StereoInformation_O::getPdbNamesAsStringList()
{
core::StringList_sp	names;
gctools::Vec0<Stereoisomer_sp>::iterator	mnpi;
    names = core::StringList_O::create();
    for ( mnpi=this->_Stereoisomers.begin();
		mnpi!=this->_Stereoisomers.end(); mnpi++ ) {
	names->append((*mnpi)->getPdb()->symbolNameAsString());
    }
    return names;
}

    void StereoConfiguration_O::lisp_initGlobals(core::Lisp_sp lisp)
    {_G();
	core::enum_<ConfigurationEnum>(_sym__PLUS_configurationEnumConverter_PLUS_,"Configuration")
	    .value(_lisp->intern("undefinedConfiguration"), _lisp->intern("U"), undefinedConfiguration )
	    .value(_lisp->intern("S"), S_Configuration )
	    .value(_lisp->intern("R"), R_Configuration )
	    ;
	core::enum_<StereochemistryType>(_sym__PLUS_stereochemistryTypeConverter_PLUS_,"Stereochemistry")
	    .value(_lisp->intern("UNDEFINED-CHIRALITY"), _lisp->intern("u"),undefinedCenter)
	    .value(_lisp->intern("CHIRAL"), _lisp->intern("C"), chiralCenter)
	    .value(_lisp->intern("PROCHIRAL"), _lisp->intern("P"), prochiralCenter)
	    ;
    }

void StereoConfiguration_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<StereoConfiguration_O>()
//	    .def_raw("core:__init__",&StereoConfiguration_O::__init__,"(self &key (atomName \"\") (config \"\"))")
	.def("getAtomName",&StereoConfiguration_O::getAtomName)
	.def("setAtomName",&StereoConfiguration_O::setAtomName)
	.def("getConfiguration",&StereoConfiguration_O::getConfiguration)
	.def("setConfiguration",&StereoConfiguration_O::setConfiguration)
	.def("getMoeConfiguration",&StereoConfiguration_O::getMoeConfiguration)
	;
	Defun(stereochemicalPermutations);
	Defun(StereoConfiguration_create_multiple);
//	def("stereochemicalPermutations",&StereoConfiguration_O::stereochemicalPermutations,lisp->lisp());

    }

    void StereoConfiguration_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,StereoConfiguration,"","",_lisp)
	.def("getAtomName",&StereoConfiguration_O::getAtomName)
	.def("setAtomName",&StereoConfiguration_O::setAtomName)
	.def("getConfiguration",&StereoConfiguration_O::getConfiguration)
	.def("setConfiguration",&StereoConfiguration_O::setConfiguration)
	.def("getMoeConfiguration",&StereoConfiguration_O::getMoeConfiguration)
	;
#endif
//	defNoWrapPackage(ChemPkg,"stereochemicalPermutations",&prim_stereochemicalPermutations,lisp->lisp());
//	defNoWrapPackage(ChemPkg,"StereoConfiguration_createMultiple",&prim_StereoConfiguration_create_multiple,lisp->lisp());
//	def("stereochemicalPermutations",&StereoConfiguration_O::stereochemicalPermutations,lisp->lisp());

    }

;



void Stereoisomer_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Stereoisomer_O>()
//	    .def_raw("core:__init__",&Stereoisomer_O::__init__,"(self &key name pdb configs)")
	    .def("getName",&Stereoisomer_O::getName,"","","",false)
	    .def("setName",&Stereoisomer_O::setName,"","","",false)
	.def("getPdb",&Stereoisomer_O::getPdb)
	.def("setPdb",&Stereoisomer_O::setPdb)
	.def("getConfigurationForCenter",&Stereoisomer_O::getConfigurationForCenter)
	;
    }

    void Stereoisomer_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Stereoisomer,"","",_lisp)
	.def("getName",&Stereoisomer_O::getName)
	.def("setName",&Stereoisomer_O::setName)
	.def("getPdb",&Stereoisomer_O::getPdb)
	.def("setPdb",&Stereoisomer_O::setPdb)
	.def("getConfigurationForCenter",&Stereoisomer_O::getConfigurationForCenter)
	;
#endif
    }

;

void StereoInformation_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<StereoInformation_O>()
//	    .def_raw("core:__init__",&StereoInformation_O::__init__,"(self &key stereoisomers)")
	.def("stereoisomersAsCons",&StereoInformation_O::stereoisomersAsCons)
	.def("validate",&StereoInformation_O::validate)
	;
    }

    void StereoInformation_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,StereoInformation,"","",_lisp)
	.def("stereoisomersAsCons",&StereoInformation_O::stereoisomersAsCons)
	.def("validate",&StereoInformation_O::validate)
	;
#endif
    }

;


    EXPOSE_CLASS_AND_GLOBALS(chem, StereoConfiguration_O );
    EXPOSE_CLASS(chem, Stereoisomer_O);
    EXPOSE_CLASS(chem, StereoInformation_O );



}; //namespace chem



