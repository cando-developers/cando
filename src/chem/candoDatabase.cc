/*
    File: candoDatabase.cc
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


#include <cando/chem/candoDatabase.h>

#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <set>
#include <stdio.h>
#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/environment.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/matter.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/specificContext.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/loop.h>
#include <cando/chem/representedEntityNameSet.h>
//#include "monomerCoordinates.h"
//#include <clasp/core/objRef.h>
//#include "core/xmlSaveArchive.h"
#include <cando/chem/coupling.h>
#include <cando/chem/frameRecognizer.h>
#include <cando/chem/monomerPack.h>
#include <cando/chem/constitution.h>
#include <cando/chem/monomerCoordinates.h>
#include <cando/chem/topology.h>
#include <cando/chem/plug.h>
#include <cando/chem/specificContext.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/candoDatabaseDependent.h>
#include <clasp/core/wrappers.h>


namespace chem 
{


    /*!Store a list of all CandoDatabases so that they can be consolidated
     */
//    core::WeakVector0<CandoDatabase_O>	AllCandoDatabases;




    void	CandoDatabase_O::initialize()
    {_OF();
	this->Base::initialize();
	this->_Entities.clear();
	LOG(BF("initialized _Entities@%p") % &(this->_Entities)  );
//    AllCandoDatabases.append(this->sharedThis<CandoDatabase_O>());
    }


    core::List_sp CandoDatabase_O::entitiesSubClassOfAsCons(core::Class_sp mc)
    {_OF();
	core::List_sp result = _Nil<core::T_O>();
	for ( Entities::const_iterator it=this->_Entities.begin();
	      it!=this->_Entities.end(); it++ )
	{
//	    if ( it->second->isAssignableToByClassSID(classSID) )
	    if ( it->second->isAInstanceOf(mc) )
	    {
		result = core::Cons_O::create(it->second,result);
	    }
	}
	return result;
    }


    bool CandoDatabase_O::recognizesEntitySubClassOf(core::Symbol_sp name, core::Class_sp mc)
    {_OF();
	Entities::const_iterator it = this->_Entities.find(name);
	if ( it == this->_Entities.end() ) return false;
	if ( it->second->isAInstanceOf(mc) ) return true;
	return false;
    }

    bool CandoDatabase_O::recognizesEntityOfClass(core::Symbol_sp name, core::Class_sp aclass)
    {_OF();
	Entities::iterator it = this->_Entities.find(name);
	if ( it == this->_Entities.end() ) return false;
	if ( it->second->__class() == aclass) return true;
	return false;
    }




CL_LISPIFY_NAME("constitutionsAsCons");
CL_DEFMETHOD     core::List_sp	CandoDatabase_O::constitutionsAsCons() 
    {_OF();
      return this->entitiesSubClassOfAsCons(core::cl__find_class(Constitution_O::static_classSymbol()));
    };


    bool CandoDatabase_O::recognizesMonomerName(core::Symbol_sp name)
    {_OF();
	return this->recognizesEntitySubClassOf(name,core::cl__find_class(Stereoisomer_O::static_classSymbol()));
    }


CL_LISPIFY_NAME("addEntity");
CL_DEFMETHOD     Entity_sp CandoDatabase_O::addEntity(Entity_sp entity)
    {
	this->_Entities.set(entity->getName(),entity);
	if ( entity->isAssignableTo<Constitution_O>() )
	{
	    entity.as<Constitution_O>()->addStereoisomersToCandoDatabase(this->sharedThis<CandoDatabase_O>());
	}
	return entity;
    }


CL_LISPIFY_NAME("representedEntityNameSetsAsCons");
CL_DEFMETHOD     core::List_sp CandoDatabase_O::representedEntityNameSetsAsCons()
    {_OF();
	return this->entitiesSubClassOfAsCons(core::cl__find_class(RepresentedEntityNameSet_O::static_classSymbol()));
    }







CL_LISPIFY_NAME("recognizesRepresentedEntityNameSet");
CL_DEFMETHOD     bool CandoDatabase_O::recognizesRepresentedEntityNameSet(core::Symbol_sp name)
    {_OF();
	return this->recognizesEntityOfClass(name,core::cl__find_class(RepresentedEntityNameSet_O::static_classSymbol()));
    }


    bool CandoDatabase_O::recognizesMonomerPack(core::Symbol_sp name)
    {_OF();
	if ( this->_Entities.count(name) == 0 ) return false;
	Entity_sp entity = this->_Entities.get(name);
	return entity->isAssignableTo<MonomerPack_O>();
    }


    SYMBOL_EXPORT_SC_(ChemKwPkg,constitution);
    SYMBOL_EXPORT_SC_(ChemKwPkg,entity);
    SYMBOL_EXPORT_SC_(ChemKwPkg,frameRecognizer);

    core::T_sp	CandoDatabase_O::oGetResource(core::Symbol_sp resource, core::Symbol_sp name )
    {_OF();
	if ( resource == chemkw::_sym_constitution )
	{
	    return this->constitutionForNameOrPdb(name);
	}
	if ( resource == chemkw::_sym_entity )
	{
	    return this->getEntity(name);
	}
	if ( resource == chemkw::_sym_frameRecognizer )
	{
	    return this->getFrameRecognizer(name);
	}
	SIMPLE_ERROR(BF("Unknown CandoDatabase Resource %s %s") % _rep_(resource) % _rep_(name));
    }



CL_LISPIFY_NAME("addFrameRecognizer");
CL_DEFMETHOD     void CandoDatabase_O::addFrameRecognizer(FrameRecognizer_sp rec)
    { 
	this->_frameRecognizers.set(rec->getRecognizerName(),rec);
    }


CL_LISPIFY_NAME("recognizesFrameRecognizerName");
CL_DEFMETHOD     bool	CandoDatabase_O::recognizesFrameRecognizerName(core::Symbol_sp nm)
    { 
	bool rec;
	rec = this->_frameRecognizers.contains(nm);
	LOG(BF("found = %d") % rec  );
	return rec;
    }




#if 0
    core::T_sp CandoDatabase_O::oGetReference(core::ObjRef_sp ref)
    {
	if ( ref->getSelector() == "Constitution" )
	{
	    return this->_Constitutions.getDefaultNil(ref->getName());
	}
	if ( ref->getSelector() == "FrameRecognizer" )
	{
	    return this->getFrameRecognizer(ref->getName());
	}
	SIMPLE_ERROR(BF("Unknown ObjRef selector: "+ref->getSelector() ));
    }
#endif



    bool	CandoDatabase_O::recognizesEntityNameSetName(core::Symbol_sp nm)
    {
	return this->recognizesEntityOfClass(nm,core::cl__find_class(EntityNameSet_O::static_classSymbol()));
    }

CL_LISPIFY_NAME("getEntityNameSet");
CL_DEFMETHOD     EntityNameSet_sp CandoDatabase_O::getEntityNameSet( core::Symbol_sp nm )
    {
	return this->getEntityOfClass(nm,core::cl__find_class(EntityNameSet_O::static_classSymbol())).as<EntityNameSet_O>();
    };



CL_LISPIFY_NAME("getRepresentedEntityNameSet");
CL_DEFMETHOD     RepresentedEntityNameSet_sp CandoDatabase_O::getRepresentedEntityNameSet( core::Symbol_sp nm )
    {
	return this->getEntityOfClass(nm,core::cl__find_class(RepresentedEntityNameSet_O::static_classSymbol())).as<RepresentedEntityNameSet_O>();
    };



CL_LISPIFY_NAME("recognizesSetOrConstitutionOrMonomerName");
CL_DEFMETHOD     bool	CandoDatabase_O::recognizesSetOrConstitutionOrMonomerName(core::Symbol_sp nm)
    {_OF();
	return this->_Entities.count(nm)>0;
    }

CL_LISPIFY_NAME("getEntity");
CL_DEFMETHOD     Entity_sp CandoDatabase_O::getEntity(core::Symbol_sp nm) const
    {
	ASSERTF(this->_Entities.contains(nm),
		BF("CandoDatabase does not contain entity[%s]") % nm->__repr__() );
	return this->_Entities.get(nm);
    }


    Entity_sp CandoDatabase_O::getEntityOfClass(core::Symbol_sp nm, core::Class_sp mc)
    {
	ASSERTF(this->_Entities.contains(nm),
		BF("CandoDatabase does not contain entity[%s]")
		% nm->__repr__() );
	Entity_sp entity = this->_Entities.get(nm);
	if ( entity->__class() != mc )
	{
	    SIMPLE_ERROR(BF("You asked for Entity[%s] of class[%s] and there is one with that name but it has the wrong class[%s]")
			       % nm->__repr__() % mc->className() % entity->className() );
	}
	return entity;
    }

CL_LISPIFY_NAME("expandEntityNameToTerminals");
CL_DEFMETHOD     adapt::SymbolSet_sp CandoDatabase_O::expandEntityNameToTerminals(core::Symbol_sp nm) const
    {_OF();
	Entity_sp entity = this->getEntity(nm);
	return entity->expandedNameSet();
    }

CL_LISPIFY_NAME("expandEntityNamesToTerminals");
CL_DEFMETHOD     adapt::SymbolSet_sp CandoDatabase_O::expandEntityNamesToTerminals(adapt::SymbolSet_sp nm) const
    {_OF();
	adapt::SymbolSet_sp names = adapt::SymbolSet_O::create();
        nm->map( [this,&names] (core::Symbol_sp sym) {
                Entity_sp entity = this->getEntity(sym);
                names->insertSymbolSet(entity->expandedNameSet());
            } );
	return names;
    }


#if 0
    bool	CandoDatabase_O::monomersAreDefinedForSetOrConstitutionOrMonomerName(core::Symbol_sp nm)
    {
	LOG(BF("looking for(%s)") % nm->__repr__()  );
	EntityNameSet_sp group;
	if ( this->recognizesEntityNameSetName(nm) )
	{
	    LOG(BF("Found group") );
	    group = this->getEntityNameSet(nm);
	    return group->isFinalized();
	}
    	// For constitutions and monomer names their monomer names are always defined
	return true;
    }
#endif

CL_LISPIFY_NAME("getMonomersForSetOrConstitutionOrMonomerName");
CL_DEFMETHOD     adapt::SymbolSet_sp CandoDatabase_O::getMonomersForSetOrConstitutionOrMonomerName(core::Symbol_sp nm)
    { _OF();
	// Use expandEntityNameToTerminals
	DEPRECIATED();
#if 0
	adapt::SymbolSet_sp	names;
	Entity_sp obj;

#ifdef	DEBUG_ON

	if ( !this->_Entities.contains(nm) )
	{
	    LOG(BF("The namespace contains %d names: ") % this->_Entities.size()  );
//        for ( adapt::SymbolMap<CandoDatabase_ODependent>::iterator it=this->_Entities.begin();
	    for ( adapt::SymbolMap<Entity_O>::iterator it=this->_Entities.begin();
		  it!=this->_Entities.end();
		  it++ )
	    {
		LOG(BF("Namespace contains name: %s") % it->first->__repr__() );
	    }
	}
#endif
//    ASSERTP(this->_Entities.contains(nm),"CandoDatabase Environment does not contain: "+nm->__repr__());
	if ( this->_Entities.contains(nm) )
	{
	    obj = this->_Entities.get(nm);
	    names = obj->expandedNameSet();
	} else
	{
	    names = _Nil<adapt::SymbolSet_O>();
	}
	return names;
#endif
    }

#if 0
    RepresentativeList_sp CandoDatabase_O::expandEntityNameToListOfRepresentatives(core::Symbol_sp nm)
    { 
	Entity_sp obj;
	geom::ObjectList_sp 	ll;
	LOG(BF("Looking for name(%s)") % nm->__repr__()  );
	ASSERT(this->recognizesEntityName(nm));
	obj = this->getEntity(nm);
	ll = obj->expandedRepresentativeList();
	return ll;
    }
#endif



CL_LISPIFY_NAME("recognizesNameOrPdb");
CL_DEFMETHOD     bool	CandoDatabase_O::recognizesNameOrPdb(core::Symbol_sp nm)
    {
	if ( !this->recognizesMonomerName(nm) )
//    ASSERTNOTNULL(this->_Names);
//    ASSERTvalue(this->_Entities.contains(nm), nm); // REMOVE
	if ( !this->_Entities.contains(nm) ) return false;
	Entity_sp obj;
	obj = this->_Entities.get(nm);
	if ( obj->isTerminalName() ) return true;
	LOG(BF("The name(%s) is in the namespace but its not a terminal") % nm);
	return false;
    };


#ifdef XML_ARCHIVE
    void	CandoDatabase_O::archiveBase(core::ArchiveP node)
    {
	this->Base::archiveBase(node);
	node->attribute("name",this->_Name);
	node->attributeIfNotDefault<string>("dateCreated",this->_DateCreated,"");
	node->attributeIfNotDefault<string>("dateUpdated",this->_DateUpdated,"");
	node->archiveSymbolMap( "entities", this->_Entities);
	node->archiveSymbolMap("frameRecognizers",this->_frameRecognizers);
//	node->archiveMapOfObjectsSubClassOf( this->_MonomerCoordinates );
	node->needsFinalization();
    }
#endif



CL_LISPIFY_NAME("constitutionForNameOrPdb");
CL_DEFMETHOD     Constitution_sp	CandoDatabase_O::constitutionForNameOrPdb(core::Symbol_sp name)
    {
	Entity_sp obj;
	ASSERTP(this->_Entities.contains(name),"Could not find: "+name->__repr__()+" in Environment");
	obj = this->_Entities.get(name);
	ASSERTP(obj->hasConstitution(), "The name: "+name->__repr__()+" is not part of a Constitution" );
	return obj->constitution();
    }

CL_LISPIFY_NAME("constitutionNameForNameOrPdb");
CL_DEFMETHOD     core::Symbol_sp CandoDatabase_O::constitutionNameForNameOrPdb(core::Symbol_sp name)
    {
	Constitution_sp con;
	con = this->constitutionForNameOrPdb(name);
	return con->getName();
    }

    core::Symbol_sp CandoDatabase_O::getMonomerNameForNameOrPdb(core::Symbol_sp name)
    {
	Entity_sp dep;
	ASSERTP(this->_Entities.contains(name), "Could not find name("+name->__repr__()+") in Environment");
	dep = this->_Entities.get(name);
	ASSERTP(dep->isTerminalName(), "Asked for MonomerName of non-terminal name("+name->__repr__()+")");
	return dep->getName();
    }

    core::Symbol_sp CandoDatabase_O::getPdbNameForNameOrPdb(core::Symbol_sp name)
    {
	Entity_sp 			dep;
	Stereoisomer_sp			si;
	ASSERTP(this->_Entities.contains(name), "Could not find name("+name->__repr__()+") in Environment");
	dep = this->_Entities.get(name);
	ASSERTP(dep->isTerminalName(), "Asked for MonomerName of non-terminal name("+name->__repr__()+")");
	si = downcast<Stereoisomer_O>(dep);
	return si->getPdb();
    }





#if 0
//
//	allMonomerNamesAsStringSet
//
//	Return a StringSet with all of the names of Monomers
//
    adapt::SymbolSet_sp CandoDatabase_O::allMonomerNamesAsSymbolSet()
    {
	adapt::SymbolMap<Constitution_O>::iterator	susi;
	adapt::SymbolSet_sp	names;
	adapt::SymbolSet_sp	monomers;
	names = adapt::SymbolSet_O::create();
	for ( susi=this->_Constitutions.begin();
	      susi!=this->_Constitutions.end(); susi++ ) {
	    monomers = susi->second->getMonomerNamesAsSymbolSet();
	    names->insertSymbolSet(monomers);
	}
	return names;
    }
#endif


#if 0
//
//	allMonomerNames
//
//	Return a StringList with all of the names of Monomers
//
    adapt::StringList_sp CandoDatabase_O::allMonomerNamesAsStringList()
    {
	adapt::SymbolMap<Constitution_O>::iterator	susi;
	adapt::StringList_sp	names;
	adapt::StringList_sp	monomers;
	names = adapt::StringList_O::create();
	for ( susi=this->_Constitutions.begin();
	      susi!=this->_Constitutions.end(); susi++ )
	{
	    monomers = susi->second->getMonomerNamesAsStringList();
	    names->appendStringList(monomers);
	}
	return names;
    }





/*!
  Add the Constitution to the CandoDatabase
*/
    void	CandoDatabase_O::addConstitution( Constitution_sp sum )
    {
//    sum->setCandoDatabase(this->sharedThis<CandoDatabase_O>());
//    this->_Names->fillFromConstitution(sum);
//    sum->setCandoDatabase(this->sharedThis<CandoDatabase_O>());
	LOG(BF("Adding constitution(%s)") % sum->getName() );
	this->_Constitutions.set(sum->getName(),sum);
	LOG(BF("CandoDatabase_Entities_add(%s) a Constitution") % (sum)->getName().c_str()  );
	this->_Entities.addUnique((sum)->getName(), sum);
	Constitution_O::stereoisomerIterator si;
	for ( si=(sum)->begin_Stereoisomers(); si!=(sum)->end_Stereoisomers(); si++ )
	{
	    LOG(BF("CandoDatabase_Entities_add(%s) a MonomerName") % (*si)->getName().c_str()  );
	    if ( (*si)->getName() == "" )
	    {
		SIMPLE_ERROR(BF("While adding constitution(%s) a stereoisomer had a blank name!!") % sum->getName() );
	    }
	    this->_Entities.addUnique((*si)->getName(), (*si), _lisp );
	    LOG(BF("CandoDatabase_Entities_add(%s) a PdbName") % (*si)->getPdb().c_str()  );
	    if ( (*si)->getPdb() == "" )
	    {
		SIMPLE_ERROR(BF("While adding constitution(%s) a stereoisomer had a blank pdb name!!") % sum->getName() );
	    }
	    this->_Entities.addUnique((*si)->getPdb(), (*si), _lisp );
	}
    }
#endif


#if 0
/*!
  Add the CouplingRule to the CandoDatabase
*/
    void	CandoDatabase_O::addCouplingRule( Coupling_spRule rule )
    {
	LOG(BF("Adding couplingRule with name=%s to CandoDatabase") % (rule->getName().c_str() ) );
	this->couplingMap.insert(make_pair(rule->getName(),rule));
	this->couplings.push_back(rule);
    }
#endif

/*!
  Add the MonomerFragment to the CandoDatabase
*/
#if 0 //[
    void	CandoDatabase_O::addCoreFragmentCoordinatesTable( FragmentCoordinates_spTable& frag )
    {
	LOG(BF("Adding fragment with name=%s to CandoDatabase") % (frag->getName().c_str() ) );
	frag->setGlobal(true);
	if ( this->_CoreFragmentMap.count(frag->getName()) != 0 ) {
	    SIMPLE_ERROR(BF("The global fragment: ", frag->getName()," has already been defined"));
	}
	this->_CoreFragmentMap[frag->getName()] = frag;
    }

    FragmentCoordinates_spTable	CandoDatabase_O::getCoreFragmentCoordinatesTableWithName(string nm)
    {
	if ( this->_CoreFragmentMap.count(nm) == 0 ) {
	    SIMPLE_ERROR(BF("Could not find global fragment with name: ", nm ));
	}
	return this->_CoreFragmentMap[nm];
    }

/*!
  Add the MonomerFragment to the CandoDatabase
*/
    void	CandoDatabase_O::addFinishFragmentCoordinatesTable( FragmentCoordinates_spTable& frag )
    {
	LOG(BF("Adding fragment with name=%s to CandoDatabase") % (frag->getName().c_str() ) );
	frag->setGlobal(true);
	if ( this->_FinishFragmentMap.count(frag->getName()) != 0 ) {
	    SIMPLE_ERROR(BF("The global fragment: ", frag->getName()," has already been defined"));
	}
	this->_FinishFragmentMap[frag->getName()] = frag;
    }

    FragmentCoordinates_spTable	CandoDatabase_O::getFinishFragmentCoordinatesTableWithName(string nm)
    {
	if ( this->_FinishFragmentMap.count(nm) == 0 ) {
	    SIMPLE_ERROR(BF("Could not find global fragment with name: ", nm ));
	}
	return this->_FinishFragmentMap[nm];
    }
#endif //]


#if 0
    void	CandoDatabase_O::addSystemRepresentedEntityNameSet(RepresentedEntityNameSet_sp group)
    {
	if ( this->recognizesEntityNameSetName(group->getName()) ) {
	    SIMPLE_ERROR(BF("You tried to create a group with a name that is already used: %s") % group->getName() );
	}
	this->_SystemRepresentedEntityNameSets.set(group->getName(),group);
    }

    void	CandoDatabase_O::defineUserMonomerPack(MonomerPack_sp group)
    {
	this->_Entities.set(group->getName(),group);
	this->_UserMonomerPacks->set(group->getName(),group);
	this->_WeakLastUserMonomerPackToChange = group;
    }



    adapt::StringSet_sp	CandoDatabase_O::allUniqueCouplingNames()
    {
	adapt::StringSet_sp			names;
	Constitutions::value_iterator	mi;
	Constitution_O::PlugMap::value_iterator	pi;
	Constitution_sp			cc;
	core::Symbol_sp			couplingName;
	names = adapt::StringSet_O::create();
	names->clear();
	for ( mi=this->begin_Constitutions();
	      mi!=this->end_Constitutions(); mi++ ) {
	    cc = *mi;
            cc->mapPlugs( [&names] (Plug
	    for ( pi=cc->begin_Plugs(); pi!=cc->end_Plugs(); pi++ )
	    {
		couplingName = DirectionalCoupling_O::couplingName((*pi)->getName());
		LOG(BF("Adding couplingName: %s") % couplingName->__repr__() );
		names->insert(couplingName->__repr__());
	    }
	}
	return names;
    }
#endif












// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//
// CoordinateDatabase stuff
//
//
//




CL_LISPIFY_NAME("addMonomerCoordinates");
CL_DEFMETHOD     uint CandoDatabase_O::addMonomerCoordinates(MonomerCoordinates_sp mc)
    {
	SpecificContextSet_sp	allContexts;
	MonomerContext_sp	context;
	SpecificContextSet_O::iterator	si;
	LOG(BF("about to generate specific keys") );
	context = mc->getContext();
//    _lisp->print(BF("%s %u addMonomerCoordinates context:\n%s") % __FILE__% __LINE__ % context->asXmlString().c_str() );
	allContexts = context->getAllSpecificContexts();
	LOG(BF("There are %u specific contexts") % (allContexts->size() ) );
//    _lisp->print(BF("%s %u There are %u specific contexts") % __FILE__% __LINE__ % allContexts->size() );
	uint count = 0;
	for ( si=allContexts->begin(); si!=allContexts->end(); si++ )
	{
//        _lisp->print(BF("    %s %u Looking at specific contexts: %s") % __FILE__% __LINE__ % si->first.c_str() );
	    if ( si->second->allMonomersInDatabase(this->sharedThis<CandoDatabase_O>()) )
	    {
//	    _lisp->print(BF("      %s %u  Adding it") % __FILE__ % __LINE__ );
		LOG(BF("Adding context: %s") % (si->first).c_str()  );
		this->_MonomerCoordinates.set(si->first,mc);
		count++;
	    } else
	    {
//	    _lisp->print(BF("      %s %u  Ignoring it") % __FILE__ % __LINE__ );
	    }
	}
	LOG(BF("Only %u of the specific contexts were used") % count  );
	return count;
    }


CL_LISPIFY_NAME("uniqueMonomerCoordinatesAsCons");
CL_DEFMETHOD     core::List_sp CandoDatabase_O::uniqueMonomerCoordinatesAsCons()
    {
	core::Cons_sp first, cur;
	first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
	cur = first;
	monomerCoordinatesIterator it;
        gctools::SmallOrderedSet<MonomerCoordinates_sp> unique;
	for ( it=this->begin_MonomerCoordinates(); it!=this->end_MonomerCoordinates(); it++ )
	{
	    unique.insert(it->second);
	}
	for ( gctools::SmallOrderedSet<MonomerCoordinates_sp>::iterator ui=unique.begin(); ui!=unique.end(); ui++ )
	{
	    core::Cons_sp one = core::Cons_O::create((*ui),_Nil<core::T_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	return first->cdr();
    }

CL_LISPIFY_NAME("monomerCoordinatesKeysAsCons");
CL_DEFMETHOD     core::List_sp CandoDatabase_O::monomerCoordinatesKeysAsCons()
    {
	core::Cons_sp first, cur;
	first = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::T_O>());
	cur = first;
	monomerCoordinatesIterator it;
	for ( it=this->begin_MonomerCoordinates(); it!=this->end_MonomerCoordinates(); it++ )
	{
          core::Cons_sp one = core::Cons_O::create(it->first,_Nil<core::T_O>());
	    cur->setCdr(one);
	    cur = one;
	}
	return first->cdr();
    }




CL_LISPIFY_NAME("recognizesMonomerCoordinatesKey");
CL_DEFMETHOD     bool	CandoDatabase_O::recognizesMonomerCoordinatesKey(core::Symbol_sp key)
    {
	return this->_MonomerCoordinates.contains(key);
    }

CL_LISPIFY_NAME("getMonomerCoordinatesWithKey");
CL_DEFMETHOD     MonomerCoordinates_sp	CandoDatabase_O::getMonomerCoordinatesWithKey(core::Symbol_sp key )
    {
	if ( this->_MonomerCoordinates.contains(key) )
	{
	    return this->_MonomerCoordinates.get(key);
	}
	return _Nil<MonomerCoordinates_O>();
    }


CL_LISPIFY_NAME("get");
CL_DEFMETHOD     MonomerCoordinates_sp	CandoDatabase_O::get(MonomerContext_sp context)
    {
        core::Symbol_sp key = context->getKey();
	if ( this->_MonomerCoordinates.contains(key) )
	{
	    return this->_MonomerCoordinates.get(key);
	}
	return _Nil<MonomerCoordinates_O>();
    }

CL_LISPIFY_NAME("recognizesContext");
CL_DEFMETHOD     bool	CandoDatabase_O::recognizesContext(MonomerContext_sp context)
    {
	bool	foundIt;
        core::Symbol_sp key = context->getKey();
	LOG(BF("Looking for context: %s") % key.c_str()  );
	foundIt = this->_MonomerCoordinates.contains(key);
	LOG(BF("Found it = %d") % foundIt  );
	return foundIt;
    }

CL_LISPIFY_NAME("saveAs");
CL_DEFMETHOD     void	CandoDatabase_O::saveAs(const string& fn)
    {
	IMPLEMENT_ME();
#ifdef XML_ARCHIVE
	core::XmlSaveArchive_sp	xml;
	xml = core::XmlSaveArchive_O::create();
	xml->put("candoDatabase",this->sharedThis<CandoDatabase_O>());
	xml->saveAs(fn);
#endif
    }
#if 0
    void CandoDatabase_O::removeMonomerCoordinatesNotRequiredByAlchemists(core::List_sp alchemists)
    {
	SpecificContextSet_sp necessary = SpecificContextSet_O::create();
	while ( alchemists.notnilp() )
	{
	    Alchemist_sp alchemist = alchemists->car<Alchemist_O>();
	    SpecificContextSet_sp one = alchemist->allSpecificMonomerContexts();
	    necessary->merge(one);
	    alchemists = alchemists->cdr();
	}
//    SpecificContextSet_sp necessary = scorer->allSpecificMonomerContexts();
	SpecificContextSet_sp available = this->allSpecificMonomerContexts();
	LOG(BF("There are %d necessary MonomerContexts") % necessary->size() );
	LOG(BF("There are %d available MonomerContexts") % available->size() );
	available->remove(necessary);
	LOG(BF("After removing necessary there are %d useless MonomerContexts") % available->size() );
	SpecificContextSet_O::iterator it;
	for ( it=available->begin(); it!=available->end(); it++ )
	{
	    this->_MonomerCoordinates.remove(it->second->asString());
	}
    	//
	// Check for missing contexts
	//
	// Get what is now available and remove it from what is necessary
	// There should be nothing left in necessary once this is done
	//
	available = this->allSpecificMonomerContexts();
	LOG(BF("Now there are %d available MonomerContexts") % available->size() );
	necessary->remove(available);
	LOG(BF("After removing available MonomerContexts from the necessary ones ther are %d missing") % necessary->size() );
	if ( necessary->size() > 0 )
	{
	    //
	    // There are missing contexts, print them and throw an exception
	    //
	    stringstream serr;
	    serr << "There are missing MonomerContexts, the database can not build all of these molecules" << std::endl;
	    for ( SpecificContextSet_O::iterator mt =necessary->begin(); mt!=necessary->end(); mt++ )
	    {
		serr << "  missing: " << (mt->second)->asString() << std::endl;
	    }
	    serr << "Total missing MonomerContexts: " << necessary->size() << std::endl;
	    SIMPLE_ERROR(BF(serr.str()));
	}
    }
#endif

CL_LISPIFY_NAME("allSpecificMonomerContexts");
CL_DEFMETHOD     SpecificContextSet_sp CandoDatabase_O::allSpecificMonomerContexts()
    {
	SpecificContextSet_sp all = SpecificContextSet_O::create();
	monomerCoordinatesIterator mi;
        gctools::SmallOrderedSet<MonomerCoordinates_sp>	uniqueMonomerCoordinates;
	for ( mi=this->_MonomerCoordinates.begin(); mi!=this->_MonomerCoordinates.end(); mi++ )
	{
	    uniqueMonomerCoordinates.insert(mi->second);
	}
	for (gctools::SmallOrderedSet<MonomerCoordinates_sp>::iterator si=uniqueMonomerCoordinates.begin();
	     si!=uniqueMonomerCoordinates.end(); si++ )
	{
	    SpecificContextSet_sp one = (*si)->getContext()->getAllSpecificContexts();
	    all->merge(one);
	}
	return all;
    }

#if INIT_TO_FACTORIES

#define ARGS_CandoDatabase_O_make "(name frame_recognizers represented_entity_name_sets constitutions)"
#define DECL_CandoDatabase_O_make ""
#define DOCS_CandoDatabase_O_make "make AtomReference args: name frame_recognizers represented_entity_name_sets constitutions"
    CandoDatabase_sp CandoDatabase_O::make(core::Symbol_sp name, core::List_sp frameRecognizers, core::List_sp representedEntityNameSets, core::List_sp constitutions)
    {
        GC_ALLOCATE(CandoDatabase_O, me );
	me->_Name = name;
	me->_frameRecognizers.clear();
        for ( auto cur : frameRecognizers ) {
	    FrameRecognizer_sp fr = cur->car<FrameRecognizer_O>();
	    core::Symbol_sp name = fr->getRecognizerName();
	    me->_frameRecognizers.set(name,fr);
	    IMPLEMENT_MEF(BF("Handle setOwner"));
//	    fr->setOwner(this->sharedThis<CandoDatabase_O>());
	}
#if 0
	// System Monomer groups
	core::List_sp monomerPacks = representedEntityNameSets;
	me->_SystemRepresentedEntityNameSets.clear();
        for ( auto cur : monomerPacks ) {
	    RepresentedEntityNameSet_sp group = cur->car<RepresentedEntityNameSet_O>();
	    me->_SystemRepresentedEntityNameSets.set(group->getName(),group);
	}
#endif
	// Constitutions
        for ( auto cur : constitutions ) {
	    Constitution_sp constitution = cur->car<Constitution_O>();
	    me->addEntity(constitution);
	}
	return me;
    };

#else

    core::T_sp	CandoDatabase_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->_Name = translate::from_object<string>::convert(env->lookup(Package(),"name"));
	core::Cons_sp frameRecognizers = translate::from_object<core::Cons_O>::convert(env->lookup(Package(),"frameRecognizers"));
	this->_frameRecognizers.clear();
	for ( auto cur : frameRecognizers ) {
	    FrameRecognizer_sp fr = cur->car<FrameRecognizer_O>();
	    core::Symbol_sp name = fr->getRecognizerName();
	    this->_frameRecognizers.set(name,fr);
	    fr->setOwner(this->sharedThis<CandoDatabase_O>());
	}
#if 0
	// System Monomer groups
	core::List_sp monomerPacks = translate::from_object<core::Cons_O>::convert(env->lookup(Package(),"representedEntityNameSets"));
	this->_SystemRepresentedEntityNameSets.clear();
	for ( auto cur : monomerPacks ) {
	    RepresentedEntityNameSet_sp group = cur->car<RepresentedEntityNameSet_O>();
	    this->_SystemRepresentedEntityNameSets.set(group->getName(),group);
	}
#endif
	// Constitutions
	core::Cons_sp constitutions =  translate::from_object<core::Cons_O>::convert(env->lookup(Pkg(),"constitutions"));
        for ( auto cur : constitutions ) {
	    Constitution_sp constitution = cur->car<Constitution_O>();
	    this->addEntity(constitution);
	}
	return _Nil<core::T_O>();
    }

#endif

























}; //
