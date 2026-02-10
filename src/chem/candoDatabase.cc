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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE

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
#include <clasp/core/array.h>
#include <clasp/core/package.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/evaluator.h>
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

namespace chem {
/*!Store a list of all CandoDatabases so that they can be consolidated
 */
//    core::WeakVector0<CandoDatabase_O>	AllCandoDatabases;

void CandoDatabase_O::initialize() {
  
  this->Base::initialize();
  this->_Topologys = core::HashTable_O::createEq();
  this->_Topologys->setupThreadSafeHashTable();
  //    AllCandoDatabases.append(this->sharedThis<CandoDatabase_O>());
}

void CandoDatabase_O::fields(core::Record_sp node) {
  node->field_if_not_nil(INTERN_(kw, Name), this->_Name);
  node->field(INTERN_(kw, Topologys), this->_Topologys);
}


#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("constitutionsAsList");
CL_DEFMETHOD core::List_sp CandoDatabase_O::constitutionsAsList() {
  
  return this->entitiesSubClassOfAsList(gc::As<core::Instance_sp>(core::cl__find_class(Constitution_O::static_classSymbol())));
};
#endif

#ifdef USE_TOPOLOGY
bool CandoDatabase_O::recognizesMonomerName(core::Symbol_sp name) {
  
  return this->recognizesEntitySubClassOf(name,
                                          gc::As<core::Instance_sp>(core::cl__find_class(Stereoisomer_O::static_classSymbol())));
}


CL_LISPIFY_NAME("addEntity");
CL_DEFMETHOD Entity_sp CandoDatabase_O::addEntity(Entity_sp entity) {
  IMPLEMENT_ME();
#if 0
  this->_Entities.set(entity->getName(),entity);
  if ( entity.isA<Constitution_O>() )
  {
    entity.as<Constitution_O>()->addStereoisomersToCandoDatabase(this->sharedThis<CandoDatabase_O>());
  }
  return entity;
#endif
}
#endif




#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("representedEntityNameSetsAsList");
CL_DEFMETHOD core::List_sp CandoDatabase_O::representedEntityNameSetsAsList() {
  
  return this->entitiesSubClassOfAsList(
      gc::As<core::Instance_sp>(core::cl__find_class(RepresentedEntityNameSet_O::static_classSymbol())));
}

CL_LISPIFY_NAME("recognizesRepresentedEntityNameSet");
CL_DEFMETHOD bool CandoDatabase_O::recognizesRepresentedEntityNameSet(core::Symbol_sp name) {
  
  return this->recognizesEntityOfClass(
      name, gc::As_unsafe<core::Instance_sp>(core::cl__find_class(RepresentedEntityNameSet_O::static_classSymbol())));
}

bool CandoDatabase_O::recognizesMonomerPack(core::Symbol_sp name) {
  
  auto found = this->_Entities->find(name);
  return found.has_value() && gc::As<Entity_sp>(*found).isA<MonomerPack_O>();
}
#endif

SYMBOL_EXPORT_SC_(ChemKwPkg, constitution);
SYMBOL_EXPORT_SC_(ChemKwPkg, entity);
SYMBOL_EXPORT_SC_(ChemKwPkg, frameRecognizer);

#ifdef USE_TOPOLOGY
core::T_sp CandoDatabase_O::oGetResource(core::Symbol_sp resource, core::Symbol_sp name) {
  
  if (resource == chemkw::_sym_constitution) {
    return core::eval::funcall(_sym_constitutionForNameOrPdb, this->asSmartPtr(), name);
  }
  if (resource == chemkw::_sym_entity) {
    return this->getEntity(name);
  }
  if (resource == chemkw::_sym_frameRecognizer) {
    return this->getFrameRecognizer(name);
  }
  SIMPLE_ERROR("Unknown CandoDatabase Resource {} {}", _rep_(resource), _rep_(name));
}
#endif

#ifdef USE_TOPOLOGY
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
  SIMPLE_ERROR(("Unknown ObjRef selector: "+ref->getSelector() ));
}
#endif

#ifdef USE_TOPOLOGY
bool CandoDatabase_O::recognizesEntityNameSetName(core::Symbol_sp nm) {
  return this->recognizesEntityOfClass(
      nm, gc::As_unsafe<core::Instance_sp>(core::cl__find_class(EntityNameSet_O::static_classSymbol())));
}

CL_LISPIFY_NAME("getEntityNameSet");
CL_DEFMETHOD EntityNameSet_sp CandoDatabase_O::getEntityNameSet(core::Symbol_sp nm) {
  return gc::As<EntityNameSet_sp>(
      this->getEntityOfClass(nm, gc::As_unsafe<core::Instance_sp>(core::cl__find_class(EntityNameSet_O::static_classSymbol()))));
};

CL_LISPIFY_NAME("getRepresentedEntityNameSet");
CL_DEFMETHOD RepresentedEntityNameSet_sp CandoDatabase_O::getRepresentedEntityNameSet(core::Symbol_sp nm) {
  return this
      ->getEntityOfClass(nm,
                         gc::As_unsafe<core::Instance_sp>(core::cl__find_class(RepresentedEntityNameSet_O::static_classSymbol())))
      .as<RepresentedEntityNameSet_O>();
};

CL_LISPIFY_NAME("recognizesSetOrConstitutionOrMonomerName");
CL_DEFMETHOD bool CandoDatabase_O::recognizesSetOrConstitutionOrMonomerName(core::Symbol_sp nm) {
  
  return this->_Entities->contains(nm);
}
#endif


#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("getEntity");
CL_DEFMETHOD Entity_sp CandoDatabase_O::getEntity(core::Symbol_sp nm) const {
  return gc::As<Entity_sp>(this->_Entities->gethash(nm));
}

Entity_sp CandoDatabase_O::getEntityOfClass(core::Symbol_sp nm, core::Instance_sp mc) {
  Entity_sp entity = this->getEntity(nm);
  if (core::cl__class_of(entity) != mc) {
    SIMPLE_ERROR("You asked for Entity[{}] of class[{}] and there is one with that name but it has the wrong class[{}]",
                 _rep_(nm), mc->className(), entity->className());
  }
  return entity;
}

CL_LISPIFY_NAME("expandEntityNameToTerminals");
CL_DEFMETHOD adapt::SymbolSet_sp CandoDatabase_O::expandEntityNameToTerminals(core::Symbol_sp nm) const {
  
  Entity_sp entity = this->getEntity(nm);
  return entity->expandedNameSet();
}

CL_LISPIFY_NAME("expandEntityNamesToTerminals");
CL_DEFMETHOD adapt::SymbolSet_sp CandoDatabase_O::expandEntityNamesToTerminals(adapt::SymbolSet_sp nm) const {
  
  adapt::SymbolSet_sp names = adapt::SymbolSet_O::create();
  nm->map([this, &names](core::Symbol_sp sym) {
    Entity_sp entity = this->getEntity(sym);
    names->insertSymbolSet(entity->expandedNameSet());
  });
  return names;
}
#endif

#if 0
bool	CandoDatabase_O::monomersAreDefinedForSetOrConstitutionOrMonomerName(core::Symbol_sp nm)
{
  LOG("looking for({})" , _rep_(nm)  );
  EntityNameSet_sp group;
  if ( this->recognizesEntityNameSetName(nm) )
  {
    LOG("Found group" );
    group = this->getEntityNameSet(nm);
    return group->isFinalized();
  }
    	// For constitutions and monomer names their monomer names are always defined
  return true;
}
#endif

#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("getMonomersForSetOrConstitutionOrMonomerName");
CL_DEFMETHOD adapt::SymbolSet_sp CandoDatabase_O::getMonomersForSetOrConstitutionOrMonomerName(core::Symbol_sp nm) {
  
  // Use expandEntityNameToTerminals
  DEPRECATED();
#if 0
  adapt::SymbolSet_sp	names;
  Entity_sp obj;

#ifdef DEBUG_ON

  if ( !this->_Entities.contains(nm) )
  {
    LOG("The namespace contains {} names: " , this->_Entities.size()  );
//        for ( adapt::SymbolMap<CandoDatabase_ODependent>::iterator it=this->_Entities.begin();
    for ( adapt::SymbolMap<Entity_O>::iterator it=this->_Entities.begin();
          it!=this->_Entities.end();
          it++ )
    {
      LOG("Namespace contains name: {}" , _rep_(it->first) );
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
#endif

#if 0
RepresentativeList_sp CandoDatabase_O::expandEntityNameToListOfRepresentatives(core::Symbol_sp nm)
{
  Entity_sp obj;
  geom::ObjectList_sp 	ll;
  LOG("Looking for name({})" , _rep_(nm)  );
  ASSERT(this->recognizesEntityName(nm));
  obj = this->getEntity(nm);
  ll = obj->expandedRepresentativeList();
  return ll;
}
#endif

#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("recognizesNameOrPdb");
CL_DEFMETHOD bool CandoDatabase_O::recognizesNameOrPdb(core::Symbol_sp nm) {
  if (!this->recognizesMonomerName(nm))
    return false;
  //    ASSERTNOTNULL(this->_Names);
  //    ASSERTvalue(this->_Entities.contains(nm), nm); // REMOVE

  Entity_sp obj = this->getEntity(nm);
  if (obj->isTerminalName())
    return true;
  LOG("The name({}) is in the namespace but its not a terminal", nm);
  return false;
};
#endif

#ifdef XML_ARCHIVE
void CandoDatabase_O::archiveBase(core::ArchiveP node) {
  this->Base::archiveBase(node);
  node->attribute("name", this->_Name);
  node->attributeIfNotDefault<string>("dateCreated", this->_DateCreated, "");
  node->attributeIfNotDefault<string>("dateUpdated", this->_DateUpdated, "");
  node->archiveSymbolMap("entities", this->_Entities);
  node->archiveSymbolMap("frameRecognizers", this->_frameRecognizers);
  //	node->archiveMapOfObjectsSubClassOf( this->_MonomerCoordinates );
  node->needsFinalization();
}
#endif

#if 0
CL_LISPIFY_NAME("constitutionForNameOrPdb");
CL_DEFMETHOD     Constitution_sp	CandoDatabase_O::constitutionForNameOrPdb(core::Symbol_sp name)
{
  Entity_sp obj;
  ASSERTP(this->_Entities.contains(name),"Could not find: "+_rep_(name)+" in Environment");
  obj = this->_Entities.get(name);
  ASSERTP(obj->hasConstitution(), "The name: "+_rep_(name)+" is not part of a Constitution" );
  return obj->constitution();
}
#endif

#if 0
CL_LISPIFY_NAME("constitutionNameForNameOrPdb");
CL_DEFMETHOD     core::Symbol_sp CandoDatabase_O::constitutionNameForNameOrPdb(core::Symbol_sp name)
{
  Constitution_sp con;
  con = this->constitutionForNameOrPdb(name);
  return con->getName();
}
#endif
#if 0
core::Symbol_sp CandoDatabase_O::getMonomerNameForNameOrPdb(core::Symbol_sp name)
{
  Entity_sp dep;
  ASSERTP(this->_Entities.contains(name), "Could not find name("+_rep_(name)+") in Environment");
  dep = this->_Entities.get(name);
  ASSERTP(dep->isTerminalName(), "Asked for MonomerName of non-terminal name("+_rep_(name)+")");
  return dep->getName();
}

core::Symbol_sp CandoDatabase_O::getPdbNameForNameOrPdb(core::Symbol_sp name)
{
  Entity_sp 			dep;
  Stereoisomer_sp			si;
  ASSERTP(this->_Entities.contains(name), "Could not find name("+_rep_(name)+") in Environment");
  dep = this->_Entities.get(name);
  ASSERTP(dep->isTerminalName(), "Asked for MonomerName of non-terminal name("+_rep_(name)+")");
  si = (dep).as<Stereoisomer_O>();
  return si->getPdb();
}
#endif

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
/*!
  Add the CouplingRule to the CandoDatabase
*/
void	CandoDatabase_O::addCouplingRule( Coupling_spRule rule )
{
  LOG("Adding couplingRule with name={} to CandoDatabase" , (rule->getName().c_str() ) );
  this->couplingMap.insert(make_pair(rule->getName(),rule));
  this->couplings.push_back(rule);
}
#endif

/*!
  Add the MonomerFragment to the CandoDatabase
*/
#if 0  //[
void	CandoDatabase_O::addCoreFragmentCoordinatesTable( FragmentCoordinates_spTable& frag )
{
  LOG("Adding fragment with name={} to CandoDatabase" , (frag->getName().c_str() ) );
  frag->setGlobal(true);
  if ( this->_CoreFragmentMap.count(frag->getName()) != 0 ) {
    SIMPLE_ERROR(("The global fragment: ", frag->getName()," has already been defined"));
  }
  this->_CoreFragmentMap[frag->getName()] = frag;
}

FragmentCoordinates_spTable	CandoDatabase_O::getCoreFragmentCoordinatesTableWithName(string nm)
{
  if ( this->_CoreFragmentMap.count(nm) == 0 ) {
    SIMPLE_ERROR(("Could not find global fragment with name: ", nm ));
  }
  return this->_CoreFragmentMap[nm];
}

/*!
  Add the MonomerFragment to the CandoDatabase
*/
void	CandoDatabase_O::addFinishFragmentCoordinatesTable( FragmentCoordinates_spTable& frag )
{
  LOG("Adding fragment with name={} to CandoDatabase" , (frag->getName().c_str() ) );
  frag->setGlobal(true);
  if ( this->_FinishFragmentMap.count(frag->getName()) != 0 ) {
    SIMPLE_ERROR(("The global fragment: ", frag->getName()," has already been defined"));
  }
  this->_FinishFragmentMap[frag->getName()] = frag;
}

FragmentCoordinates_spTable	CandoDatabase_O::getFinishFragmentCoordinatesTableWithName(string nm)
{
  if ( this->_FinishFragmentMap.count(nm) == 0 ) {
    SIMPLE_ERROR(("Could not find global fragment with name: ", nm ));
  }
  return this->_FinishFragmentMap[nm];
}
#endif //]

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
//
// CoordinateDatabase stuff
//
//
//

#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("addMonomerCoordinates");
CL_DEFMETHOD uint CandoDatabase_O::addMonomerCoordinates(MonomerCoordinates_sp mc) {
  SpecificContextSet_sp allContexts;
  MonomerContext_sp context;
  SpecificContextSet_O::iterator si;
  LOG("about to generate specific keys");
  context = mc->getContext();
  //  core::clasp_write_string(fmt::format("{} {} addMonomerCoordinates context:\n{}\n" , __FILE__, __LINE__ ,
  //  context->asXmlString() ));
  allContexts = context->getAllSpecificContexts();
  LOG("There are {} specific contexts", (allContexts->size()));
  //    _lisp->print(BF("%s %u There are %u specific contexts") , __FILE__% __LINE__ , allContexts->size() );
  uint count = 0;
  for (si = allContexts->begin(); si != allContexts->end(); si++) {
    //        _lisp->print(BF("    %s %u Looking at specific contexts: %s") , __FILE__% __LINE__ , si->first.c_str() );
    if (si->second->allMonomersInDatabase(this->sharedThis<CandoDatabase_O>())) {
      //	    _lisp->print(BF("      %s %u  Adding it") , __FILE__ , __LINE__ );
      LOG("Adding context: {}", (si->first).c_str());
      this->_MonomerCoordinates.set(si->first, mc);
      count++;
    } else {
      //	    core::clasp_write_string(fmt::format("      {} {}  Ignoring it\n" , __FILE__ , __LINE__ ));
    }
  }
  LOG("Only {} of the specific contexts were used", count);
  return count;
}

CL_LISPIFY_NAME("uniqueMonomerCoordinatesAsList");
CL_DEFMETHOD core::List_sp CandoDatabase_O::uniqueMonomerCoordinatesAsList() {
  core::Cons_sp first, cur;
  first = core::Cons_O::create(nil<core::T_O>(), nil<core::T_O>());
  cur = first;
  monomerCoordinatesIterator it;
  gctools::SmallOrderedSet<MonomerCoordinates_sp> unique;
  for (it = this->begin_MonomerCoordinates(); it != this->end_MonomerCoordinates(); it++) {
    unique.insert(it->second);
  }
  for (gctools::SmallOrderedSet<MonomerCoordinates_sp>::iterator ui = unique.begin(); ui != unique.end(); ui++) {
    core::Cons_sp one = core::Cons_O::create((*ui), nil<core::T_O>());
    cur->setCdr(one);
    cur = one;
  }
  return first->cdr();
}

CL_LISPIFY_NAME("monomerCoordinatesKeysAsList");
CL_DEFMETHOD core::List_sp CandoDatabase_O::monomerCoordinatesKeysAsList() {
  core::Cons_sp first, cur;
  first = core::Cons_O::create(nil<core::T_O>(), nil<core::T_O>());
  cur = first;
  monomerCoordinatesIterator it;
  for (it = this->begin_MonomerCoordinates(); it != this->end_MonomerCoordinates(); it++) {
    core::Cons_sp one = core::Cons_O::create(it->first, nil<core::T_O>());
    cur->setCdr(one);
    cur = one;
  }
  return first->cdr();
}

CL_LISPIFY_NAME("recognizesMonomerCoordinatesKey");
CL_DEFMETHOD bool CandoDatabase_O::recognizesMonomerCoordinatesKey(core::Symbol_sp key) {
  return this->_MonomerCoordinates.contains(key);
}

CL_LISPIFY_NAME("getMonomerCoordinatesWithKey");
CL_DEFMETHOD MonomerCoordinates_sp CandoDatabase_O::getMonomerCoordinatesWithKey(core::Symbol_sp key) {
  if (this->_MonomerCoordinates.contains(key)) {
    return this->_MonomerCoordinates.get(key);
  }
  return nil<MonomerCoordinates_O>();
}

CL_LISPIFY_NAME("CandoDatabase-get");
CL_DEFMETHOD MonomerCoordinates_sp CandoDatabase_O::get(MonomerContext_sp context) {
  core::Symbol_sp key = context->getKey();
  if (this->_MonomerCoordinates.contains(key)) {
    return this->_MonomerCoordinates.get(key);
  }
  return nil<MonomerCoordinates_O>();
}

CL_LISPIFY_NAME("recognizesContext");
CL_DEFMETHOD bool CandoDatabase_O::recognizesContext(MonomerContext_sp context) {
  bool foundIt;
  core::Symbol_sp key = context->getKey();
  LOG("Looking for context: {}", key.c_str());
  foundIt = this->_MonomerCoordinates.contains(key);
  LOG("Found it = {}", foundIt);
  return foundIt;
}

CL_LISPIFY_NAME("saveAs");
CL_DEFMETHOD void CandoDatabase_O::saveAs(const string &fn) {
  IMPLEMENT_ME();
#ifdef XML_ARCHIVE
  core::XmlSaveArchive_sp xml;
  xml = core::XmlSaveArchive_O::create();
  xml->put("candoDatabase", this->sharedThis<CandoDatabase_O>());
  xml->saveAs(fn);
#endif
}
#endif

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
  LOG("There are {} necessary MonomerContexts" , necessary->size() );
  LOG("There are {} available MonomerContexts" , available->size() );
  available->remove(necessary);
  LOG("After removing necessary there are {} useless MonomerContexts" , available->size() );
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
  LOG("Now there are {} available MonomerContexts" , available->size() );
  necessary->remove(available);
  LOG("After removing available MonomerContexts from the necessary ones ther are {} missing" , necessary->size() );
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
    SIMPLE_ERROR((serr.str()));
  }
}
#endif

#ifdef USE_TOPOLOGY
CL_LISPIFY_NAME("allSpecificMonomerContexts");
CL_DEFMETHOD SpecificContextSet_sp CandoDatabase_O::allSpecificMonomerContexts() {
  SpecificContextSet_sp all = SpecificContextSet_O::create();
  monomerCoordinatesIterator mi;
  gctools::SmallOrderedSet<MonomerCoordinates_sp> uniqueMonomerCoordinates;
  for (mi = this->_MonomerCoordinates.begin(); mi != this->_MonomerCoordinates.end(); mi++) {
    uniqueMonomerCoordinates.insert(mi->second);
  }
  for (gctools::SmallOrderedSet<MonomerCoordinates_sp>::iterator si = uniqueMonomerCoordinates.begin();
       si != uniqueMonomerCoordinates.end(); si++) {
    SpecificContextSet_sp one = (*si)->getContext()->getAllSpecificContexts();
    all->merge(one);
  }
  return all;
}
#endif


CL_LISPIFY_NAME(make-cando-database);
CL_DEF_CLASS_METHOD
CandoDatabase_sp CandoDatabase_O::make(core::Symbol_sp name) {
  auto me = gctools::GC<CandoDatabase_O>::allocate();
  me->_Name = name;
#if 0
  me->_frameRecognizers.clear();
  for ( auto cur : frameRecognizers ) {
    FrameRecognizer_sp fr = cur->car<FrameRecognizer_O>();
    core::Symbol_sp name = fr->getRecognizerName();
    me->_frameRecognizers.set(name,fr);
//	    fr->setOwner(this->sharedThis<CandoDatabase_O>());
  }
	// System Monomer groups
  core::List_sp monomerPacks = representedEntityNameSets;
  me->_SystemRepresentedEntityNameSets.clear();
  for ( auto cur : monomerPacks ) {
    RepresentedEntityNameSet_sp group = cur->car<RepresentedEntityNameSet_O>();
    me->_SystemRepresentedEntityNameSets.set(group->getName(),group);
  }
	// Constitutions
  for ( auto cur : constitutions ) {
    Constitution_sp constitution = cur->car<Constitution_O>();
    me->addEntity(constitution);
  }
#endif
  return me;
};

CL_DEFMETHOD void CandoDatabase_O::walk_topologys(core::Function_sp func) {
  this->_Topologys->maphash([&func](core::T_sp key, core::T_sp value) { core::eval::funcall(func, value); });
}

CL_DEFMETHOD void CandoDatabase_O::setf_findTopology(core::T_sp name, core::T_sp topology) {
  this->_Topologys->setf_gethash(name, topology);
}

core::T_mv CandoDatabase_O::findTopology(core::T_sp name, bool errorp) const {
  core::T_mv result_mv = this->_Topologys->gethash(name);
  core::MultipleValues &values = core::lisp_multipleValues();
  if (errorp && values.second(result_mv.number_of_values()).nilp()) {
    SIMPLE_ERROR("Could not find topology with name {}", core::_rep_(name));
  }
  return result_mv;
}

DOCGROUP(cando);
CL_LAMBDA(name &optional errorp);
CL_DEFUN core::T_mv chem__findTopology(core::T_sp name, bool errorp) {
  CandoDatabase_sp cdb = chem::getCandoDatabase();
  return cdb->findTopology(name, errorp);
};

}; // namespace chem
