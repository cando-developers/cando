/*
    File: monomer.cc
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


//#include	<iostream>
#include <iterator>
#include <vector>
#include <set>
#include <stdio.h>
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/monomer.h>
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/loop.h>
#include <cando/chem/oligomer.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/symbolList.h>
#include <cando/chem/coupling.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/monomerPack.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <clasp/core/numerics.h>
#include <cando/chem/alias.h>
#include <cando/chem/topology.h>
#include <clasp/core/wrappers.h>



namespace chem {

CL_LISPIFY_NAME("CHEM:NAME");
CL_LAMBDA(value isoname);
CL_DOCSTRING(R"doc(isoname-name-set)doc")
CL_DEFUN_SETF core::Symbol_sp chem__isoname_name_set(core::Symbol_sp name, Isoname_sp isoname)
{
  isoname->_Name = name;
  return name;
}

CL_LISPIFY_NAME("CHEM:ISOMER");
CL_LAMBDA(value isoname);
CL_DOCSTRING(R"doc(isoname-isomer-set)doc")
CL_DEFUN_SETF Fixnum chem__isoname_isomer_set(Fixnum isomer, Isoname_sp isoname)
{
  isoname->_Isomer = isomer;
  return isomer;
}


void Isoname_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,name),this->_Name);
  node->field_if_not_default(INTERN_(kw,isomer),this->_Isomer,(Fixnum)0);
  this->Base::fields(node);
}

string Isoname_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " " << _rep_(this->name()) << ">";
  return ss.str();
}


#if 0
Residue_sp Monomer_O::createResidue()
{
  Residue_sp	orig;
  Constitution_sp	con;
  core::Symbol_sp		monName;
  monName = this->getName();
  con = core::eval::funcall(_sym_constitutionForNameOrPdb, getCandoDatabase(), monName );
  Residue_sp res = con->createResidueForStereoisomerName(monName);
  Constitution_sp residueConstitution = res->getConstitution();
  if ( residueConstitution != con )
  {
    SIMPLE_ERROR(BF("Residue created with a different constitution(%s) from the one it was created from constitution(%s)") % _rep_(residueConstitution) % _rep_(con) );
  }
    	//
	// WORKING WORKING WORKING
	//
	// Here we should set the Monomer alias names for the residue
	// and the atom alias names for the atoms of the residue
	// according to the _GroupName
	//
  res->setMonomerAliases(this->_Aliases);
  res->setUniqueLabel(this->_Id);
#if 0
  // Commented out April 2017 - I'm not sure what the point of GroupName is
  if ( getCandoDatabase()->recognizesEntityNameSetName(this->_GroupName) )
  {
    EntityNameSet_sp monomerPack = getCandoDatabase()->getEntity(this->_GroupName).as<EntityNameSet_O>();
    if ( monomerPack->supportsInterestingAtomAliases() )
    {
      core::List_sp atomAliases = monomerPack->getInterestingAtomAliases();
      AtomIndexer_sp aliasAtoms = monomerPack->getAtomIndexerForMonomerName(this->getName());
      ASSERTP(core::cl__length(atomAliases) == aliasAtoms->numberOfAtomNames(),
              "The number of atom names has to match the number of atom aliases");
      AtomIndexer_O::atomNameIterator ii;
      core::List_sp ai;
      for ( ii = aliasAtoms->begin_AtomNames(), ai=atomAliases;
            ii != aliasAtoms->end_AtomNames(); ii++ )
      {
        Atom_sp atom = res->atomWithName(*ii);
        atom->setAlias(oCar(ai));
        ai = oCdr(ai);
      }
#if 0
      {
        stringstream ss;
        ss << "    At this point you should set the aliases of the atoms of this residue"<<std::endl;
        ss << "    Source file: " << __FILE__ << " line number: " << __LINE__ << std::endl;
        ss << "    The _GroupName is ("<<this->_GroupName<<") and the monomer has the aliases("<<this->_Aliases->asString()<< ")"<< std::endl;
        ss << "    The monomer name is ("<<this->getName()<<")"<<std::endl;
        ss << "    The atoms interesting atom names are ("<<monomerPack->getInterestingAtomNamesForMonomerName(this->getName())<<")"<<std::endl;
        ss << "    The atom aliases are ("<<monomerPack->getInterestingAtomAliasesAsString() <<")"<<std::endl;
        _lisp->print(BF("%s"), ss.str());
      }
#endif
    }
  }
#endif
  this->setTempResidue(res);
  return res;
}
#endif


CL_LISPIFY_NAME("getGeneralMonomerContext");
CL_DEFMETHOD     MonomerContext_sp Monomer_O::getGeneralMonomerContext()
{
  MonomerContext_sp	context;
  EntityNameSetWithCap_sp	selfRecognizer, neighborRecognizer;
  CandoDatabase_sp	db;
  Monomer_O::Couplings::iterator	ci;
  Coupling_sp		coupling;
  Monomer_sp		neighborMonomer;
  context = nil<MonomerContext_O>();
  if ( !this->isMonomerContextValid() )
  {
    SIMPLE_ERROR(BF("Monomer context is invalid: %s") % this->sharedThis<Monomer_O>()->description() );
  }
  db = getCandoDatabase();
  context = MonomerContext_O::create();
  selfRecognizer = EntityNameSetWithCap_O::create(); //,db);
//  selfRecognizer->addGroupName(this->getGroupName());
  context->setFocus(selfRecognizer);
  for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ )
  {
    coupling = (ci->second);
	    //
	    // Skip RingClosing couplings, they aren't part of the context
	    //
    if ( coupling.isA<RingCoupling_O>() ) continue;
    neighborRecognizer = EntityNameSetWithCap_O::create(); //,db);
    neighborMonomer = coupling->getOtherSideMonomer(this->sharedThis<Monomer_O>());
//    neighborRecognizer->addGroupName(neighborMonomer->getGroupName());
    context->addNeighbor(ci->first,neighborRecognizer);
  }
  return context;
}

CL_LISPIFY_NAME("plugNamesAndCouplingsAsList");
CL_DEFMETHOD     core::List_sp	Monomer_O::plugNamesAndCouplingsAsList()
{
  core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::Cons_sp cur = first;
  Monomer_O::Couplings::iterator	ci;
  for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ )
  {
    Coupling_sp coupling = (ci->second);
    core::Cons_sp oneOne = core::Cons_O::createList(ci->first,coupling);
    core::Cons_sp one = core::Cons_O::create(oneOne,nil<core::T_O>());
    cur->setCdr(one);
    cur = one;
  }
  return first->cdr();
}


bool Monomer_O::hasRingClosingOutPlug() const
{_OF();
  Couplings::iterator it;
  for ( it=this->_Couplings.begin(); it!=this->_Couplings.end(); it++ )
  {
    if ( it->second->isRingClosing() ) return true;
  }
  return false;
}


#if 0
CL_LISPIFY_NAME("setAliasesFromCons");
CL_DEFMETHOD     void	Monomer_O::setAliasesFromCons(core::List_sp aliases)
{
  this->_Aliases->clear();
  this->_Aliases->insertConsSymbols(aliases);
}
#endif

CL_LISPIFY_NAME("checkForBadConnections");
CL_DEFMETHOD     bool	Monomer_O::checkForBadConnections()
{
  IMPLEMENT_ME();
#if 0
  Couplings::iterator	it;
  Coupling_sp		coup;
  bool			badConnections;
  badConnections = false;
  LOG(BF("Check for bad connections %s") % this->description().c_str() );
  for ( it=this->_Couplings.begin(); it!=this->_Couplings.end(); it++ )
  {
    if ( it->second.nilp() )
    {
      this->addErrorMessage("Coupling with plug("+_rep_(it->first)+") is Null");
      badConnections = true;
      continue;
    }
    coup = it->second;
    if ( !coup->containsMonomer(this->sharedThis<Monomer_O>()) )
    {
      this->addErrorMessage("The coupling through plug("+_rep_(it->first)+") doesn't contain me");
      badConnections = true;
    }
  }
  return badConnections;
#endif
}

CL_LISPIFY_NAME("throwIfBadConnections");
CL_DEFMETHOD     void	Monomer_O::throwIfBadConnections()
{
  IMPLEMENT_ME();
#if 0
  if ( this->checkForBadConnections() )
  {
    LOG(BF("Found bad connections in %s  error: %s") % this->description().c_str() % this->getStatusMessage().c_str()  );
    SIMPLE_ERROR(BF("Bad connections for monomer: %s") % this->sharedThis<Monomer_O>()->description());
  }
#endif
}





CL_LISPIFY_NAME("getSpecificMonomerContext");
CL_DEFMETHOD     MonomerContext_sp Monomer_O::getSpecificMonomerContext()
{
  MonomerContext_sp	context;
  EntityNameSet_sp	selfRecognizer, neighborRecognizer;
  CandoDatabase_sp db;
  Monomer_O::Couplings::iterator	ci;
  Coupling_sp		coupling;
  Monomer_sp		neighborMonomer;
  context = nil<MonomerContext_O>();
  if ( !this->isMonomerContextValid() )
  {
    SIMPLE_ERROR(BF("Monomer context is invalid: %s") %this->sharedThis<Monomer_O>()->description() );
  }
  db = getCandoDatabase();
  context = MonomerContext_O::create();
  selfRecognizer = EntityNameSet_O::create();
  selfRecognizer->setMonomerNameOrPdb(this->currentStereoisomerName());
  context->setFocus(selfRecognizer);
  for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ )
  {
    coupling = (ci->second);
    if ( coupling.isA<DirectionalCoupling_O>() )
    {
      neighborRecognizer = EntityNameSet_O::create();
      neighborMonomer = coupling->getOtherSideMonomer(this->sharedThis<Monomer_O>());
      neighborRecognizer->setMonomerNameOrPdb(neighborMonomer->currentStereoisomerName());
      context->addNeighbor(ci->first,neighborRecognizer);
    } else
    {
      RingCoupling_sp rc = coupling.as<RingCoupling_O>();
      (void)rc;
      LOG(BF("Ignoring out coupling for RingCoupling: %s") % _rep_(rc->getName()) );
    }
  }
  return context;
}




core::T_sp Monomer_O::getMissingRingClosingPlug(Monomer_sp mate)
{
  Constitution_sp	constitution = this->currentTopology()->_Constitution;
  return constitution->getMissingRingClosingPlug(this->sharedThis<Monomer_O>(),mate);
}

uint Monomer_O::numberOfCouplings() const {
  return this->_Couplings.size();
}

CL_DEFMETHOD Topology_sp Monomer_O::currentTopology()
{
  Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[this->_CurrentMonomerIndex]));
  return topology;
}


//
// Destructor
//


CL_LISPIFY_NAME("getInCoupling");
CL_DEFMETHOD     Coupling_sp	Monomer_O::getInCoupling()
{
  Monomer_O::Couplings::iterator	it;
  for ( it=this->_Couplings.begin(); it!=this->_Couplings.end(); it++ ) 
  {
    ASSERTNOTNULL(it->second);
    Coupling_sp coup = it->second;
    if (gc::IsA<DirectionalCoupling_sp>(coup)) {
      DirectionalCoupling_sp dc = gc::As_unsafe<DirectionalCoupling_sp>(coup);
      if ( dc->isInCouplingToMonomer(this->sharedThis<Monomer_O>()) ) return it->second;
    }
  }
  SIMPLE_ERROR(BF("There is no in coupling"));
};


/*! Return the name of the inCoupling.
 * If it doesn't exist return ""
 */
core::Symbol_sp Monomer_O::getInCouplingName()
{
  Monomer_O::Couplings::iterator	it;
  for ( it=this->_Couplings.begin(); it!=this->_Couplings.end(); it++ ) 
  {
    ASSERTNOTNULL(it->second);
    Coupling_sp coup = it->second;
    if (gc::IsA<DirectionalCoupling_sp>(coup)) {
      DirectionalCoupling_sp dc = gc::As_unsafe<DirectionalCoupling_sp>(coup);
      if ( dc->isInCouplingToMonomer(this->sharedThis<Monomer_O>()) ) return it->first;
    }
  }
  return nil<core::Symbol_O>();
};


CL_LISPIFY_NAME("hasInCoupling");
CL_DEFMETHOD     bool	Monomer_O::hasInCoupling()
{
  Monomer_O::Couplings::iterator	it;
  for ( it=this->_Couplings.begin(); it!=this->_Couplings.end(); it++ ) {
    Coupling_sp coup = it->second;
    if (gc::IsA<DirectionalCoupling_sp>(coup)) {
      DirectionalCoupling_sp dc = gc::As_unsafe<DirectionalCoupling_sp>(coup);
      if ( dc->isInCouplingToMonomer(this->sharedThis<Monomer_O>()) ) return true;
    }
  }
  return false;
};

bool	Monomer_O::hasOutCouplings()
{
  Coupling_sp	coup;
  coup = this->getInCoupling();
  return (coup.notnilp());
};



CL_DEFMETHOD Plug_sp	Monomer_O::getPlugNamed(core::Symbol_sp pn)
{
  Topology_sp topology = this->currentTopology();
  if ( !topology->hasPlugNamed(pn) )
  {
    stringstream ss;
    ss << _rep_(topology);
    ss << " does not have plug named(" << _rep_(pn) << ")";
    SIMPLE_ERROR(BF("%s")%ss.str());
  }
  Plug_sp plug = topology->plugNamed(pn);
  return plug;
}


CL_LISPIFY_NAME("getConstitution");
CL_DEFMETHOD     Constitution_sp Monomer_O::getConstitution() {
  Topology_sp topology = this->currentTopology();
  Constitution_sp con = topology->_Constitution;
  return con;
}

Coupling_sp Monomer_O::getCouplingWithPlugName( core::Symbol_sp s) {
  if ( this->_Couplings.contains(s) ) {
    Coupling_sp	coup = this->_Couplings.get(s);
    ASSERTNOTNULL(coup);
    return coup;
  }
  SIMPLE_ERROR(BF("Could not find coupling %s") % _rep_(s));
}

bool	Monomer_O::hasCouplingWithPlugName( core::Symbol_sp s ) {
  return this->_Couplings.contains(s);
}


bool	Monomer_O::hasMatchingPlugNameAndCoupling( core::Symbol_sp plugName, Coupling_sp coup)
{
  Couplings::iterator	it;
  pair<Couplings::iterator,Couplings::iterator>	range;
  int					count;
  Coupling_sp				myCoup;
  bool					foundIt;
  LOG(BF("Looking for plug name: %s and coupling: %s") % _rep_(plugName) % coup->description() );
  count = 0;
  foundIt = false;
  range = this->_Couplings.equal_range(plugName);
  for ( it=range.first; it!=range.second; it++ )
  {
    count++;
    ASSERTNOTNULLP(it->second,"Coupling is undefined");
    myCoup = it->second;
    LOG(BF("Looking at coupling: %s") % myCoup->description().c_str()  );
    if ( myCoup == coup )
    {
      LOG(BF("Found match,setting foundIt to true") );
      foundIt = true;
      break;
    }
  }
  LOG(BF("There were %d couplings with that plug name") % count );
  return foundIt;
}

Monomer_sp Monomer_O::copyAndRemoveCouplings() const {
  auto  newMonomer = gctools::GC<Monomer_O>::copy( *this);
  newMonomer->_Couplings.clear();
  return newMonomer;
}

void Monomer_O::eraseCouplings()
{
  this->_Couplings.clear();
}

CL_LISPIFY_NAME("addCoupling");
CL_DEFMETHOD     void	Monomer_O::addCoupling( core::Symbol_sp plugName,  Coupling_sp coup)
{
  LOG(BF("add coupling to %s") % this->description() );
  LOG(BF("Adding plug(%s) coupling: %s") % _rep_(plugName) % coup->description() );
  this->_Couplings.insert2(plugName, coup->sharedThis<Coupling_O>());
  LOG(BF("After add monomer is %s") % this->description() );
};



void	Monomer_O::resetInCoupling( )
{
  LOG(BF("resetting in coupling of %s") % this->description() );
  Coupling_sp inCoup;
  if ( this->hasInCoupling() ) {
    LOG(BF("Had to remove old coupling") );
	    // Remove the current in coupling
    inCoup = this->getInCoupling();
    this->removeCoupling(inCoup);
  }
}


CL_LISPIFY_NAME("setInCoupling");
CL_DEFMETHOD     void	Monomer_O::setInCoupling( Coupling_sp coup)
{
  core::Symbol_sp inCoupName;
  LOG(BF("setting in coupling of %s") % this->description() );
  ASSERTNOTNULL(coup);
  this->resetInCoupling();
  LOG(BF("Adding in coupling: %s") % coup->description() );
  inCoupName = DirectionalCoupling_O::inPlugName(coup->getName());
  this->addCoupling( inCoupName, coup );
}



CL_LISPIFY_NAME("addOutCoupling");
CL_DEFMETHOD     void	Monomer_O::addOutCoupling( Coupling_sp coup )
{
  ASSERTNOTNULL(coup);
  LOG(BF("adding out coupling to %s") % this->description()  );
  core::Symbol_sp outPlugName = DirectionalCoupling_O::outPlugName(coup->getName());
  LOG(BF("Adding out coupling: %s") % coup->description() );
  this->addCoupling( outPlugName, coup );
}


CL_LISPIFY_NAME("removeCouplingToMonomer");
CL_DEFMETHOD     void	Monomer_O::removeCouplingToMonomer(Monomer_sp mon )
{
  Couplings::iterator		wci;
  bool				foundIt;
  Coupling_sp			coup;
  Monomer_sp			otherMonomer;
  foundIt = false;
  for ( wci=this->_Couplings.begin(); wci!=this->_Couplings.end(); wci++ )
  {
    coup = wci->second;
    otherMonomer = coup->getOtherSideMonomer(this->sharedThis<Monomer_O>());
    if ( otherMonomer == mon )
    {
      foundIt = true;
      break;
    }
  }
  if ( !foundIt )
  {
    SIMPLE_ERROR(BF("Could not find coupling to monomer: %s %s") % mon->description() % this->sharedThis<Monomer_O>()->description() );
  }
  LOG(BF("Removing coupling: %s") % coup->description() );
  this->_Couplings.erase(wci);
};




void	Monomer_O::fixPlugNameForCoupling(Coupling_sp coup)
{
  Couplings::iterator		wci;
  bool				foundIt, isInPlug;
  Coupling_sp			myCoup;
  Monomer_sp			otherMonomer;
  core::Symbol_sp			plugName;
  foundIt = false;
  for ( wci=this->_Couplings.begin(); wci!=this->_Couplings.end(); wci++ )
  {
    myCoup = (wci->second);
    if ( coup == myCoup )
    {
      foundIt = true;
      isInPlug = DirectionalCoupling_O::isInPlugName(wci->first);
      LOG(BF("Original plug name: %s") % _rep_(wci->first) );
      break;
    }
  }
  if ( !foundIt )
  {
    SIMPLE_ERROR(BF("Could not find coupling : %s %s") % coup->description() % this->sharedThis<Monomer_O>()->description() );
  }
  LOG(BF("Fixing the PlugName key for coupling: %s") % coup->description().c_str()  );
  this->_Couplings.erase(wci);
  if ( isInPlug)
  {
    plugName = DirectionalCoupling_O::inPlugName(coup->getName());
  } else
  {
    plugName = DirectionalCoupling_O::outPlugName(coup->getName());
  }
  LOG(BF("Changing plug name to: %s") % _rep_(plugName) );
  this->addCoupling( plugName, coup );
}












CL_LISPIFY_NAME("removeCoupling");
CL_DEFMETHOD     void	Monomer_O::removeCoupling(Coupling_sp coup)
{
  bool				foundIt = false;
  Coupling_sp			myCoup;
  Couplings::iterator		wci;
  for ( wci=this->_Couplings.begin(); wci!=this->_Couplings.end(); wci++ )
  {
    myCoup = wci->second;
    if ( coup == myCoup  )
    {
      foundIt = true;
      break;
    }
  }
  if ( !foundIt )
  {
    SIMPLE_ERROR(BF("Could not find coupling: %s %s") % coup->description() % this->sharedThis<Monomer_O>()->description() );
  }
  LOG(BF("Removing coupling: %s") % myCoup->description().c_str()  );
  this->_Couplings.erase(wci);
}



//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------


string Monomer_O::__repr__() const {
  stringstream ss;
  ss << "#<MONOMER :index " << this->_CurrentMonomerIndex << " :stereoisomer " << this->_CurrentStereoisomerOffset << " :monomer " << _rep_(this->currentStereoisomerName()) << ">";
  return ss.str();
}
void Monomer_O::fields(core::Record_sp node)
{
  node->field_if_not_nil(INTERN_(kw,id),this->_Id);
  node->field(INTERN_(kw,seqnum),this->_SequenceNumber);
  gctools::Vec0<core::T_sp> pairs;
  switch (node->stage()) {
  case core::Record_O::saving:
    {
	    // Accumulate intraresidue bonds into a vector
      for ( Monomer_O::Couplings::iterator ci = this->_Couplings.begin(); ci!=this->_Couplings.end();ci++ ) {
        pairs.push_back(core::Cons_O::create(ci->first,ci->second));
      }
      node->field_if_not_empty(INTERN_(kw,couplings),pairs);
    }
    break;
  case core::Record_O::initializing:
  case core::Record_O::loading:
    {
      this->_Couplings.clear();
      node->field_if_not_empty( INTERN_(kw,couplings), pairs);
      for (size_t i=0; i<pairs.size(); i++ ) {
        core::T_sp pair = pairs[i];
        this->_Couplings.insert2(gc::As<core::Symbol_sp>(oCar(pair)),gc::As<Coupling_sp>(oCdr(pair)));
      }
    }
    break;
  case core::Record_O::patching: {
    // Nothing should need to be done
  }
      break;
  }
  node->field(INTERN_(kw,index),this->_CurrentMonomerIndex);
  node->field(INTERN_(kw,monomers),this->_Monomers);
//  node->field(INTERN_(kw,groupName),this->_GroupName);
  this->Base::fields(node);
}

void Monomer_O::initialize()
{
  this->Base::initialize();
//        this->_Couplings = core::HashTableEq_O::create_default();
  this->_Id = nil<core::Symbol_O>();
  this->_SequenceNumber = 0;
  this->_Monomers.clear();
  this->_CurrentMonomerIndex = -1;
}


CL_DOCSTRING(R"doc(Create a monomer using a list of topology objects)doc")
CL_LAMBDA(&optional topology-list)
CL_DEF_CLASS_METHOD 
Monomer_sp Monomer_O::makeMonomer(core::List_sp topology_list)
{
  auto  me  = gctools::GC<Monomer_O>::allocate_with_default_constructor();
  for ( auto cur : topology_list ) {
    core::Symbol_sp top = gc::As<core::Symbol_sp>(CONS_CAR(cur));
    me->addTopologyName(top);
  }
  return me;
};


CL_DEFMETHOD core::Symbol_sp Monomer_O::currentStereoisomerName() const
{_OF();
  if ( this->_Monomers.size() == 0 ) {
    return nil<core::Symbol_O>();
  }
  Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[this->_CurrentMonomerIndex]));
  core::Symbol_sp name = topology->getStereoisomerName(this->_CurrentStereoisomerOffset);
//  printf("%s:%d In Monomer_O::getName() -> %s name->className() -> %s\n", __FILE__, __LINE__, _rep_(name).c_str(), name->className().c_str());
  return name;
}

#if 0
Fixnum Monomer_O::getIsomer() const
{
  if ( this->_Monomers.size() == 0 ) {
    SIMPLE_ERROR(BF("There are no monomers"));
  }
  return this->_CurrentStereoisomerOffset;
}
#endif

/*! Check for errors and if any of the specific MonomerContexts that
 * this Monomer represents are not found in the CandoDatabase
 */
CL_LISPIFY_NAME("checkForErrorsAndUnknownContexts");
CL_DEFMETHOD     void	Monomer_O::checkForErrorsAndUnknownContexts(CandoDatabase_sp cdb)
{
  IMPLEMENT_ME();
#if 0
  adapt::SymbolSet_sp	allSpecificContextKeys;
  MonomerContext_sp	context;
  bool			allContextsRecognized = true;
  ASSERTNOTNULLP(cdb,"CandoDatabase is undefined");
  this->Monomer_O::checkForErrorsAndUnknownContexts(cdb);
  if ( this->numberOfStereoisomers() == 0 )
  {
    SIMPLE_ERROR(BF("Empty MultiMonomer"));
  }
  if ( !this->isMonomerContextValid() )
  {
    SIMPLE_ERROR(BF("The MonomerContext is invalid"));
  }
  context = this->getGeneralMonomerContext();
  allSpecificContextKeys = context->getAllSpecificKeys();
  stringstream ss;
  ss.str("");
  ss << "There are " << allSpecificContextKeys->size() << " specific contexts" << std::endl;
  this->addStatusMessage(ss.str());
  if ( allSpecificContextKeys->size() == 0 )
  {
    this->addErrorMessage("There are no legal specific contexts");
  } else
  {
    allSpecificContextKeys->map_while_true( [&cdb,this,&allContextsRecognized] (core::Symbol_sp sit) {
        if ( !cdb->recognizesMonomerCoordinatesKey(sit) ) {
          this->addErrorMessage("Monomer->Unrecognized context: "+_rep_(sit));
          allContextsRecognized = false;
          return false;
        }
        return true;
      } );
    if ( allContextsRecognized ) {
      this->addStatusMessage("All contexts are recognized.");
    }
  }
#endif
}

#if 0
CL_DEFMETHOD void	Monomer_O::addIsoname(Isoname_sp name)
{
  Topology_sp topology = name->topology();
  core::List_sp plugs = topology->plugsAsList();
  if (this->_Monomers.size()!=0) {
    Topology_sp topology0 = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[0]));
    if (!topology0->matchesPlugs(topology)) {
      SIMPLE_ERROR(BF("The first topology %s in this monomer does not match all of the plugs in the one you are trying to add %s") % _rep_(topology0) % _rep_(topology));
    }
  }
  this->_Monomers.push_back(name);
}
#endif

CL_DEFMETHOD void Monomer_O::addTopologyName(core::Symbol_sp name)
{
  Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(name));
  core::List_sp plugs = topology->plugsAsList();
  if (this->_Monomers.size()!=0) {
    Topology_sp topology0 = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[0]));
    if (!topology0->matchesPlugs(topology)) {
      SIMPLE_ERROR(BF("The first topology %s in this monomer does not match all of the plugs in the one you are trying to add %s") % _rep_(topology0) % _rep_(topology));
    }
  }
  this->_Monomers.push_back(name);
}


CL_DEFMETHOD core::List_sp Monomer_O::plugNamesAsList() const
{
  if (this->_Monomers.size()==0) {
    return nil<core::T_O>();
  }
  Topology_sp topology0 = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[0]));
  core::List_sp plugs0 = topology0->plugsAsList();
  ql::list result;
  for ( auto ocur : plugs0 ) {
    Plug_sp plug = gc::As<Plug_sp>(CONS_CAR(ocur));
    result << plug->getName();
  }
  return result.cons();
}

CL_DEFMETHOD size_t Monomer_O::numberOfStereoisomers() const
{
  size_t number = 0;
  for ( auto cur : this->_Monomers ) {
    Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(cur));
    number += topology->numberOfStereoisomers();
  }
  return number;
}

CL_DEFMETHOD void Monomer_O::setMonomerIndex(size_t index) {
  size_t offset = index;
  size_t mindex = 0;
  for ( auto cur : this->_Monomers ) {
    Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(cur));
    size_t numberOfStereoisomers = topology->numberOfStereoisomers();
    if (offset<numberOfStereoisomers) {
      this->_CurrentMonomerIndex = mindex;
      this->_CurrentStereoisomerOffset = offset;
      return;
    }
    offset -= numberOfStereoisomers;
    ++mindex;
  }
  SIMPLE_ERROR(BF("The index %lu is beyond the number of stereoisomers in this monomer %lu") % index % this->numberOfStereoisomers());
}

CL_DEFMETHOD size_t Monomer_O::getMonomerIndex() const {
  size_t index = 0;
  for ( size_t idx(0), end(this->_CurrentMonomerIndex); idx < end; ++idx ) {
    Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[idx]));
    index += topology->numberOfStereoisomers();
  }
  index += this->_CurrentStereoisomerOffset;
  return index;
}

  

string Monomer_O::description() const
{//
  stringstream			ss;
  ss.str("");
  ss << "Monomer(";
  ss << "[" <<" id:" << _rep_(this->_Id) << "]=";
  ss << " plugs: ";
  Couplings::iterator	ci;
  for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ ) {
    ss << _rep_(ci->first) << ",";
  }
  ss << ")";
  return ss.str();
}


CL_LISPIFY_NAME("isMonomerContextValid");
CL_DEFMETHOD     bool Monomer_O::isMonomerContextValid()
{
  LOG(BF("status") );
  adapt::SymbolSet_sp keys = adapt::SymbolSet_O::create();
  for ( Couplings::const_iterator it = this->_Couplings.begin(); it!=this->_Couplings.end(); ++it ) {
    keys->insert(it->first);
  }
  LOG(BF("status") );
    	// If a plug name is duplicated then we don't have a valid context and we can't return the context
  if ( keys->size() < this->_Couplings.size() )
  {
    return false;
  }
  return true;
}


#if 0
CL_LISPIFY_NAME("getInterestingAtomIndexer");
CL_DEFMETHOD     AtomIndexer_sp	Monomer_O::getInterestingAtomIndexer()
{
  AtomIndexer_sp		atomIndexer;
  EntityNameSet_sp		entityNameSet;
  CandoDatabase_sp	bdb;
  bdb = getCandoDatabase();
  if ( bdb->recognizesEntityNameSetName(this->getGroupName()) )
  {
    entityNameSet = bdb->getEntityNameSet(this->getGroupName());
    if ( entityNameSet->supportsInterestingAtomAliases() )
    {
      atomIndexer = entityNameSet->getAtomIndexerForMonomerName(this->getStereoisomerName());
      return atomIndexer;
    }
  }
  atomIndexer = AtomIndexer_O::create();
  return atomIndexer;
}
#endif

#if 0
core::List_sp Monomer_O::allAtomAliases()
{
  EntityNameSet_sp		entityNameSet;
  CandoDatabase_sp	bdb;
  bdb = getCandoDatabase();
  ASSERTNOTNULL(bdb);
  ASSERTP(bdb.notnilp(),"Can't get atom aliases because the CandoDatabase is undefined");
  if ( bdb->recognizesEntityNameSetName(this->getGroupName()) )
  {
    LOG(BF("Recognized the EntityNameSetName(%s)")% this->getGroupName() );
    entityNameSet = bdb->getEntityNameSet(this->getGroupName());
    core::List_sp aliases = entityNameSet->getInterestingAtomAliases();
    LOG(BF("Returning with aliases(%s)")% aliases->asString() );
    return aliases;
  } else 
  {
    LOG(BF("Did not recognize the EntityNameSetName(%s)")% this->getGroupName());
  }
  LOG(BF("Returning with nothing"));
  return _Nil<core::T_O>();
}
#endif


class AliasWrapper : public adapt::SymbolSetCartesianProductWrapper
{
public:
  core::T_sp operator()(core::Symbol_sp obj1, core::Symbol_sp obj2) const
  {
    Alias_sp alias = Alias_O::create(obj1,obj2);
    return alias;
  }
};

#if 0
adapt::ObjectSet_sp Monomer_O::getAllAliases()
{
  ASSERTNOTNULL(this->_Aliases);
  if ( this->_Aliases->size() == 0 )
  {
    if ( core::cl__length(this->allAtomAliases()) != 0 )
    {
      SIMPLE_ERROR(BF("The monomer[%s] doesn't have monomer aliases but it has atom aliases(%s) - this should never happen") % this->getStereoisomerName() % _rep_(this->allAtomAliases()) );
    }
    return _Nil<adapt::ObjectSet_O>();
  }
  ASSERT_gt(this->_Aliases->size(),0);
  adapt::SymbolSet_sp atomAliases = adapt::SymbolSet_O::create();
  core::List_sp atomAliasesList = this->allAtomAliases();
  atomAliases->insertConsSymbols(atomAliasesList);
  AliasWrapper wrapper;
  adapt::ObjectSet_sp allAliases = this->_Aliases->cartesianProductWrapped(atomAliases,wrapper);
  return allAliases;
}
#endif



#if 0
bool Monomer_O::recognizesAlias(Alias_sp alias)
{
  AtomIndexer_sp		atomIndexer;
  EntityNameSet_sp		entityNameSet;
  CandoDatabase_sp	bdb;
  LOG(BF("Checking if %s recognizes alias(%s)") % this->description().c_str() % _rep_(alias).c_str()  );
  if ( !this->_Aliases->contains(alias->getMonomerAlias()) )
  {
    LOG(BF("It does not") );
    LOG(BF("The monomer aliases that it recognizes are: %s") % this->_Aliases->asString().c_str()  );
    return false;
  }
  LOG(BF("%s does recognize monomer alias(%s)") % this->description().c_str() % _rep_(alias).c_str() );
  bdb = getCandoDatabase();
  if ( bdb->recognizesEntityNameSetName(this->getGroupName()) )
  {
    LOG(BF("Checking if represents a EntityNameSet with interesting atom alias(%s)") % _rep_(alias->getAtomAlias())  );
    entityNameSet = bdb->getEntityNameSet(this->getGroupName());
    if ( entityNameSet->hasInterestingAtomAlias( alias ) )
    {
      LOG(BF("EntityNameSet does have interesting atom with alias(%s)") % _rep_(alias->getAtomAlias())  );
      return true;
    } else
    {
      LOG(BF("EntityNameSet does not have interesting atom with alias(%s)") % _rep_(alias->getAtomAlias())  );
    }

  }
  LOG(BF("Fail,does not recognize alias") );
  return false;
}
#endif



bool Monomer_O::incrementMonomerIndex()
{
  this->_CurrentStereoisomerOffset++;
  Topology_sp topology = gc::As<Topology_sp>(chem__findTopology(this->_Monomers[this->_CurrentMonomerIndex]));
  if (this->_CurrentStereoisomerOffset < topology->numberOfStereoisomers()) {
    return true;
  }
  this->_CurrentStereoisomerOffset = 0;
  this->_CurrentMonomerIndex++;
  if (this->_CurrentMonomerIndex < this->_Monomers.size()) {
    return true;
  }
  this->_CurrentMonomerIndex = 0;
  return false;
}



gctools::Vec0<Coupling_sp>	Monomer_O::getOutCouplings()
{_OF();
  gctools::Vec0<Coupling_sp>	outCouplings;
  Monomer_O::Couplings::iterator wci;
  Coupling_sp		coup;
  for ( wci=this->_Couplings.begin();
        wci!=this->_Couplings.end(); wci++ ) {
    if ( DirectionalCoupling_O::isOutPlugName(wci->first) )
    {
      ASSERTNOTNULL(wci->second);
      coup = wci->second;
      outCouplings.push_back(coup);
    }
  }
  return outCouplings;
}


Monomer_O::Monomer_O(const Monomer_O& other)
{
  this->_Id = other._Id;
  this->_SequenceNumber = other._SequenceNumber;
  this->_Couplings.clear();
  this->_CurrentMonomerIndex = other._CurrentMonomerIndex;
  this->_CurrentStereoisomerOffset = other._CurrentStereoisomerOffset;
  this->_Monomers.resize(other._Monomers.size());
  for ( size_t ii(0),iEnd(other._Monomers.size()); ii<iEnd; ++ii ) {
    this->_Monomers[ii] = other._Monomers[ii];
  }
}

};
