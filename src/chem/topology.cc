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
//#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/common.h>
#include <clasp/core/bformat.h>
#include <clasp/core/lispStream.h>
#include <clasp/core/ql.h>
#include <clasp/core/evaluator.h>
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
#include <cando/chem/stereochemistry.h>
#include <cando/chem/complexRestraints.h>
#include <clasp/core/wrappers.h>

namespace chem
{


void Topology_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,name),this->_Name);
  node->field(INTERN_(kw,constitution),this->_Constitution);
  node->field(INTERN_(kw,stereoisomer_coding), this->_StereoisomerCoding);
  node->field(INTERN_(kw,stereo_atom_props),this->_StereoisomerAtomProperties);
  node->field(INTERN_(kw,properties),this->_Properties);
  node->field(INTERN_(kw,defaultoutplug),this->_DefaultOutPlugName);
  node->field(INTERN_(kw,plugs),this->_Plugs);
  // this->Base::fields(node); Base is CxxObject_O
}


CL_LISPIFY_NAME(makeTopologyFromResidue);
CL_LAMBDA(residue topology-name &optional constitution);
CL_DOCSTRING("Create a topology from a residue. The constitution may be NIL or a constitution to use to define the atoms");
CL_DEF_CLASS_METHOD Topology_mv Topology_O::makeTopologyFromResidue(chem::Residue_sp residue, core::Symbol_sp topologyName, core::T_sp tconstitution )
{
  core::Symbol_sp stereoisomerName = topologyName;
  Constitution_sp constitution;
  if ( constitution.nilp() ) {
    SIMPLE_ERROR(BF("The constitution was nil - you need to pass a constitution"));
  }
  constitution = gc::As<Constitution_sp>(tconstitution);
  ConstitutionAtoms_sp ca = constitution->getConstitutionAtoms();
  Topology_sp topology = Topology_O::make(topologyName, constitution, _Nil<core::T_O>() );
  core::VectorTNs_sp vec = core::VectorTNs_O::make(residue->numberOfAtoms(),_Nil<core::T_O>());
  int atomIndex = 0;
  StereoisomerAtoms_sp stereoisomerAtoms = StereoisomerAtoms_O::create(stereoisomerName);
  for ( auto ai = residue->begin_atoms(); ai!=residue->end_atoms(); ++ai, ++atomIndex) {
    Atom_sp atom = *ai;
    ConstitutionAtomIndex0N constitutionIndex = constitution->getConstitutionAtoms()->index(atom->getName());
    StereoisomerAtom_sp sai = StereoisomerAtom_O::make(atom->getName(),
                                                       atom->getType(),
                                                       atom->getCharge(),
                                                       constitutionIndex);
    stereoisomerAtoms->addStereoisomerAtom(sai);
  }
  topology->_StereoisomerAtomProperties.push_back(stereoisomerAtoms);
  topology->_StereoisomerCoding = kw::_sym_coded;
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

string Topology_O::__repr__() const {
  stringstream ss;
  ss << "#<TOPOLOGY ";
  ss << " :name " << _rep_(this->_Name);
#ifdef USE_BOEHM
  ss << " @" << (void*)this;
#endif
  ss << ">";
  return ss.str();
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueForIsomer(size_t isomer) const
{
//  core::write_bf_stream(BF("%s:%d Topology_O::buildResidueForIsomer\n") % __FILE__ % __LINE__ );
  StereoisomerAtoms_sp info = this->_StereoisomerAtomProperties[isomer];
//  printf("%s:%d buildResidueForIsomer isomer = %lu  stereoisomerAtoms = %s\n", __FILE__, __LINE__, isomer, _rep_(info).c_str());
  LOG(BF("creating residue\n"));
  core::Symbol_sp residueName = info->getName();
  if (residueName.unboundp()) {
    SIMPLE_ERROR(BF("residueName for %s was unbound") % _rep_(info));
  }
  Residue_sp res = Residue_O::make(residueName);
  LOG(BF("created residue\n"));
  ASSERT(atomInfo.notnilp());
  ConstitutionAtoms_sp constitutionAtoms = this->_Constitution->getConstitutionAtoms();
  size_t numAtoms = constitutionAtoms->numberOfAtoms();
  gctools::Vec0<Atom_sp> atoms;
  atoms.resize(numAtoms);
  res->resizeContents(numAtoms);
  size_t idx = 0;
  for ( size_t idx=0, idxEnd(numAtoms); idx<idxEnd; ++idx ) {
    StereoisomerAtom_sp ai = (*info)[idx];
    Atom_sp atom = Atom_O::create();
    ConstitutionAtom_sp ca = (*constitutionAtoms)[ai->constitutionAtomIndex()];
    atom->setName(ca->_AtomName);
    atom->setElement(ca->_Element);
    atom->setType(ai->_AtomType);
    atom->setCharge(ai->_AtomCharge);
    atom->turnOnFlags(needsBuild);
//    printf("%s:%d  Creating atom@%d -> %s\n", __FILE__, __LINE__, ai->_ConstitutionAtomIndex, _rep_(atom).c_str());
    atoms[ai->_ConstitutionAtomIndex] = atom;
    res->putMatter(idx,atom); // order atoms as in Topology
  }
  for ( size_t i=0, iEnd(constitutionAtoms->numberOfAtoms()); i<iEnd; ++i ) {
    Atom_sp fromAtom = atoms[i];
    ConstitutionAtom_sp ca = (*constitutionAtoms)[i];
//    printf("%s:%d @%zu fromAtom -> %s\n", __FILE__, __LINE__, i, _rep_(fromAtom).c_str());
    for ( auto bi=ca->_Bonds.begin(); bi!=ca->_Bonds.end(); ++bi )
    {
      if ((*bi)->_ToAtomIndex>=atoms.size()) {
        SIMPLE_ERROR(BF("Atom index %d out of bounds (num-atoms %d)") % (*bi)->_ToAtomIndex % atoms.size());
      }
      Atom_sp toAtom = atoms[(*bi)->_ToAtomIndex];
//      printf("%s:%d     @%d toAtom -> %s\n", __FILE__, __LINE__, (*bi)->_ToAtomIndex, _rep_(toAtom).c_str());
      if ( fromAtom->atLowerUniqueAtomOrderThan(toAtom) ) {
        BondOrder order = (*bi)->_BondOrder;
        fromAtom->bondTo(toAtom,order);
      }
    }
  }
  // Now add stereochemical restraints

    //
    // Set chiral restraints
    //
  Stereoisomer_sp si = this->_Constitution->getStereoisomerWithName(residueName);
  gctools::Vec0<StereoConfiguration_sp>::iterator	sci;
  for (sci=si->_Configurations_begin();sci!=si->_Configurations_end();sci++){
    core::T_sp name = (*sci)->getAtomName();
    Atom_sp aa = res->atomWithName(name);
    LOG(BF("Setting the configuration of atom(%s) to(%s)") % aa->description().c_str() % _rep_((*sci)->getConfiguration())  ); //
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
  for ( tpi=this->_Constitution->_StereoInformation->begin_ComplexRestraints(); 
        tpi!=this->_Constitution->_StereoInformation->end_ComplexRestraints(); tpi++ ) {
    (*tpi)->fillRestraints(res);
  }
  return res;
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueForMonomerName(core::Symbol_sp monomerName) const {
  for ( size_t isomer(0); isomer<this->_StereoisomerAtomProperties.size(); ++isomer) {
    if (this->_StereoisomerAtomProperties[isomer]->getName() == monomerName) {
      return this->buildResidueForIsomer(isomer);
    }
  }
  stringstream snames;
  for ( size_t isomer(0); isomer<this->_StereoisomerAtomProperties.size(); ++isomer) {
    snames << _rep_(this->_StereoisomerAtomProperties[isomer]->getName()) << " ";
  }
  SIMPLE_ERROR(BF("Could not find monomer named %s in topology with names %s") % _rep_(monomerName) % snames.str());
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueForIsoname(Isoname_sp isoname) const {
  if (isoname->_Isomer<this->_StereoisomerAtomProperties.size()) {
    return this->buildResidueForIsomer(isoname->_Isomer);
  }
  stringstream snames;
  for ( size_t isomer(0); isomer<this->_StereoisomerAtomProperties.size(); ++isomer) {
    snames << _rep_(this->_StereoisomerAtomProperties[isomer]->getName()) << " ";
  }
  SIMPLE_ERROR(BF("Could not find isoname %s in topology with names %s") % _rep_(isoname) % snames.str());
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueSingleName() const {
  if (this->_StereoisomerAtomProperties.size()==1) {
    return this->buildResidueForIsomer(0);
  }
  stringstream snames;
  for ( size_t isomer(0); isomer<this->_StereoisomerAtomProperties.size(); ++isomer) {
    snames << _rep_(this->_StereoisomerAtomProperties[isomer]->getName()) << " ";
  }
  SIMPLE_ERROR(BF("This topology does not have a single name - it has names %s") % snames.str());
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

StereoisomerAtoms_sp Topology_O::getStereoisomerAtoms(core::Symbol_sp stereoisomerName) const
{
  for (int i(0); i<this->_StereoisomerAtomProperties.size(); ++i ) {
    if (stereoisomerName == this->_StereoisomerAtomProperties[i]->getName() ) {
      return this->_StereoisomerAtomProperties[i];
    }
  }
  SIMPLE_ERROR(BF("Could not find StereoisomerAtoms named %s") % _rep_(stereoisomerName));
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

CL_DOCSTRING("Return true if this topology has identically named plugs of the other");
CL_DEFMETHOD bool Topology_O::matchesPlugs(Topology_sp other)
{
  if (this->_Plugs.size() != other->_Plugs.size()) {
    return false;
  } else {
    for ( auto ocur : this->_Plugs ) {
      Plug_sp oplug = ocur.second;
      if (!other->_Plugs.contains(oplug->getName())) return false;
    }
  }
  return true;
}

bool	Topology_O::matchesMonomerEnvironment( Monomer_sp mon )
{
  LOG(BF("Checking if monomer[%s] matches the topology environment") % mon->description() );
  uint numPlugsWithMates = 0;
  LOG(BF("%s:%d monomer->%s  number of plugs: %d\n") % __FILE__ % __LINE__ % _rep_(mon) % this->_Plugs.size());
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
      SIMPLE_WARN(BF("The monomer %s doesn't have a coupling with plug name[%s]") % _rep_(mon) % _rep_(i->second->getName()) );
      return false;
    }
    numPlugsWithMates++;
  }
  if ( numPlugsWithMates != mon->numberOfCouplings() ) 
  {
    SIMPLE_WARN(BF("There is a mismatch with the number of plugs %d in the topology %s and the number of couplings %d for the monomer %s") % numPlugsWithMates % _rep_(this->asSmartPtr()) % mon->numberOfCouplings() % _rep_(mon));
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


CL_DEFMETHOD core::Symbol_sp Topology_O::getStereoisomerName(size_t idx) const {
  return this->_StereoisomerAtomProperties[idx]->_StereoisomerName;
}




void	Topology_O::initialize()
{
  this->Base::initialize();
}


SYMBOL_EXPORT_SC_(KeywordPkg,coded);
SYMBOL_EXPORT_SC_(KeywordPkg,absolute);


CL_LISPIFY_NAME("clearProperty");
CL_DEFMETHOD void Topology_O::clearProperty(core::Symbol_sp symbol)
{
  this->_Properties = core::core__rem_f(this->_Properties,symbol);
}
CL_DOCSTRING("Set the property **symbol** of **this** (a chem:topology) to **value**.");
CL_DEFMETHOD void Topology_O::setProperty(core::Symbol_sp symbol, core::T_sp value)
{
  this->_Properties = core::core__put_f(this->_Properties,value,symbol);
}

CL_DOCSTRING("Return the property **symbol** of **this** (a chem:topology) - if it isn't defined return NIL.");
CL_LISPIFY_NAME("getProperty");
CL_DEFMETHOD core::T_sp Topology_O::getProperty(core::Symbol_sp symbol)
{
  core::T_sp res = core::cl__getf(this->_Properties,symbol,_Unbound<core::T_O>());
  if (res.unboundp()) {
    stringstream props;
    props << _rep_(this->_Properties);
    SIMPLE_ERROR(BF("You asked for an unknown property[%s] for topology[%s@%p] - the available properties are[%s]") % _rep_(symbol) % this->__repr__() % this % props.str()  );
  }
  return res;
}

CL_DOCSTRING("Return the property **symbol** of **this** (a chem:topology) - if it isn't defined return **defval**.");
CL_LISPIFY_NAME("getPropertyOrDefault");
CL_DEFMETHOD core::T_sp Topology_O::getPropertyOrDefault(core::Symbol_sp prop,core::T_sp defval)
{
  return core::cl__getf(this->_Properties,prop,defval);
}

CL_DOCSTRING("Return T if the property **symbol** of **this** (a chem:topology) is defined.");
CL_LISPIFY_NAME("hasProperty");
CL_DEFMETHOD bool Topology_O::hasProperty(core::Symbol_sp symbol)
{
  return !core::cl__getf(this->_Properties,symbol,_Unbound<core::T_O>()).unboundp();
}



CL_DOCSTRING(R"(Coding can be either :absolute or :coded. :absolute means the bitvector/integer stereoisomer index is used to determine which stereoisomer is specified.  Coded means the stereoisomer index is looked up in the _StereoisomerAtomProperties.)");
CL_DEFMETHOD void Topology_O::setStereoisomerAtoms(core::Symbol_sp coding, core::List_sp stereoisomer_atoms) {
  this->_StereoisomerCoding = coding;
  this->_StereoisomerAtomProperties.clear();
  for ( auto tentry : stereoisomer_atoms ) {
    StereoisomerAtoms_sp entry = gc::As<StereoisomerAtoms_sp>(CONS_CAR(tentry));
    this->_StereoisomerAtomProperties.push_back(entry);
  }
}



CL_LISPIFY_NAME("lookupStereoisomerAtoms");
CL_DOCSTRING("Return (values stereoisomer-atoms coded-or-absolute isomer)")
CL_DEFMETHOD core::T_mv Topology_O::lookupStereoisomerAtoms(Fixnum isomer)
{
  if (this->_StereoisomerCoding == kw::_sym_coded) {
    return Values(this->_StereoisomerAtomProperties[isomer],this->_StereoisomerCoding,core::make_fixnum(isomer));
  }
  if (this->_StereoisomerAtomProperties.size() < 1 ) {
    return Values(this->_StereoisomerAtomProperties[0],this->_StereoisomerCoding,core::make_fixnum(isomer));
  }
  return Values(_Nil<core::T_O>(),this->_StereoisomerCoding,core::make_fixnum(isomer));
}


};
