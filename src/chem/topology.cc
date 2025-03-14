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
#include <clasp/core/hashTableEq.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/atom.h>
#include <cando/chem/molecule.h>
#include <cando/chem/topology.h>
#include <cando/chem/coupling.h>
#include <cando/chem/monomer.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/complexRestraints.h>
#include <clasp/core/wrappers.h>

#ifdef USE_TOPOLOGY
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
CL_DOCSTRING(R"dx(Create a topology from a residue. The constitution may be NIL or a constitution to use to define the atoms)dx");
CL_DEF_CLASS_METHOD Topology_mv Topology_O::makeTopologyFromResidue(chem::Residue_sp residue, core::Symbol_sp topologyName, core::T_sp tconstitution )
{
  core::Symbol_sp stereoisomerName = topologyName;
  Constitution_sp constitution;
  if ( constitution.nilp() ) {
    SIMPLE_ERROR("The constitution was nil - you need to pass a constitution");
  }
  constitution = gc::As<Constitution_sp>(tconstitution);
  ConstitutionAtoms_sp ca = constitution->getConstitutionAtoms();
  Topology_sp topology = Topology_O::make(topologyName, constitution, nil<core::T_O>() );
  core::ComplexVector_T_sp vec = core::ComplexVector_T_O::make(residue->numberOfAtoms(),nil<core::T_O>());
  int atomIndex = 0;
  StereoisomerAtoms_sp stereoisomerAtoms = StereoisomerAtoms_O::create(stereoisomerName);
  for ( auto ai = residue->begin_atoms(); ai!=residue->end_atoms(); ++ai, ++atomIndex) {
    Atom_sp atom = gc::As_unsafe<Atom_sp>(*ai);
    ConstitutionAtomIndex0N constitutionIndex = constitution->getConstitutionAtoms()->index(atom->getName());
    StereoisomerAtom_sp sai = StereoisomerAtom_O::make(atom->getName(),
                                                       atom->getType(atomTypes),
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
  ss << "#<" << this->_instanceClass()->_classNameAsString();
  ss << " :name " << this->_Name->formattedName(true);
#ifdef USE_BOEHM
  ss << " @" << (void*)this;
#endif
  ss << ">";
  return ss.str();
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueForIsomer(size_t isomer) const
{
//  core::clasp_write_string(fmt::format("{}:{} Topology_O::buildResidueForIsomer\n" , __FILE__ , __LINE__ ));
  StereoisomerAtoms_sp info = this->_StereoisomerAtomProperties[isomer];
//  printf("%s:%d buildResidueForIsomer isomer = %lu  stereoisomerAtoms = %s\n", __FILE__, __LINE__, isomer, _rep_(info).c_str());
  LOG("creating residue\n");
  core::Symbol_sp residueName = info->getName();
  if (residueName.unboundp()) {
    SIMPLE_ERROR("residueName for {} was unbound" , _rep_(info));
  }
  Residue_sp res = Residue_O::make(residueName);
  LOG("created residue\n");
//  ASSERT(atomInfo.notnilp());
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
    atom->setProperties(cl__copy_seq(ca->_Properties));
    atom->turnOnFlags(needsBuild);
//    printf("%s:%d  Creating atom@%d -> %s\n", __FILE__, __LINE__, ai->_ConstitutionAtomIndex, _rep_(atom).c_str());
    atoms[ai->_ConstitutionAtomIndex] = atom;
    if (ai->_ConstitutionAtomIndex != idx) {
      SIMPLE_ERROR("The atom {} _ConstitutionAtomIndex {} does not match the StereoisomerAtoms idx {}" , _rep_(atom) , ai->_ConstitutionAtomIndex , idx );
    }
    res->putMatter(idx,atom); // order atoms as in Topology
  }
  for ( size_t i=0, iEnd(constitutionAtoms->numberOfAtoms()); i<iEnd; ++i ) {
    Atom_sp fromAtom = atoms[i];
    ConstitutionAtom_sp ca = (*constitutionAtoms)[i];
//    printf("%s:%d @%zu fromAtom -> %s\n", __FILE__, __LINE__, i, _rep_(fromAtom).c_str());
    for ( auto bi=ca->_Bonds.begin(); bi!=ca->_Bonds.end(); ++bi )
    {
      if ((*bi)->_ToAtomIndex>=atoms.size()) {
        SIMPLE_ERROR("Atom index {} out of bounds (num-atoms {})" , (*bi)->_ToAtomIndex , atoms.size());
      }
      Atom_sp toAtom = atoms[(*bi)->_ToAtomIndex];
//      printf("%s:%d     @%d toAtom -> %s\n", __FILE__, __LINE__, (*bi)->_ToAtomIndex, _rep_(toAtom).c_str());
      if ( !fromAtom->isBondedTo(toAtom) ) {
        BondOrder order = (*bi)->_BondOrder;
        Atom_sp tempFromAtom = fromAtom;
        Atom_sp tempToAtom = toAtom;
        Bond_O::canonicalizeBondOrder(tempFromAtom,tempToAtom,order);
        tempFromAtom->bondTo(tempToAtom,order);
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
    core::T_mv aa_mv = res->atomWithName(name);
    Atom_sp aa = gc::As<Atom_sp>(aa_mv);
    LOG("Setting the configuration of atom({}) to({})" , aa->description().c_str() , _rep_((*sci)->getConfiguration())  ); //
    auto trans = translate::from_object<ConfigurationEnum,std::true_type>((*sci)->getConfiguration())._v;
    aa->setStereochemistryType(chiralCenter);
    aa->setConfiguration(trans);
  }
    //
    //
    // Set dihedral restraints for E/Z pi bonds
    //
    //
  Molecule_sp tempMolecule = Molecule_O::create();
  tempMolecule->addMatter(res);
  core::HashTable_sp cip_priorities = CipPrioritizer_O::assignPrioritiesHashTable(tempMolecule);
  gctools::Vec0<ComplexRestraint_sp>::iterator	tpi;
  for ( tpi=this->_Constitution->_StereoInformation->begin_ComplexRestraints(); 
        tpi!=this->_Constitution->_StereoInformation->end_ComplexRestraints(); tpi++ ) {
    (*tpi)->fillRestraints(res,cip_priorities);
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
  SIMPLE_ERROR("Could not find monomer named {} in topology with names {}" , _rep_(monomerName) , snames.str());
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueForIsoname(Isoname_sp isoname) const {
  if (isoname->_Isomer<this->_StereoisomerAtomProperties.size()) {
      return this->buildResidueForIsomer(isoname->_Isomer);
  }
  stringstream snames;
  for ( size_t isomer(0); isomer<this->_StereoisomerAtomProperties.size(); ++isomer) {
    snames << _rep_(this->_StereoisomerAtomProperties[isomer]->getName()) << " ";
  }
  SIMPLE_ERROR("Could not find isoname {} in topology with names {}" , _rep_(isoname) , snames.str());
}

CL_DEFMETHOD Residue_sp Topology_O::buildResidueSingleName() const {
  core::HashTable_sp cip = core::HashTable_O::createEq(); // dummy hash table
  if (this->_StereoisomerAtomProperties.size()==1) {
      return this->buildResidueForIsomer(0);
  }
  stringstream snames;
  for ( size_t isomer(0); isomer<this->_StereoisomerAtomProperties.size(); ++isomer) {
    snames << _rep_(this->_StereoisomerAtomProperties[isomer]->getName()) << " ";
  }
  SIMPLE_ERROR("This topology does not have a single name - it has names {}" , snames.str());
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
    SIMPLE_ERROR("{}" , ss.str());
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
    SIMPLE_ERROR((se.str()));
  }
}
#endif

struct OverlappingFragments {
  uint	_TimesSeen;
  set<string>	_Fragments;
};

CL_LISPIFY_NAME(connect_residues);
DOCGROUP(cando);
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
    SIMPLE_ERROR("Could not find b0 in plug {} for topology {}" , _rep_(out_plug) , _rep_(prev_topology));
  }
  if ( in_plug_atom_name.nilp() ) {
    SIMPLE_ERROR("Could not find b0 in plug {} for topology {}" , _rep_(in_plug) , _rep_(next_topology));
  }
//  printf("%s:%d  out_plug_atom_name = %s  in_plug_atom_name = %s\n", __FILE__, __LINE__, _rep_(out_plug_atom_name).c_str(), _rep_(in_plug_atom_name).c_str());
  Atom_sp out_atom = gc::As_unsafe<Atom_sp>(prev_residue->atomWithName(out_plug_atom_name));
  Atom_sp in_atom = gc::As_unsafe<Atom_sp>(next_residue->atomWithName(in_plug_atom_name));
//  printf("%s:%d  out_atom = %s  in_atom = %s\n", __FILE__, __LINE__, _rep_(out_atom).c_str(), _rep_(in_atom).c_str());
  BondOrder boIn0 = in_plug->getBondOrder0();
  BondOrder boOut0 = out_plug->getBondOrder0();
  Bond_O::canonicalizeBondOrder(in_atom,out_atom,boIn0);
  if (boIn0!=boOut0) {
    SIMPLE_ERROR("bond order of in-plug {}/{} does not match bond order of out-plug {}/{}" , _rep_(in_plug_name) , bondOrderToString(boIn0) , _rep_(out_plug_name) , bondOrderToString(boOut0));
  }
  in_atom->bondTo(out_atom, boIn0);
  if (in_plug->getB1().notnilp()) {
    Atom_sp out_atom = gc::As_unsafe<Atom_sp>(prev_residue->atomWithName(out_plug->getB1()));
    Atom_sp in_atom = gc::As_unsafe<Atom_sp>(next_residue->atomWithName(in_plug->getB1()));
    BondOrder boIn1 = in_plug->getBondOrder1();
    BondOrder boOut1 = out_plug->getBondOrder1();
    Bond_O::canonicalizeBondOrder(in_atom,out_atom,boIn1);
    if (boIn1!=boOut1) {
      SIMPLE_ERROR("bond order 1 of in-plug {}/{} does not match bond order of out-plug {}/{}" , _rep_(in_plug_name) , bondOrderToString(boIn1) , _rep_(out_plug_name) , bondOrderToString(boOut1));
    }
    in_atom->bondTo(out_atom, boIn1);
  }
}

CL_DEFMETHOD StereoisomerAtoms_sp Topology_O::getStereoisomerAtoms(core::Symbol_sp stereoisomerName) const
{
  for (int i(0); i<this->_StereoisomerAtomProperties.size(); ++i ) {
    if (stereoisomerName == this->_StereoisomerAtomProperties[i]->getName() ) {
      return this->_StereoisomerAtomProperties[i];
    }
  }
  SIMPLE_ERROR("Could not find StereoisomerAtoms named {}" , _rep_(stereoisomerName));
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
  return nil<plugOType>();
}


CL_NAME(CHEM:TOPOLOGY/PLUGS-AS-LIST);
CL_DEFMETHOD     core::List_sp Topology_O::plugsAsList()
{
  core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::Cons_sp cur = first;
  LOG("The number of plugs = {}" , this->_Plugs.size()  );
  for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
  {
    LOG("Adding plug: {}" , i->second->getName());
    core::Cons_sp one = core::Cons_O::create(i->second,nil<core::T_O>());
    cur->setCdr(one);
    cur = one;
  }
  LOG("Returning the plugs");
  return first->cdr();
}


CL_LISPIFY_NAME("plugsWithMatesAsList");
CL_DEFMETHOD     core::List_sp Topology_O::plugsWithMatesAsList()
{
  core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::Cons_sp cur = first;
  LOG("The number of plugs = {}" , this->_Plugs.size()  );
  for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
  {
	    // skip origin plugs
    if ( !i->second.isA<PlugWithMates_O>() ) continue;
    LOG("Adding plug: {}" , _rep_(i->second) );
    core::Cons_sp one = core::Cons_O::create(i->second,nil<core::T_O>());
    cur->setCdr(one);
    cur = one;
  }
  LOG("Returning the plugs" );
  return first->cdr();
}

CL_LISPIFY_NAME("outPlugsAsList");
CL_DEFMETHOD     core::List_sp Topology_O::outPlugsAsList()
{
  core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(),nil<core::T_O>());
  core::Cons_sp cur = first;
  for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
  {
    Plug_sp plug = i->second;
    if ( !plug->getIsIn() )
    {
      core::Cons_sp one = core::Cons_O::create(i->second,nil<core::T_O>());
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
    ss << fmt::format("{}@{} " , _rep_(i->second->getName()) , (void*)i->second->getName().raw_() );
  }
  ss << ">";
  return ss.str();
}


bool Topology_O::hasMatchingPlugsWithMates(adapt::SymbolSet_sp plugSet)
{
  adapt::SymbolSet_sp myPlugSet = adapt::SymbolSet_O::create();
  for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
  {
    LOG("Looking at plug[{}]" , _rep_(i->second->getName()) );
    if (i->second.isA<PlugWithMates_O>() )
    {
      myPlugSet->insert(i->first);
    }
  }
  LOG("Topology plugs are: {}" , myPlugSet->asString() );
  bool match = plugSet->equal(myPlugSet);
  LOG("Do they match the plugs passed as arguments[{}] --> {}" , plugSet->asString() , match );
  return match;
}

CL_DOCSTRING(R"dx(Return true if this topology has identically named plugs of the other)dx");
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
  LOG("Checking if monomer[{}] matches the topology environment" , mon->description() );
  uint numPlugsWithMates = 0;
  LOG("{}:{} monomer->{}  number of plugs: {}\n" , __FILE__ , __LINE__ , _rep_(mon) , this->_Plugs.size());
  for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++)
  {
    LOG("Looking at plug[{}]" , _rep_(i->second->getName()) );
    if (!i->second.isA<PlugWithMates_O>() )
    {
      LOG("It's not a PlugWithMates");
      continue;
    }
    if ( !mon->hasCouplingWithPlugName(i->second->getName()) ) 
    {
      SIMPLE_WARN("The monomer {} doesn't have a coupling with plug name[{}]" ,_rep_(mon) ,_rep_(i->second->getName()) );
      return false;
    }
    numPlugsWithMates++;
  }
  if ( numPlugsWithMates != mon->numberOfCouplings() ) 
  {
    SIMPLE_WARN("There is a mismatch with the number of plugs {} in the topology {} and the number of couplings {} for the monomer {}" ,numPlugsWithMates ,_rep_(this->asSmartPtr()) ,mon->numberOfCouplings() ,_rep_(mon));
    return false;
  }
  LOG("They match");
  return true;
}



core::T_sp Topology_O::provideMissingRingClosingPlug( Monomer_sp mon )
{
  core::T_sp missingRingClosingPlug = nil<core::T_O>();
  uint numPlugsWithMates = 0;
  for ( Plugs::iterator i=this->_Plugs.begin();
        i!= this->_Plugs.end(); i++)
  {
    if (!i->second.isA<PlugWithMates_O>() ) continue;
    if ( gc::IsA<OutPlug_sp>(i->second) && gc::As_unsafe<OutPlug_sp>(i->second)->getIsRingClosing())
    {
      missingRingClosingPlug = gc::As<OutPlug_sp>(i->second);
      continue;
    }
    if ( !mon->hasCouplingWithPlugName(i->second->getName()) ) return nil<core::T_O>();
    numPlugsWithMates++;
  }
  if ( numPlugsWithMates != mon->numberOfCouplings() ) return nil<core::T_O>();
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
  LOG("Result = {}" , res );
  return res;
};

CL_LISPIFY_NAME("plugNamed");
CL_DEFMETHOD     Plug_sp Topology_O::plugNamed(core::Symbol_sp name)
{
  bool res = this->_Plugs.contains(name);
  if ( !res ) {
    SIMPLE_ERROR("Could not find plug with name {}" , _rep_(name));
  }
  return this->_Plugs.get(name);
};

CL_DEFMETHOD core::List_sp Topology_O::allOutPlugNamesThatMatchInPlugName(core::Symbol_sp inPlugName)
{
  ql::list result;
  for ( Plugs::iterator i=this->_Plugs.begin(); i!= this->_Plugs.end(); i++) {
    Plug_sp plug = i->second;
    if ( plug->getIsOut() ) {
      if (DirectionalCoupling_O::outPlugNameMatchesInPlugName(plug->getName(),inPlugName)) {
        result << plug->getName();
      }
    }
  }
  return result.cons();
}

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
CL_DOCSTRING(R"dx(Set the property **symbol** of **this** (a chem:topology) to **value**.)dx");
CL_DEFMETHOD void Topology_O::setProperty(core::Symbol_sp symbol, core::T_sp value)
{
  this->_Properties = core::core__put_f(this->_Properties,value,symbol);
}

CL_DOCSTRING(R"dx(Return the property **symbol** of **this** (a chem:topology) - if it isn't defined return NIL.)dx");
CL_NAME(CHEM:TOPOLOGY/GET-PROPERTY);
CL_DEFMETHOD core::T_sp Topology_O::getProperty(core::Symbol_sp symbol)
{
  core::T_sp res = core::cl__getf(this->_Properties,symbol,unbound<core::T_O>());
  if (res.unboundp()) {
    stringstream props;
    props << _rep_(this->_Properties);
    SIMPLE_ERROR("You asked for an unknown property[{}] for topology[{}@{}] - the available properties are[{}]" , _rep_(symbol) , this->__repr__() , (void*)this , props.str()  );
  }
  return res;
}

CL_DOCSTRING(R"dx(Return the property **symbol** of **this** (a chem:topology) - if it isn't defined return **defval**.)dx");
CL_LISPIFY_NAME("getPropertyOrDefault");
CL_DEFMETHOD core::T_sp Topology_O::getPropertyOrDefault(core::Symbol_sp prop,core::T_sp defval)
{
  return core::cl__getf(this->_Properties,prop,defval);
}

CL_DOCSTRING(R"dx(Return T if the property **symbol** of **this** (a chem:topology) is defined.)dx");
CL_LISPIFY_NAME("hasProperty");
CL_DEFMETHOD bool Topology_O::hasProperty(core::Symbol_sp symbol)
{
  return !core::cl__getf(this->_Properties,symbol,unbound<core::T_O>()).unboundp();
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
CL_DOCSTRING(R"dx(Return (values stereoisomer-atoms coded-or-absolute isomer))dx");CL_DEFMETHOD core::T_mv Topology_O::lookupStereoisomerAtoms(Fixnum isomer)
{
  if (this->_StereoisomerCoding == kw::_sym_coded) {
    return Values(this->_StereoisomerAtomProperties[isomer],this->_StereoisomerCoding,core::make_fixnum(isomer));
  }
  if (this->_StereoisomerAtomProperties.size() < 1 ) {
    return Values(this->_StereoisomerAtomProperties[0],this->_StereoisomerCoding,core::make_fixnum(isomer));
  }
  return Values(nil<core::T_O>(),this->_StereoisomerCoding,core::make_fixnum(isomer));
}

CL_DEFMETHOD void Topology_O::walkStereoisomerAtoms(core::Function_sp func)
{
  for ( auto entry : this->_StereoisomerAtomProperties ) {
    core::eval::funcall(func,entry);
  }
};

};

#endif
