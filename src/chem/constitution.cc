/*
    File: constitution.cc
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
#include <cando/adapt/adapters.h>
#include <clasp/core/symbol.h>
#include <cando/adapt/stringList.h>
#include <clasp/core/evaluator.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/array.h>
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


#ifdef USE_TOPOLOGY
namespace chem {


Constitution_sp Constitution_O::make(core::Symbol_sp name, core::String_sp comment, ConstitutionAtoms_sp constitutionAtoms, StereoInformation_sp stereoInformation, core::List_sp plugs, core::List_sp topologies)
  {
    auto  me  = gctools::GC<Constitution_O>::allocate_with_default_constructor();
      me->_Name = name;
      me->_Comment = comment;
      me->_ConstitutionAtoms = constitutionAtoms;
      me->_StereoInformation = stereoInformation;
      stereoInformation->validate();
      if ( me->_StereoInformation.notnilp() ) {
        gctools::As<StereoInformation_sp>(me->_StereoInformation)->validate();
      }
      {
	  me->_PlugsByName.clear();
          for ( auto cur : plugs ) {
            Plug_sp p = oCar(cur).as<Plug_O>();
	      
	      if ( me->_PlugsByName.count(p->getName())>0 )
	      {
                SIMPLE_ERROR("There is already a plug named: {}" , _rep_(p->getName()) );
	      }
	      core::Symbol_sp key = p->getName();
	      ASSERTF(!key->isKeywordSymbol(),("Don't use keyword symbols for plug names"));
	      me->_PlugsByName.set(key,p);
	  }
      }
      {
	  me->_Topologies.clear();
          for ( auto cur : topologies ) {
            Topology_sp t = oCar(cur).as<Topology_O>();
	      
	      if ( me->_Topologies.count(t->getName())>0 )
	      {
                SIMPLE_ERROR(("There is already a topology named: "+ _rep_(t->getName()) ));
	      }
	      me->_Topologies.set(t->getName(),t);
	  }
      }
      return me;
  };

void Constitution_O::addStereoisomersToCandoDatabase(CandoDatabase_sp db)
{
	for ( stereoisomerIterator si=this->begin_Stereoisomers(); 
	      si!=this->end_Stereoisomers(); si++ )
	{
	    db->addEntity(*si);
	}
    }


void Constitution_O::add_topology(Topology_sp topology) {
  // Add the named topology to this constitution
  this->_Topologies.set(topology->getName(),topology);
}
  


/*!
	Return a copy of the residue that this constitution defines
*/
void	Constitution_O::makeResidueConsistentWithStereoisomerNamed(Residue_sp res,
								   core::Symbol_sp stereoisomerName,
                                                                   core::HashTable_sp cip)
{
  DEPRECATED();
  Stereoisomer_sp				si;
  string					atomName;
  Atom_sp					aa;
  core::T_sp bdb = getCandoDatabase();
  Constitution_sp residueConstitution = this->asSmartPtr();
  core::Symbol_sp fullName = gc::As<core::Symbol_sp>(core::eval::funcall(_sym_monomerNameForNameOrPdb,bdb,stereoisomerName));
  core::Symbol_sp pdbName = gc::As<core::Symbol_sp>(core::eval::funcall(_sym_pdbNameForNameOrPdb,bdb,stereoisomerName));
  res->setName(fullName);
  res->setPdbName(pdbName);
    //
    // Set chiral restraints
    //
  si = this->_StereoInformation->getStereoisomer(stereoisomerName);
  gctools::Vec0<StereoConfiguration_sp>::iterator	sci;
  for (sci=si->_Configurations_begin();sci!=si->_Configurations_end();sci++){
    aa = gc::As_unsafe<Atom_sp>(res->atomWithName((*sci)->getAtomName()));
    LOG("Setting the configuration of atom({}) to({})" , aa->description().c_str() , _rep_((*sci)->getConfiguration())  ); //
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
      (*tpi)->fillRestraints(res,cip);
  }
}

    bool	Constitution_O::recognizesStereoisomerName(core::Symbol_sp nm)
{ 
    return this->_StereoInformation->recognizesNameOrPdb(nm);
};



RepresentativeList_sp	Constitution_O::expandedRepresentativeList() const
{
    gctools::Vec0<Stereoisomer_sp>::const_iterator	si;
    RepresentativeList_sp allRepresentatives  = RepresentativeList_O::create();
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++)
    {
	RepresentativeList_sp oneList = (*si)->expandedRepresentativeList();
	allRepresentatives->vectorPushExtend(oneList);
    }
    return allRepresentatives;
}



core::T_sp Constitution_O::getMissingRingClosingPlug(Monomer_sp mon, Monomer_sp mate)
{
  core::T_sp plug;
  core::T_sp missing = nil<core::T_O>();
  TopologyMap::iterator	ti;
  gctools::SmallOrderedSet<Topology_sp>	candidateTopologies;
  for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ ) 
  {
    plug = (ti->second)->provideMissingRingClosingPlug(mon);
    if ( plug.notnilp() ) 
    {
      if ( gc::As<OutPlug_sp>(plug)->recognizesRingClosingMate(mate->currentStereoisomerName()) )
      {
        candidateTopologies.insert(ti->second);
        missing = plug;
      }
    }
  }
  if ( candidateTopologies.size() > 1 )
  {
    stringstream ss;
    ss << "In Constitution(" << this->constitutionName() << ")" << std::endl;
    ss << "there are more than one topologies that match the current monomer environment: " << mon->description() << std::endl;
    ss << "that are missing a ring closing plug"<<std::endl;
    SIMPLE_ERROR("{}" , ss.str() );
  }
  if ( candidateTopologies.size() == 0 )
  {
    SIMPLE_ERROR("There are no Topologies with missing ring closing plugs");
  }
  return missing;
}

    core::Symbol_sp Constitution_O::pdbFromNameOrPdb(core::Symbol_sp nm)
    {
	return this->_StereoInformation->pdbFromNameOrPdb(nm);
    };


CL_LISPIFY_NAME("getMonomerNameAsStringSet");
CL_DEFMETHOD adapt::SymbolSet_sp	Constitution_O::getMonomerNamesAsSymbolSet()
{
  IMPLEMENT_ME();
#if 0
    return this->_StereoInformation->getMonomerNamesAsSymbolSet();
#endif
};


bool Constitution_O::hasStereoisomerWithName(core::Symbol_sp stereoisomerName)
{
    const_stereoisomerIterator	si;
    LOG("Looking for Stereoisomer with name({})" , _rep_(stereoisomerName) );
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++)
    {
      LOG("Looking at stereoisomer name({})" , _rep_((*si)->getName()) );
	if ( (*si)->getName() == stereoisomerName )
	{
	    LOG("Found it!!!!" ); //
	    return true;
	}
    }
    LOG("Could not find it" ); //
    return false;
}

    Stereoisomer_sp Constitution_O::getStereoisomerWithName(core::Symbol_sp stereoisomerName) const
{
    const_stereoisomerIterator	si;
    LOG("Looking for Stereoisomer with name({})" , _rep_(stereoisomerName) );
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++)
    {
      LOG("Looking at stereoisomer name({})" , _rep_((*si)->getName()) );
	if ( (*si)->getName() == stereoisomerName )
	{
	    LOG("Found it!!!!" ); //
	    return *si;
	}
    }
    LOG("Could not find it" ); //
    return nil<Stereoisomer_O>();
}



CL_DEFMETHOD core::List_sp Constitution_O::stereoisomersAsList() { return this->_StereoInformation->stereoisomersAsList(); };
CL_DEFMETHOD core::List_sp Constitution_O::topologiesAsList() {
    core::List_sp result = nil<core::T_O>();
    for ( TopologyMap::iterator it = this->_Topologies.begin(); it!=this->_Topologies.end(); ++it ) {
        result = core::Cons_O::create(it->second,result);
    }
    return result;
}

core::List_sp Constitution_O::plugsAsList() {
    core::List_sp result = nil<core::T_O>();
    for ( PlugMap::iterator it = this->_PlugsByName.begin(); it!=this->_PlugsByName.end(); ++it ) {
        result = core::Cons_O::create(it->second,result);
    }
    return result;
}

CL_DEFMETHOD core::List_sp Constitution_O::plugsWithMatesAsList()
{
    core::List_sp first = nil<core::T_O>();
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


#if 0
// This is all done in oligomer.cc now
/*!
	Return a copy of the residue that this constitution defines
*/
CL_LISPIFY_NAME("createResidueForStereoisomerName");
CL_DEFMETHOD     Residue_sp	Constitution_O::createResidueForStereoisomerName(core::Symbol_sp stereoisomerName)
{
  Residue_sp res = this->_ConstitutionAtoms->buildResidue();
  this->makeResidueConsistentWithStereoisomerNamed(res,stereoisomerName);
  return res;
}
#endif


core::Symbol_sp	Constitution_O::nameFromNameOrPdb(core::Symbol_sp nm)
{
  return this->_StereoInformation->nameFromNameOrPdb(nm);
};



string	Constitution_O::description() const
{
stringstream	ss;
    ss << this->className() << "[";
    ss << "Name(";
    ss << _rep_(this->_Name);
    ss <<") ";
    ss << "]";
    return ss.str();
}

core::List_sp Constitution_O::getPlugNames()
{
  ql::list result;
  PlugMap::iterator pi;
  for ( pi=this->_PlugsByName.begin(); pi!=this->_PlugsByName.end(); pi++ )
  {
    result << pi->first;
  }
  return result.result();
}


CL_LISPIFY_NAME("simplestTopologyWithPlugNamed");
CL_DEFMETHOD     Topology_sp	Constitution_O::simplestTopologyWithPlugNamed(core::Symbol_sp name)
{
    TopologyMap::iterator	ti;
Topology_sp			tres;
int				connects,temp;
tres = nil<Topology_O>();
    connects = 9999999;
    for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ )
    {
      LOG("Looking for plug[{}]@{} in topology[{}]" , _rep_(name) , name.get() , ti->second->description());
	if ( ti->second->hasPlugNamed(name) )
	{
          LOG("Found the plug[{}]" , _rep_(name) );
	    temp = ti->second->numberOfPlugs();
	    if ( temp < connects )
	    {
		LOG("Found one with only {} plugs" , temp );
		connects = temp;
		tres = ti->second;
	    } else
	    {
              LOG("Did not find the plug[{}]" , _rep_(name) );
	    }
	}
    }
#ifdef	DEBUG_ON
    if ( tres.notnilp() )
    {
      LOG("simplest topology is named: {}" , _rep_(tres->getName())  );
    } else
    {
      LOG("there is no topology that has the plug[{}]" , _rep_(name) );
    }
#endif
    return tres;
}


CL_LISPIFY_NAME("topologyWithName");
CL_DEFMETHOD     Topology_sp	Constitution_O::topologyWithName(core::Symbol_sp name) const
    {
	Topology_sp			tres;
	ASSERTF(this->_Topologies.contains(name),("There is no topology with name[{}] in constitution[{}]") , _rep_(name) , _rep_(this->getName()) );
	return this->_Topologies.get(name);
    }



string Constitution_O::__repr__() const {
  stringstream ss;
  ss << "#<CONSTITUTION ";
  ss << " :name " << _rep_(this->_Name);
#ifdef USE_BOEHM
  ss << " @" << (void*)this;
#endif
  ss << ">";
  return ss.str();
}



Topology_sp	Constitution_O::getTopologyForMonomerEnvironment(Monomer_sp mon )
{
  LOG("constitution: {}  monomer: {}\n" , this->__repr__() , _rep_(mon));
  TopologyMap::iterator	ti;
  Topology_sp tres;
  for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ ) 
  {
    LOG("Checking topology({})"% _rep_(ti->first) );
    if ( ti->second->matchesMonomerEnvironment(mon) ) {
      LOG("About to return topology {}" , _rep_(ti->second));
      return ti->second;
    }
  }
  if (this->_Topologies.size() == 0 ) {
    SIMPLE_ERROR("There are no topology(s) in the constiutition {} - so we won't find any topology(s) that match the monomer {} environment - did you fail to add topology(s) to the constitution?" , _rep_(this->asSmartPtr()) , _rep_(mon));
  }
  SIMPLE_ERROR("There is no topology in constitution {} that matches the monomer {} in its environment - there are {} topology(s) in the constitution" , _rep_(this->asSmartPtr()) , _rep_(mon) , this->_Topologies.size() );
}


CL_LISPIFY_NAME("getTopologyForContext");
CL_DEFMETHOD Topology_sp	Constitution_O::getTopologyForContext(MonomerContext_sp cont )
{
    TopologyMap::iterator	ti;
Topology_sp			tres;
tres = nil<Topology_O>();
    for ( ti=this->_Topologies.begin(); ti!=this->_Topologies.end(); ti++ ) {
	if ( ti->second->matchesContext(cont) ) {
	    return ti->second;
	}
    }
    return tres;
}




#if 0
CL_LISPIFY_NAME("asGroup");
CL_DEFMETHOD RepresentedEntityNameSet_sp	Constitution_O::asGroup()
{
RepresentedEntityNameSet_sp				group;
gctools::Vec0<Stereoisomer_sp>::iterator	si;
    group = RepresentedEntityNameSet_O::create(); //,getCandoDatabase());

//    group->setName(this->getName());
    for (si=this->begin_Stereoisomers(); si!=this->end_Stereoisomers(); si++){
	group->addMonomerName((*si)->getName());
    }
    return group;
}
#endif





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
{
    gctools::Vec0<Topology_sp>::iterator	i;

#ifdef	DEBUG_ON
    LOG("==========  Comparing coupling map VVVVVVV" ); //
    cm->dump();
#endif

    for ( i=this->links.begin(); i!= this->links.end(); i++ ) {
#ifdef	DEBUG_ON
    LOG("----------  To constitution coupling map" ); //
    (*i)->dump();
#endif

	if ( (*i)->matchesTopology(cm) ) {
	    return true;
            LOG("      YES!!!" ); //
        }
    }
    LOG("      Nope!!!" ); //
    return false;
}

#else //] [


#endif //]



void	Constitution_O::initialize()
{
    this->Base::initialize();
    this->_Comment = core::Str_O::create("");
    this->_Name = nil<core::Symbol_O>();
    this->_ConstitutionAtoms = ConstitutionAtoms_O::create();
}

//
// Destructor
//


void Constitution_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,name),this->_Name);
  node->field(INTERN_(kw,comment),this->_Comment);
  node->field(INTERN_(kw,constitution_atoms),this->_ConstitutionAtoms);
  node->field(INTERN_(kw,plugs),this->_PlugsByName);
  node->field(INTERN_(kw,topologies),this->_Topologies);
  node->field(INTERN_(kw,stereo_information),this->_StereoInformation);
  this->Base::fields(node);
}

#ifdef XML_ARCHIVE
void	Constitution_O::archiveBase(core::ArchiveP node)
{
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















}; // namespace chem
#endif
