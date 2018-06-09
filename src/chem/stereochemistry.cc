/*
    File: stereochemistry.cc
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
       
//
// (C) 2004 Christian E. Schafmeister
//

#define	TURN_DEBUG_OFF

#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/environment.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/stringList.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/constitution.h>
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/restraint.h>
#include <cando/chem/complexRestraints.h>
#include <cando/chem/entity.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/wrappers.h>




namespace chem {

/*! Return +1 for "R" and -1 for "S"
*/
CL_LISPIFY_NAME("getMoeConfiguration");
CL_DEFMETHOD int	StereoConfiguration_O::getMoeConfiguration()
{
  if ( this->_Configuration == chem::_sym_R ) return 1;
  else if (this->_Configuration == chem::_sym_S ) return -1;
  SIMPLE_ERROR(BF("Configuration can only be chem:S or chem:R - it was %s") % _rep_(this->_Configuration));
}


SYMBOL_EXPORT_SC_(ChemKwPkg,R);
SYMBOL_EXPORT_SC_(ChemKwPkg,S);

#ifdef XML_ARCHIVE
void	StereoConfiguration_O::archive(core::ArchiveP node)
{
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

StereoConfiguration_sp StereoConfiguration_O::make(core::Symbol_sp atomName, core::Symbol_sp config)
{
  GC_ALLOCATE(StereoConfiguration_O, me );
  me->_AtomName = atomName;
  me->_Configuration = config;
  return me;
};

    
    
core::List_sp StereoConfiguration_O::stereochemical_permutations(uint numberOfCenters)
{
  core::List_sp list = _Nil<core::T_O>();
  uint pow = 1 << numberOfCenters;
  for ( uint i=0; i<pow; i++ )
  {
    uint map = i;
    core::List_sp one = _Nil<core::T_O>();
    for ( uint ix=0; ix<numberOfCenters; ix++ )
    {
      if ( map & 1 )
      {
        one = core::Cons_O::create(chem::_sym_R,one);
      } else
      {
        one = core::Cons_O::create(chem::_sym_S,one);
      }
      map = map >> 1;
    }
    list = core::Cons_O::create(one,list);
  }
  return list;
}

core::List_sp StereoConfiguration_O::create_multiple(core::List_sp atomNames, core::List_sp configurations)
{
  ASSERT_eq(core::cl__length(atomNames), core::cl__length(configurations) );
    core::List_sp list = _Nil<core::T_O>();
    core::List_sp curName = atomNames;
    core::List_sp curConfig = configurations;
    while ( curName.notnilp() )
    {
	GC_ALLOCATE(StereoConfiguration_O, one );
	one->setAtomName(oCar(curName).as<core::Symbol_O>());
	one->setConfiguration(oCar(curConfig).as<core::Symbol_O>());
	list = core::Cons_O::create(one,list);
	curName = oCdr(curName);
	curConfig = oCdr(curConfig);
    }
    return list;
}

void StereoConfiguration_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,name),this->_AtomName);
  node->field(INTERN_(kw,config),this->_Configuration);
  this->Base::fields(node);
}


void Stereoisomer_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,name),this->_Name);
  node->field(INTERN_(kw,pdb),this->_Pdb);
  node->field(INTERN_(kw,ent),this->_Enantiomer);
  node->field(INTERN_(kw,configs),this->_Configurations);
  this->Base::fields(node);
}

adapt::SymbolSet_sp	Stereoisomer_O::expandedNameSet()
{
  adapt::SymbolSet_sp	ss;
  ss = adapt::SymbolSet_O::create();
  ss->insert(this->getName());
  return ss;
}

	/*! Build and return an ObjectList that contains one RepresentedEntityNameSet
	 * that has this Stereoisomers name as its representative and
	 * a ContainedNames list with only this Stereoisomers name
	 */
RepresentativeList_sp Stereoisomer_O::expandedRepresentativeList() const
{
  RepresentativeList_sp ss = RepresentativeList_O::create();
  RepresentedEntityNameSet_sp gr = RepresentedEntityNameSet_O::create();
  gr->setName(this->getName());
  gr->setRepresentative(this->getName());
  gr->addEntityName(this->getName());
  ss = RepresentativeList_O::create();
  ss->vectorPushExtend(gr);
  return ss;
}


string Stereoisomer_O::__repr__() const
{
  stringstream ss;
  ss << "Stereoisomer-name(" << _rep_(this->_Name) << ")-pdbName(" << _rep_(this->_Pdb) << ")";
  return ss.str();
}

CL_LISPIFY_NAME("getConfigurationForCenter");
CL_DEFMETHOD core::Symbol_sp Stereoisomer_O::getConfigurationForCenter( core::Symbol_sp centerName )
{
  for ( gctools::Vec0<StereoConfiguration_sp>::iterator it=this->_Configurations.begin();
        it!=this->_Configurations.end(); it++ )
  {
    if ( (*it)->getAtomName() == centerName )
    {
      return (*it)->getConfiguration();
    }
  }
  return chem::_sym_UnknownConfiguration;
}


/*
__BEGIN_DOC(classes.Stereoisomer.!class.Stereoisomer)
\scriptCmdRet{Stereoisomer}{name:(Text::name) pdb:(Text::pdb) configs:(List::)}{StereoIsomer::}
__END_DOC
*/

Stereoisomer_sp Stereoisomer_O::make(core::Symbol_sp name, core::Symbol_sp pdb, core::List_sp configs)
{
  GC_ALLOCATE_VARIADIC(Stereoisomer_O, me, name, pdb);
  for ( auto cur : configs ) me->_Configurations.push_back(oCar(cur).as<StereoConfiguration_O>() );
  return me;
};


/*
__BEGIN_DOC(classes.Stereoisomer.!class.MultiStereoisomers)
\scriptCmdRet{MultiStereoisomers}{nameTemplate:(Text::template) centers:(List::) configs:(List::)}{Cons::stereoisomers}
__END_DOC
*/




#define ARGS_af_multiStereoisomers "(name-template centers configs)"
#define DECL_af_multiStereoisomers ""
#define DOCS_af_multiStereoisomers "multiStereoisomers"
core::T_sp af_multiStereoisomers(const string& nameTemplate, core::List_sp centers, core::List_sp configs)
{
  IMPLEMENT_MEF("Implement multiStereoisomers");
};



void StereoInformation_O::fields(core::Record_sp node) {
  node->field(INTERN_(kw,stereoisomers),this->_Stereoisomers);
  node->field(INTERN_(kw,namemap),this->_NameOrPdbToStereoisomer);
  node->field(INTERN_(kw,restraints),this->_ComplexRestraints);
  this->Base::fields(node);
}

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
  this->_Stereoisomers.push_back(s);
}


CL_LISPIFY_NAME("validate");
CL_DEFMETHOD void StereoInformation_O::validate()
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


Stereoisomer_sp	StereoInformation_O::getStereoisomer(core::T_sp  nameOrPdb )
{_OF();
  return this->_NameOrPdbToStereoisomer.get(nameOrPdb);
}

#ifdef XML_ARCHIVE
void	StereoInformation_O::archive(core::ArchiveP node)
{
  node->archiveVector0OfObjectsSubClassOf<Stereoisomer_O>(this->_Stereoisomers);
//    node->archiveVector0OfObjectsSubClassOf<O_ProChiralCenter>( this->_ProChiralCenters);
  node->archiveVector0( "complexRestraints", this->_ComplexRestraints);
  node->needsFinalization();
}
#endif


/*
__BEGIN_DOC(classes.StereoInformation.!class.StereoInformation)
\requiredKeyed{stereoisomers:}{List::stereoisomers}
\optionalKeyed{proChiralCenters:}{List::}
\optionalKeyed{constrainedPiBonds:}{List::}
__END_DOC
*/


StereoInformation_sp StereoInformation_O::make(core::List_sp stereoisomers, core::List_sp restraints)
{
  GC_ALLOCATE(StereoInformation_O, me );
  core::fillVec0(stereoisomers,me->_Stereoisomers);
  me->_NameOrPdbToStereoisomer.clear();
  for ( gctools::Vec0<Stereoisomer_sp>::iterator it = me->_Stereoisomers.begin();
        it != me->_Stereoisomers.end(); it++ )
  {
    if ( me->_NameOrPdbToStereoisomer.count((*it)->getName()) >0 )
    {
      SIMPLE_ERROR(BF("Stereoisomer name: "+_rep_((*it)->getName())+" has already been defined"));
    }
    if ( me->_NameOrPdbToStereoisomer.count((*it)->getPdb()) >0 )
    {
      SIMPLE_ERROR(BF("Stereoisomer pdb: "+_rep_((*it)->getPdb())+" has already been defined"));
    }
    me->_NameOrPdbToStereoisomer.set((*it)->getName(), (*it));
    me->_NameOrPdbToStereoisomer.set((*it)->getPdb(), (*it));
  }
  return me;
};


//
//	getMonomerNames
//
//!	Return a adapt::StringList_sp of monomer names for this mold
adapt::StringList_sp StereoInformation_O::getMonomerNamesAsStringList()
{
  adapt::StringList_sp				names;
  gctools::Vec0<Stereoisomer_sp>::iterator	mnpi;
  names = adapt::StringList_O::create();
  for ( mnpi=this->_Stereoisomers.begin();
        mnpi!=this->_Stereoisomers.end(); mnpi++ ) {
    names->append((*mnpi)->getName()->symbolNameAsString());
  }
  return names;
}




//!	Return a adapt::StringSet_sp of monomer names for this Constitution
adapt::SymbolSet_sp StereoInformation_O::getMonomerNamesAsSymbolSet()
{
  adapt::SymbolSet_sp				names;
  gctools::Vec0<Stereoisomer_sp>::iterator	mnpi;
  names = adapt::SymbolSet_O::create();
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
adapt::StringList_sp StereoInformation_O::getPdbNamesAsStringList()
{
  adapt::StringList_sp	names;
  gctools::Vec0<Stereoisomer_sp>::iterator	mnpi;
  names = adapt::StringList_O::create();
  for ( mnpi=this->_Stereoisomers.begin();
        mnpi!=this->_Stereoisomers.end(); mnpi++ ) {
    names->append((*mnpi)->getPdb()->symbolNameAsString());
  }
  return names;
}
SYMBOL_EXPORT_SC_(ChemKwPkg,undefinedConfiguration);
SYMBOL_EXPORT_SC_(ChemKwPkg,S);
SYMBOL_EXPORT_SC_(ChemKwPkg,R);
CL_BEGIN_ENUM(ConfigurationEnum,_sym__PLUS_configurationEnumConverter_PLUS_,"Configuration");
CL_VALUE_ENUM(chemkw::_sym_undefinedConfiguration, undefinedConfiguration );
CL_VALUE_ENUM(chemkw::_sym_S, S_Configuration );
CL_VALUE_ENUM(chemkw::_sym_R, R_Configuration );
CL_END_ENUM(_sym__PLUS_configurationEnumConverter_PLUS_);
;

SYMBOL_EXPORT_SC_(ChemKwPkg,undefinedCenter);
SYMBOL_EXPORT_SC_(ChemKwPkg,chiral);
SYMBOL_EXPORT_SC_(ChemKwPkg,prochiral);
CL_BEGIN_ENUM(StereochemistryType,_sym__PLUS_stereochemistryTypeConverter_PLUS_,"Stereochemistry");
CL_VALUE_ENUM(chemkw::_sym_undefinedCenter,undefinedCenter);
CL_VALUE_ENUM(chemkw::_sym_chiral, chiralCenter);
CL_VALUE_ENUM(chemkw::_sym_prochiral, prochiralCenter);
CL_END_ENUM(_sym__PLUS_stereochemistryTypeConverter_PLUS_);
;

}; //namespace chem
