#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	forceField.cc
 *
 *	Maintains a parameter set for the Moe general force field
 *	from these parameters we can derive any force field parameters
 *
 */

#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/chem/forceField.h>
//#include "core/archive.h"
//#include "core/xmlLoadArchive.h"
//#include "core/xmlSaveArchive.h"
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>



namespace chem {

//extern  bool    moeForceFieldCompile(ForceField_O* ff, string fileName );



void	InfoDb_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,k),this->_database);
}

    void InfoDb_O::addInfo( core::Symbol_sp key, const string& data )
    {
	if ( this->_database.count(key)==0 ) {
	    this->_database.set(key,core::Str_O::create(data));
	} else {
	    this->_database.set(key,core::Str_O::create(this->_database.get(key)->get()+"\n"+data));
	}
    };




void	ForceField_O::initialize()
{
    this->Base::initialize();
    LOG(BF("Initializing ForceField") );
    this->_Info = InfoDb_O::create();
    LOG(BF("Initializing ForceField") );
    this->_Types = FFTypesDb_O::create();
    LOG(BF("Initializing ForceField") );
    this->_Stretches = FFStretchDb_O::create();
    LOG(BF("Initializing ForceField angle") );
    this->_Angles = FFAngleDb_O::create();
    LOG(BF("Initializing ForceField itor") );
    this->_Itors = FFItorDb_O::create();
    LOG(BF("Initializing ForceField ptor") );
    this->_Ptors = FFPtorDb_O::create();
    LOG(BF("Initializing ForceField nonbonded") );
    this->_Nonbonds = FFNonbondDb_O::create();
    LOG(BF("Initializing ForceField vdw") );
    this->_Vdws = FFVdwDb_O::create();
};


void	ForceField_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,title),this->_Title);
  node->field( INTERN_(kw,ref),this->_Ref);
  node->field( INTERN_(kw,bondDistinctions), this->_SingleBondMultiBondDistinctions );
  node->field( INTERN_(kw,Info), this->_Info );
  node->field( INTERN_(kw,Types), this->_Types );
  node->field( INTERN_(kw,Stretches), this->_Stretches );
  node->field( INTERN_(kw,Angles), this->_Angles );
  node->field( INTERN_(kw,Itors), this->_Itors );
  node->field( INTERN_(kw,Ptors), this->_Ptors );
  node->field( INTERN_(kw,Nonbonds), this->_Nonbonds );
  node->field( INTERN_(kw,Vdws), this->_Vdws );
}


#ifdef XML_ARCHIVE
void	ForceField_O::saveAs(const string& fileName)
{
    core::XmlSaveArchive_sp	xml;
    xml = core::XmlSaveArchive_O::create();
    xml->put("forceField",this->sharedThis<ForceField_O>());
    xml->saveAs(fileName);
}
#endif

CL_LISPIFY_NAME("assignTypes");
CL_DEFMETHOD void	ForceField_O::assignTypes(Matter_sp matter)
{
FFTypesDb_sp	types;
    types = this->getTypes();
    types->assignTypes(matter);
}



void	ForceField_O::setTitle(const string& title)
{
    this->_Title = title;
}

void	ForceField_O::setInfoDb( InfoDb_sp Info )
{
    this->_Info = Info;
}

void	ForceField_O::setFFTypeDb( FFTypesDb_sp Types)
{
  this->_Types = Types;
  if (this->_Types.notnilp()) this->_Types->setForceField(this->sharedThis<ForceField_O>());
}


void	ForceField_O::setFFStretchDb( FFStretchDb_sp Stretches)
{
    this->_Stretches = Stretches;
    if (this->_Stretches.notnilp()) this->_Stretches->setForceField(this->asSmartPtr());
}

void	ForceField_O::setFFAngleDb( FFAngleDb_sp Angles)
{
    this->_Angles = Angles;
    if (this->_Angles.notnilp()) this->_Angles->setForceField(this->asSmartPtr());
}
void	ForceField_O::setFFItorDb( FFItorDb_sp Itors)
{
    this->_Itors = Itors;
    if (this->_Itors.notnilp()) this->_Itors->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFPtorDb( FFPtorDb_sp Ptors)
{
    this->_Ptors = Ptors;
    if (this->_Ptors.notnilp()) this->_Ptors->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFNonbondDb(FFNonbondDb_sp Nonbonds )
{
    this->_Nonbonds = Nonbonds;
    if (this->_Nonbonds.notnilp()) this->_Nonbonds->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFVdwDb(FFVdwDb_sp Vdws )
{
    this->_Vdws = Vdws;
    if (this->_Vdws.notnilp()) this->_Vdws->setForceField(this->sharedThis<ForceField_O>());
}



};


