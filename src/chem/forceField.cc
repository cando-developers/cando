/*
    File: forceField.cc
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


/*
 *	forceField.cc
 *
 *	Maintains a parameter set for the Moe general force field
 *	from these parameters we can derive any force field parameters
 *
 */

#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <clasp/core/ql.h>
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

void InfoDb_O::addInfo( core::Symbol_sp key, core::String_sp data )
{
  if ( this->_database.count(key)==0 ) {
    this->_database.set(key,data);
  } else {
    this->_database.set(key,core::SimpleBaseString_O::make(this->_database.get(key)->get_std_string()+"\n"+data->get_std_string()));
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


CL_DEFMETHOD void ForceField_O::forceFieldMerge(ForceField_sp other)
{
  this->_Info->forceFieldMerge(other->_Info);
  this->_Types->forceFieldMerge(other->_Types);
  this->_Stretches->forceFieldMerge(other->_Stretches);
  this->_Angles->forceFieldMerge(other->_Angles);
  this->_Itors->forceFieldMerge(other->_Itors);
  this->_Ptors->forceFieldMerge(other->_Ptors);
  this->_Nonbonds->forceFieldMerge(other->_Nonbonds);
  this->_Vdws->forceFieldMerge(other->_Vdws);
}

CL_LISPIFY_NAME("assignTypes");
CL_DEFMETHOD void	ForceField_O::assignTypes(Matter_sp matter)
{
  FFTypesDb_sp types = this->getTypes();
  types->assignTypes(matter);
}


string ForceField_O::__repr__() const {
  stringstream ss;
  ss << "#<FORCE-FIELD ";
  ss << this->_Title;
  ss << ">";
  return ss.str();
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
}


void	ForceField_O::setFFStretchDb( FFStretchDb_sp Stretches)
{
  this->_Stretches = Stretches;
}

void	ForceField_O::setFFAngleDb( FFAngleDb_sp Angles)
{
  this->_Angles = Angles;
}
void	ForceField_O::setFFItorDb( FFItorDb_sp Itors)
{
  this->_Itors = Itors;
}
void	ForceField_O::setFFPtorDb( FFPtorDb_sp Ptors)
{
  this->_Ptors = Ptors;
}
void	ForceField_O::setFFNonbondDb(FFNonbondDb_sp Nonbonds )
{
  this->_Nonbonds = Nonbonds;
}
void	ForceField_O::setFFVdwDb(FFVdwDb_sp Vdws )
{
  this->_Vdws = Vdws;
}


CL_LISPIFY_NAME("make-CombinedForceField");
CL_DEF_CLASS_METHOD
CombinedForceField_sp CombinedForceField_O::make() {
  GC_ALLOCATE(CombinedForceField_O,ff);
  return ff;
};

void	CombinedForceField_O::fields(core::Record_sp node)
{
  node->field( INTERN_(kw,force_fields),this->_ForceFields);
}



CL_LISPIFY_NAME(CombinedForceField_addShadowingForceField);
CL_DEFMETHOD
void CombinedForceField_O::addShadowingForceField(ForceField_sp forceField, core::T_sp info)
{
  core::Cons_sp pair = core::Cons_O::create(forceField,info);
  this->_ForceFields = core::Cons_O::create(pair,this->_ForceFields);
}


CL_LISPIFY_NAME(CombinedForceField_forceFieldsAsList);
CL_DEFMETHOD
core::List_sp CombinedForceField_O::forceFieldsAsList() const {
  ql::list result;
  for ( auto cur : this->_ForceFields) {
    core::Cons_sp pair = gc::As_unsafe<core::Cons_sp>(CONS_CAR(cur));
    result << CONS_CAR(pair);
  }
  return result.cons();
}


CL_LISPIFY_NAME(CombinedForceField_assignForceFieldTypes);
CL_DEFMETHOD
void CombinedForceField_O::assignForceFieldTypes(Matter_sp molecule) {
  FFTypesDb_sp fftypes = FFTypesDb_O::create();
  core::List_sp parts = this->forceFieldsAsList();
  for ( auto cur : parts ) {
    ForceField_sp other = gc::As<ForceField_sp>(CONS_CAR(cur));
    FFTypesDb_sp other_fftypes = other->getTypes();
    fftypes->forceFieldMerge(other_fftypes);
  }
  fftypes->assignTypes(molecule);
}

};
