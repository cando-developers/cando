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
{_G();
    core::XmlSaveArchive_sp	xml;
    xml = core::XmlSaveArchive_O::create();
    xml->put("forceField",this->sharedThis<ForceField_O>());
    xml->saveAs(fileName);
}
#endif

CL_NAME("assignTypes");
CL_DEFMETHOD void	ForceField_O::assignTypes(Matter_sp matter)
{_G();
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
{_G();
    this->_Types = Types;
    if ( this->_Types.notnilp() ) {
	this->_Types->setForceField(this->sharedThis<ForceField_O>());
    }
}


void	ForceField_O::setFFStretchDb( FFStretchDb_sp Stretches)
{_G();
    this->_Stretches = Stretches;
    if (this->_Stretches.generalp()) this->_Stretches->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFAngleDb( FFAngleDb_sp Angles)
{_G();
    this->_Angles = Angles;
    if (this->_Angles.generalp()) this->_Angles->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFItorDb( FFItorDb_sp Itors)
{_G();
    this->_Itors = Itors;
    if (this->_Itors.generalp() ) this->_Itors->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFPtorDb( FFPtorDb_sp Ptors)
{_G();
    this->_Ptors = Ptors;
    if (this->_Ptors.generalp()) this->_Ptors->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFNonbondDb(FFNonbondDb_sp Nonbonds )
{_G();
    this->_Nonbonds = Nonbonds;
    if (this->_Nonbonds.generalp())this->_Nonbonds->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFVdwDb(FFVdwDb_sp Vdws )
{_G();
    this->_Vdws = Vdws;
    if (this->_Vdws.generalp())this->_Vdws->setForceField(this->sharedThis<ForceField_O>());
}










    void ForceField_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<ForceField_O>()
    	.def("getTypes",&ForceField_O::getTypes)
	.def("getStretchDb",&ForceField_O::getStretchDb)
	.def("getAngleDb",&ForceField_O::getAngleDb)
	.def("getItorDb",&ForceField_O::getItorDb)
	.def("getPtorDb",&ForceField_O::getPtorDb)
    	.def("getNonbondDb",&ForceField_O::getNonbondDb)
	.def("getVdwDb",&ForceField_O::getVdwDb)
	.def("assignTypes", &ForceField_O::assignTypes)
    ;
//    defInPackage(ChemPkg,"loadForceField",&ForceField_O::open_ForceField);
}

void ForceField_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef	USEBOOSTPYTHON //[
#ifdef USEBOOSTPYTHON
    PYTHON_CLASS(ChemPkg,ForceField,"","",_lisp)
	.def("assignTypes", &ForceField_O::assignTypes)
	.def("setInfoDb",&ForceField_O::setInfoDb)
	.def("setTitle",&ForceField_O::setTitle)
	.def("setFFTypeDb",&ForceField_O::setFFTypeDb)
	.def("setFFStretchDb",&ForceField_O::setFFStretchDb)
	.def("setFFAngleDb",&ForceField_O::setFFAngleDb)
	.def("setFFPtorDb",&ForceField_O::setFFPtorDb)
	.def("setFFItorDb",&ForceField_O::setFFItorDb)
	.def("setFFNonbondDb",&ForceField_O::setFFNonbondDb)
	.def("setFFVdwDb",&ForceField_O::setFFVdwDb)

	.def_readonly("_Info",&ForceField_O::_Info)
	.def_readonly("_Title",&ForceField_O::_Title)
	.def_readonly("_Types",&ForceField_O::_Types)
	.def_readonly("_Stretches",&ForceField_O::_Stretches)
	.def_readonly("_Angles",&ForceField_O::_Angles)
	.def_readonly("_Ptors",&ForceField_O::_Ptors)
	.def_readonly("_Itors",&ForceField_O::_Itors)
	.def_readonly("_Nonbonds",&ForceField_O::_Nonbonds)
	.def_readonly("_Vdws",&ForceField_O::_Vdws)
//	.def("saveAs", &ForceField_O::saveAs,"","","",false)
	;
#endif
//    boost::python::def("open_ForceField",&ForceField_O::open_ForceField);
#endif //]
}

REGISTER_CLASS(chem, InfoDb_O );
EXPOSE_CLASS(chem, ForceField_O );
};


