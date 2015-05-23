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



#ifdef CONSPACK
    void	InfoDb_O::archiveBase(core::ArchiveP node)
{
    node->attributeStringMapIfDefined(KW("k"),this->_database);
}
#endif

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


#ifdef CONSPACK
    void	ForceField_O::archiveBase(core::ArchiveP node)
{
    node->attribute("title",this->_Title);
    node->attribute("ref",this->_Ref);
    node->attributeVector0( KW("bondDistinctions"), this->_SingleBondMultiBondDistinctions );
    node->attribute( "_Info", this->_Info );
    node->attribute( "_Types", this->_Types );
    node->attribute( "_Stretches", this->_Stretches );
    node->attribute( "_Angles", this->_Angles );
    node->attribute( "_Itors", this->_Itors );
    node->attribute( "_Ptors", this->_Ptors );
    node->attribute( "_Nonbonds", this->_Nonbonds );
    node->attribute( "_Vdws", this->_Vdws );
}
#endif


#ifdef XML_ARCHIVE
void	ForceField_O::saveAs(const string& fileName)
{_G();
    core::XmlSaveArchive_sp	xml;
    xml = core::XmlSaveArchive_O::create();
    xml->put("forceField",this->sharedThis<ForceField_O>());
    xml->saveAs(fileName);
}
#endif

void	ForceField_O::assignTypes(Matter_sp matter)
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
    if ( this->_Types.pointerp() ) {
	this->_Types->setForceField(this->sharedThis<ForceField_O>());
    }
}


void	ForceField_O::setFFStretchDb( FFStretchDb_sp Stretches)
{_G();
    this->_Stretches = Stretches;
    if (this->_Stretches.pointerp()) this->_Stretches->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFAngleDb( FFAngleDb_sp Angles)
{_G();
    this->_Angles = Angles;
    if (this->_Angles.pointerp()) this->_Angles->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFItorDb( FFItorDb_sp Itors)
{_G();
    this->_Itors = Itors;
    if (this->_Itors.pointerp() ) this->_Itors->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFPtorDb( FFPtorDb_sp Ptors)
{_G();
    this->_Ptors = Ptors;
    if (this->_Ptors.pointerp()) this->_Ptors->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFNonbondDb(FFNonbondDb_sp Nonbonds )
{_G();
    this->_Nonbonds = Nonbonds;
    if (this->_Nonbonds.pointerp())this->_Nonbonds->setForceField(this->sharedThis<ForceField_O>());
}
void	ForceField_O::setFFVdwDb(FFVdwDb_sp Vdws )
{_G();
    this->_Vdws = Vdws;
    if (this->_Vdws.pointerp())this->_Vdws->setForceField(this->sharedThis<ForceField_O>());
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


