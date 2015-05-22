

#include <openmm/include/openmm/serialization/XmlSerializer.h>
#include <boost/filesystem/fstream.hpp>
#include <foundation.h>
#include <object.h>
#include <lisp.h>
#include <core/vector3.h>
#include <cando/openmmPackage/openmmPackage.h>
#include <clasp/core/fileSystem.h>
#include <cando/openmmPackage/openmmExpose.h>
#include <cando/openmmPackage/serializer.h>
#include <clasp/core/wrappers.h>




namespace omm
{



    
    
#define ARGS_af_systemSerialize "(system path)"
#define DECL_af_systemSerialize ""
#define DOCS_af_systemSerialize "systemSerialize"
    core::T_sp af_systemSerialize(System_sp system, core::Path_sp path)
    {_G();
	boost::filesystem::ofstream fout;
	fout.open(path->getPath(),ios::out);
	if (fout.good() )
	    OpenMM::XmlSerializer::serialize(system->wrappedPtr(),"OpenMMSystem",fout);
	fout.close();
	return _lisp->onil();
    }




    
    
#define ARGS_af_systemDeserialize "(path)"
#define DECL_af_systemDeserialize ""
#define DOCS_af_systemDeserialize "systemDeserialize"
    System_sp af_systemDeserialize(core::Path_sp path)
    {_G();
	OpenMM::System* rawSystem = NULL;
	boost::filesystem::ifstream fin;
	fin.open(path->getPath(),ios::in);
	if (fin.good() )
	{
	    rawSystem = OpenMM::XmlSerializer::deserialize<OpenMM::System>(fin);
	} else
	{
	    SIMPLE_ERROR(BF("Could not open file[%s] for reading") % path->asString() );
	}
	fin.close();
	ASSERTF(rawSystem!=NULL,BF("The system was NULL"));
	System_sp system = System_O::create(_lisp);
	system->takeOwnershipOfSystem(rawSystem);
	return system;
    }


    void initializeSerializers(core::Lisp_sp lisp)
    {_G();
	Defun(systemSerialize);
	Defun(systemDeserialize);
    }

};
