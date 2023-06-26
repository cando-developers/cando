/*
    File: serializer.cc
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
    {
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
    {
	OpenMM::System* rawSystem = NULL;
	boost::filesystem::ifstream fin;
	fin.open(path->getPath(),ios::in);
	if (fin.good() )
	{
	    rawSystem = OpenMM::XmlSerializer::deserialize<OpenMM::System>(fin);
	} else
	{
	    SIMPLE_ERROR("Could not open file[{}] for reading" , path->asString() );
	}
	fin.close();
	ASSERTF(rawSystem!=NULL,("The system was NULL"));
	System_sp system = System_O::create(_lisp);
	system->takeOwnershipOfSystem(rawSystem);
	return system;
    }


    void initializeSerializers(core::LispPtr lisp)
    {
	Defun(systemSerialize);
	Defun(systemDeserialize);
    }

};
