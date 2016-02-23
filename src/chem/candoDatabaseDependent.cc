/*
    File: candoDatabaseDependent.cc
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
#define	DEBUG_LEVEL_FULL

#include <cando/chem/candoDatabaseDependent.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/candoDatabase.h>
#include <cando/chem/constitution.h>
#include <clasp/core/wrappers.h>


namespace chem {


//
// Constructor
//

void	CandoDatabase_ODependent::initialize()
{
    this->Base::initialize();
}


//
// Copy Constructor
//
CandoDatabase_ODependent::CandoDatabase_ODependent(const CandoDatabase_ODependent& emr ) : O EnvironmentDependent(emr)
{
}

//
// Destructor
//

    void	CandoDatabase_ODependent::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->needsFinalization();
}


/*! Every dependent of a CandoDatabase must listen for the CandoDatabase_replaceMe signal
 */
#ifdef XML_ARCHIVE
    bool	CandoDatabase_ODependent::loadFinalize(core::ArchiveP node)
{
core::Archive_sp		archive;
    if ( !this->Base::loadFinalize(node) ) return false;
    return true;
}
#endif


Constitution_sp CandoDatabase_ODependent::constitution()
{
    Constitution_sp x = _Nil<Constitution_O>();
    return x;
}


#if 0
void	CandoDatabase_ODependent::setCandoDatabase(CandoDatabase_sp db)
{
CandoDatabase_sp	bdb;
    ASSERTNOTNULL(db);
    core::Lisp_sp env = db->getEnvironment();
    this->setEnvironment(env);
}

bool	CandoDatabase_ODependent::hasCandoDatabase()
{
    core::Lisp_sp env = this->getEnvironment();
    return env->hasCandoDatabase();
}



CandoDatabase_sp	CandoDatabase_ODependent::getCandoDatabase()
{
    return getCandoDatabase();
}

#endif








SEBOOSTPYTHON //[
    PYTHON_CLASS(ChemPkg,CandoDatabase,"","",_lisp)
	.def("setCandoDatabase",&CandoDatabase_ODependent::setCandoDatabase)
	.def("getCandoDatabase",&CandoDatabase_ODependent::getCandoDatabase)
    ;
#endif //]
}

};
