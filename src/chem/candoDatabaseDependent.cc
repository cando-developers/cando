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


