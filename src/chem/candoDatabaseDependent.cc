#define	DEBUG_LEVEL_FULL

#include "candoDatabaseDependent.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "candoDatabase.h"
#include "constitution.h"
#include "core/wrappers.h"


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
{_G();
    this->Base::archiveBase(node);
    node->needsFinalization();
}


/*! Every dependent of a CandoDatabase must listen for the CandoDatabase_replaceMe signal
 */
#ifdef XML_ARCHIVE
    bool	CandoDatabase_ODependent::loadFinalize(core::ArchiveP node)
{_G();
core::Archive_sp		archive;
    if ( !this->Base::loadFinalize(node) ) return false;
    return true;
}
#endif


Constitution_sp CandoDatabase_ODependent::constitution()
{_G();
    Constitution_sp x = _Nil<Constitution_O>();
    return x;
}


#if 0
void	CandoDatabase_ODependent::setCandoDatabase(CandoDatabase_sp db)
{_G();
CandoDatabase_sp	bdb;
    ASSERTNOTNULL(db);
    core::Lisp_sp env = db->getEnvironment();
    this->setEnvironment(env);
}

bool	CandoDatabase_ODependent::hasCandoDatabase()
{_G();
    core::Lisp_sp env = this->getEnvironment();
    return env->hasCandoDatabase();
}



CandoDatabase_sp	CandoDatabase_ODependent::getCandoDatabase()
{_G();
    core::Lisp_sp env = this->getEnvironment();
    ASSERTP(env.notnilp(),"My environment is nil");
    return getCandoDatabase(env);
}

#endif






    void CandoDatabase_ODependent::exposeCando(core::Lisp_sp lisp)
{
    core::class_<CandoDatabase_ODependent>()
	.def("setCandoDatabase",&CandoDatabase_ODependent::setCandoDatabase)
	.def("getCandoDatabase",&CandoDatabase_ODependent::getCandoDatabase)
    ;
}

    void CandoDatabase_ODependent::exposePython()
    {_FG);
#ifdef	USEBOOSTPYTHON //[
    PYTHON_CLASS(ChemPkg,CandoDatabase,"","",_lisp)
	.def("setCandoDatabase",&CandoDatabase_ODependent::setCandoDatabase)
	.def("getCandoDatabase",&CandoDatabase_ODependent::getCandoDatabase)
    ;
#endif //]
}

    EXPOSE_CLASS(chem,CandoDatabase_ODependent);
};


