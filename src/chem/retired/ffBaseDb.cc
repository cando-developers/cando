#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


/*!
 *	ffBase.cc
 *
 *	Maintain the base class for parameters
 */
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/forceField.h>
//#include "core/archiveNode.h"
#include <clasp/core/wrappers.h>

namespace chem {






string	FFParameter_O::levelDescription()
{
    switch ( this->_Level ) {
	case parameterized:
	    return "parameterized ";
	case estimated:
	    return "estimated ";
	case rough:
	    return "guessed ";
    }
    return "unknownParameterLevel ";
}


    core::NullTerminatedEnumAssociation paramLevelEnum[] = {
	{ "parameterized", parameterized },
	{ "estimated", estimated },
	{ "rough", rough },
	{ "", -1 }
};

    void	FFParameter_O::archiveBase(core::ArchiveP node)
    {
	node->attributeEnum( KW("level"), this->_Level);
    }

#ifdef CONSPACK
    void	FFBaseDb_O::archiveBase(core::ArchiveP node)
    {
	node->attributeWeakPointer( KW("forceField"), this->_ForceField );
    }
#endif

void	FFBaseDb_O::setForceField(ForceField_sp ff)
{_G();
    ASSERTNOTNULL(ff);
    this->_ForceField = ff;
};



    void FFParameter_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<FFParameter_O>()
    ;
}
    void FFParameter_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFParameter,"","",_lisp)
//        .def_readwrite("_Level",&FFParameter_O::_Level)
	.def("levelDescription",&FFParameter_O::levelDescription)
	;
#endif
#endif
}





    void FFBaseDb_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<FFBaseDb_O>()
    ;
}
    void FFBaseDb_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFBaseDb,"","",_lisp)
	.def("setForceField",&FFBaseDb_O::setForceField)
	;
#endif
}



    EXPOSE_CLASS(chem, FFParameter_O);
    EXPOSE_CLASS(chem, FFBaseDb_O);


};

