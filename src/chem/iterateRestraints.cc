#define	DEBUG_LEVEL_FULL

#include <cando/chem/iterateRestraints.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/restraint.h>
#include <clasp/core/wrappers.h>


namespace chem {


IterateRestraints_sp IterateRestraints_O::create(core::Lisp_sp e,Matter_sp matter)
{_G();
    GC_ALLOCATE(IterateRestraints_O, m );
    m->initTop(matter);
    return m;
};


void	IterateRestraints_O::initialize()
{
    this->Base::initialize();
}


void	IterateRestraints_O::_accumulateRestraints(Matter_sp m)
{_G();
    RestraintList_sp mine = m->allRestraints();
    RestraintList_O::iterator it;
    for ( it=mine->begin(); it!=mine->end(); it++ )
    {
	core::Cons_sp one = core::Cons_O::create((*it),this->_Restraints);
	this->_Restraints = one;
    }
}






void	IterateRestraints_O::initTop(Matter_sp m)
{_G();
    this->_Restraints = _Nil<core::Cons_O>();
    this->_accumulateRestraints(m);
}

void	IterateRestraints_O::first()
{_G();
    this->_CurRestraint = this->_Restraints;
}

void	IterateRestraints_O::next()
{_G();
    if ( this->_CurRestraint.notnilp() )
    {
        this->_CurRestraint = this->_CurRestraint->cdr();
    }
}


    core::T_sp IterateRestraints_O::currentObject()
{_G();
    if ( this->_CurRestraint.nilp() ) return _Nil<core::T_O>();
    return core::oCar(this->_CurRestraint);
}

bool	IterateRestraints_O::isDone()
{_G();
    ASSERTNOTNULL(this->_CurRestraint);
    return this->_CurRestraint.nilp();
}



    void IterateRestraints_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<IterateRestraints_O>()
	;
//	def("create_IterateRestraints",&IterateRestraints_O::create);
}

    void IterateRestraints_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,IterateRestraints,"","",_lisp)
    ;
//    boost::python::def("create_IterateRestraints",&IterateRestraints_O::create);

#endif //]
}



    EXPOSE_CLASS(chem,IterateRestraints_O);
};




