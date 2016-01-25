#define	DEBUG_LEVEL_FULL

#include <cando/chem/iterateRestraints.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/restraint.h>
#include <clasp/core/wrappers.h>


namespace chem {


IterateRestraints_sp IterateRestraints_O::create(Matter_sp matter)
{
    GC_ALLOCATE(IterateRestraints_O, m );
    m->initTop(matter);
    return m;
};


void	IterateRestraints_O::initialize()
{
    this->Base::initialize();
}


void	IterateRestraints_O::_accumulateRestraints(Matter_sp m)
{
    RestraintList_sp mine = m->allRestraints();
    RestraintList_O::iterator it;
    for ( it=mine->begin(); it!=mine->end(); it++ )
    {
	core::Cons_sp one = core::Cons_O::create((*it),this->_Restraints);
	this->_Restraints = one;
    }
}






void	IterateRestraints_O::initTop(Matter_sp m)
{
    this->_Restraints = _Nil<core::T_O>();
    this->_accumulateRestraints(m);
}

void	IterateRestraints_O::first()
{
    this->_CurRestraint = this->_Restraints;
}

void	IterateRestraints_O::next()
{
    if ( this->_CurRestraint.notnilp() )
    {
      this->_CurRestraint = oCdr(this->_CurRestraint);
    }
}


    core::T_sp IterateRestraints_O::currentObject()
{
    if ( this->_CurRestraint.nilp() ) return _Nil<core::T_O>();
    return core::oCar(this->_CurRestraint);
}

bool	IterateRestraints_O::isDone()
{
    ASSERTNOTNULL(this->_CurRestraint);
    return this->_CurRestraint.nilp();
}







};




