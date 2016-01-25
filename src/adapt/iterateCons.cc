#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/adapt/iterateCons.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/lisp.h>
#include <clasp/core/cons.h>
#include <clasp/core/wrappers.h>


namespace adapt {





void	IterateCons_O::initialize()
{
    this->Base::initialize();
}


void IterateCons_O::first() 
{
    this->_Cur = this->_Begin;
};
void IterateCons_O::next() 
{
  this->_Cur = oCdr(this->_Cur);
};
bool IterateCons_O::isDone() 
{_OF();
    ASSERTNOTNULL(this->_Cur);
    return (this->_Cur.nilp());
};

    core::T_sp IterateCons_O::currentObject() 
{ 
    return core::oCar(this->_Cur);
};




};




