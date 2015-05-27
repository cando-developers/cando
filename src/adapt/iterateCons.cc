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


    void IterateCons_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<IterateCons_O>()
		.def("init",&IterateCons_O::init)
	;
//	def("create_IterateCons",&IterateCons_O::create);
}

    void IterateCons_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(CorePkg,IterateCons,"","",_lisp)
    ;
#endif //]
}

    EXPOSE_CLASS(adapt,IterateCons_O);
};




