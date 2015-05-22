
#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "container.h"
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include "core/lisp.h"
#include "core/wrappers.h"



namespace candoBase {





    void Container_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<Container_O>()
	.def("size",&Container_O::size)
	.def("isEmpty",&Container_O::isEmpty)
	.def("notEmpty",&Container_O::notEmpty)
	;
}

    void Container_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(CorePkg,Container,"","",_lisp)
	.def("size",&Container_O::size)
	.def("isEmpty",&Container_O::isEmpty)
	.def("notEmpty",&Container_O::notEmpty)
    ;
#endif //]
}



    EXPOSE_CLASS(candoBase,Container_O);
};




