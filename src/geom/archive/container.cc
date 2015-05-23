
#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/geom/container.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/lisp.h>
#include <clasp/core/wrappers.h>



namespace geom {





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



    EXPOSE_CLASS(geom,Container_O);
};




