#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <cando/candoBase/range.h>
//#include "core/archiveNode.h"
#include <clasp/core/symbolTable.h>
#include <clasp/core/multipleValues.h>
//#include "core/archive.h"
#include <clasp/core/lisp.h>
#include <clasp/core/environment.h>
#include <clasp/core/wrappers.h>


namespace candoBase {

#if 0
    core::T_sp Range_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp)
    {_G();
	LOG(BF("Argument types %s %s") % core::oCar(args)->__class()->classNameAsString() % args->ocadr()->className() );
	SYMBOL_SC_(CandoBasePkg,begin);
	SYMBOL_SC_(CandoBasePkg,end);
	SYMBOL_SC_(CandoBasePkg,step);
	this->_Begin = environ->lookup(_sym_begin)->object().as<Integer_O>()->as_core::LongLongInt();
	this->_End = environ->lookup(_sym_end)->object().as<Integer_O>()->as_core::LongLongInt();
	this->_Step = environ->lookup(_sym_step)->object().as<Integer_O>()->as_core::LongLongInt();
	return _Nil<core::T_O>();
    }
#endif



    Range_sp Range_O::create(core::Lisp_sp lisp,core::LongLongInt b, core::LongLongInt e)
{_G();
    LOG(BF("Range_O::create with begin=%d, end=%d") % b % e  );
    GC_ALLOCATE(Range_O,r);
    r->_Begin = b;
    r->_End = e;
    r->_Step = 1;
    return r;
};
    Range_sp Range_O::create3(core::Lisp_sp lisp, core::LongLongInt b, core::LongLongInt e, core::LongLongInt step)
{_G();
    LOG(BF("Range_O::create with begin=%d, end=%d, step=%d") % b % e % step  );
    GC_ALLOCATE(Range_O,r);
    r->_Begin = b;
    r->_End = e;
    r->_Step = step;
    return r;
};



    core::Cons_sp Range_O::asCons()
{_G();
    core::Cons_sp dummy = core::Cons_O::create(_Nil<core::T_O>(),_Nil<core::Cons_O>());
    core::Cons_sp cur = dummy;
    for (int i=this->_Begin; i<this->_End; i+=this->_Step )
    {
	core::Cons_sp one = core::Cons_O::create(core::Integer_O::create(i),_Nil<core::Cons_O>());
	cur->setCdr(one);
	cur = one;
    }
    return dummy->cdr();
}



string	Range_O::__repr__() const
{
    stringstream ss;
    ss << "Range(";
    ss << this->_Begin;
    ss << ",";
    ss << this->_End;
    ss << ",";
    ss << this->_Step;
    ss << ")";
    return ss.str();
}

void	Range_O::initialize()
{
    this->Base::initialize();
    this->_Begin = 0;
    this->_End = 10;
    this->_Step = 1;
}

#ifdef XML_ARCHIVE
    void	Range_O::archive(core::ArchiveP node)
{
    IMPLEMENT_ME();
}
#endif


    void Range_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<Range_O>()
	.def("asCons",&Range_O::asCons)
	;
//	def("oldRange",&Range_O::create,_lisp);
}

    void Range_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(CorePkg,Range,"","",_lisp)
	.def("asCons",&Range_O::asCons)
	;
#endif
//	def("oldRange",&Range_O::create,_lisp);
}

    EXPOSE_CLASS(candoBase,Range_O);

};




