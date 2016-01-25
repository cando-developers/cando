#define DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/monomerId.h>
#include <clasp/core/wrappers.h>



namespace kinematics
{

    string MonomerId::asString() const
    {
	stringstream ss;
	ss << "Chain["<<this->_Chain;
	ss << "]Res["<<this->_Monomer <<"]";
	return ss.str();
    }







#if INIT_TO_FACTORIES

#define ARGS_MonomerId_O_make "(chain monomer)"
#define DECL_MonomerId_O_make ""
#define DOCS_MonomerId_O_make "make MonomerId"
  MonomerId_sp MonomerId_O::make(const int chain, const int monomer)
    {
	GC_ALLOCATE(MonomerId_O, me );
      me->_MonomerId._Chain = chain;
      me->_MonomerId._Monomer = monomer;
      return me;
    };

#else
    core::T_sp MonomerId_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {
	this->_MonomerId._Chain = translate::from_object<int>::convert(env->lookup(Pkg(),"chain"));
	this->_MonomerId._Monomer = translate::from_object<int>::convert(env->lookup(Pkg(),"monomer"));
	return _Nil<core::T_O>();
    };
#endif


#ifdef XML_ARCHIVE
    void MonomerId_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("Chain",this->_MonomerId._Chain);
	node->attribute("Monomer",this->_MonomerId._Monomer);
    }
#endif

};
