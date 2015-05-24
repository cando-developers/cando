#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/kinematics/chiDihedrals.h>
#include <clasp/core/multipleValues.h>
#include <clasp/core/wrappers.h>
namespace kinematics
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(kinematics,ChiDihedral_O);

    void ChiDihedral_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<ChiDihedral_O>()
//	    .def_raw("__init__",&ChiDihedral_O::__init__,"(self atom1Index atom2Index atom3Index atom4Index)")
	    .def("chi-dihedral-atom-indices",&ChiDihedral_O::atomIndices)
	    .def("chi-dihedral-atom1-index",&ChiDihedral_O::atom1Index)
	    .def("chi-dihedral-atom2-index",&ChiDihedral_O::atom2Index)
	    .def("chi-dihedral-atom3-index",&ChiDihedral_O::atom3Index)
	    .def("chi-dihedral-atom4-index",&ChiDihedral_O::atom4Index)
	;
    }

    void ChiDihedral_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,ChiDihedral,"","",_lisp)
//	.def_raw("__init__",&"(self atom1Index atom2Index atom3Index atom4Index)")
	    .def("chi-dihedral-atom-indices",&ChiDihedral_O::atomIndices)
	    .def("chi-dihedral-atom1-index",&ChiDihedral_O::atom1Index)
	    .def("chi-dihedral-atom2-index",&ChiDihedral_O::atom2Index)
	    .def("chi-dihedral-atom3-index",&ChiDihedral_O::atom3Index)
	    .def("chi-dihedral-atom4-index",&ChiDihedral_O::atom4Index)
	;
#endif
    }


#if INIT_TO_FACTORIES

#define ARGS_ChiDihedral_O_make "(atom1-index atom2-index atom3-index atom4-index)"
#define DECL_ChiDihedral_O_make ""
#define DOCS_ChiDihedral_O_make "make ChiDihedral"
  ChiDihedral_sp ChiDihedral_O::make(const int atom1Index, const int atom2Index, const int atom3Index, const int atom4Index)
    {_G();
	GC_ALLOCATE(ChiDihedral_O, me );
      me->_Atom1 = atom1Index;
      me->_Atom2 = atom2Index;
      me->_Atom3 = atom3Index;
      me->_Atom4 = atom4Index;
      return me;
    };

#else
    core::T_sp ChiDihedral_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {_G();
	this->_Atom1 = env->lookup(Pkg(),"atom1Index").as<core::Rational_O>()->as_int();
	this->_Atom2 = env->lookup(Pkg(),"atom2Index").as<core::Rational_O>()->as_int();
	this->_Atom3 = env->lookup(Pkg(),"atom3Index").as<core::Rational_O>()->as_int();
	this->_Atom4 = env->lookup(Pkg(),"atom4Index").as<core::Rational_O>()->as_int();
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void ChiDihedral_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
	node->attribute("a1",this->_Atom1);
	node->attribute("a2",this->_Atom2);
	node->attribute("a3",this->_Atom3);
	node->attribute("a4",this->_Atom4);
    }
#endif


    void ChiDihedral_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    core::T_mv ChiDihedral_O::atomIndices() const
    {_OF();
	return Values(
	    core::clasp_make_fixnum(this->_Atom1),
	    core::clasp_make_fixnum(this->_Atom2),
	    core::clasp_make_fixnum(this->_Atom3),
	    core::clasp_make_fixnum(this->_Atom4));
    }
    



    EXPOSE_CLASS(kinematics,ChiList_O);

    void ChiList_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<ChiList_O>()
//	.def_raw("__init__",&ChiList_O::__init__,"(self)")
	    .def("chi-list-append",&ChiList_O::append)
	;
    }

    void ChiList_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(KinPkg,ChiList,"","",_lisp)
//	.def_raw("__init__",&"(self)")
	    .def("chi-list-append",&ChiList_O::append)
	;
#endif
    }


#if INIT_TO_FACTORIES
#else
    core::T_sp ChiList_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = translate::from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _Nil<core::T_O>();
    }
#endif

#ifdef XML_ARCHIVE
    void ChiList_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->archiveVector0("chis",this->_Chis);
	// Archive other instance variables here
    }
#endif


    void ChiList_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    void ChiList_O::append(ChiDihedral_sp chi)
    {_OF();
	this->_Chis.append(chi);
    }

    

}; /* kinematics */
