#define DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/chem/atomId.h>
#include <clasp/core/wrappers.h>



namespace chem
{

    string AtomId::asString() const
    {
	stringstream ss;
	ss << "Mol["<<this->_Molecule;
	ss << "]Res["<<this->_Residue;
	ss << "]Atom["<<this->_Atom<<"]";
	return ss.str();
    }



    EXPOSE_CLASS(chem,AtomId_O);

    void AtomId_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<AtomId_O>()
////	    .def_raw("core:__init__",&AtomId_O::__init__,"(self &key (molecule 0) (residue 0) (atom 0))")
	    .def("asString",&AtomId_O::asString,"","","",false)
	    ;
    }

    void AtomId_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,AtomId,"","",_lisp)
//	    .def("asString",&AtomId_O::asString,"","","",false)
	    ;
#endif
    }



#if INIT_TO_FACTORIES

#define ARGS_AtomId_O_make "(molecule residue atom)"
#define DECL_AtomId_O_make ""
#define DOCS_AtomId_O_make "make AtomId args: molecule residue atom"
    AtomId_sp AtomId_O::make(int mol, int res, int atom)
    {_G();
	GC_ALLOCATE(AtomId_O, ai );
	ai->_AtomId._Molecule = mol;
	ai->_AtomId._Residue = res;
	ai->_AtomId._Atom = atom;
	return ai;
    }
#endif

#ifdef XML_ARCHIVE
    void AtomId_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("Molecule",this->_AtomId._Molecule);
	node->attribute("Residue",this->_AtomId._Residue);
	node->attribute("Atom",this->_AtomId._Atom);
    }
#endif


};
