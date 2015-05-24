#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/chem/atom.h>
#include <cando/chem/improperTorsion.h>
#include <clasp/core/wrappers.h>
namespace chem
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(chem,ImproperTorsion_O);

    void ImproperTorsion_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ImproperTorsion_O>()
	;
    }

    void ImproperTorsion_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ImproperTorsion,"","",_lisp)
	;
#endif
    }






    ImproperTorsion_sp ImproperTorsion_O::create(Atom_sp a1, Atom_sp a2, Atom_sp a3central, Atom_sp a4, core::Lisp_sp lisp)
    {_G();
	GC_ALLOCATE(ImproperTorsion_O, obj );
	obj->setup(a1,a2,a3central,a4);
	return obj;
    }





    void ImproperTorsion_O::setup(Atom_sp a1, Atom_sp a2, Atom_sp a3central, Atom_sp a4)
    {_OF();
	this->_a1 = a1;
	this->_a2 = a2;
	this->_a3central = a3central;
	this->_a4 = a4;
    }








#ifdef XML_ARCHIVE
    void ImproperTorsion_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void ImproperTorsion_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    string ImproperTorsion_O::__repr__() const
    {
	stringstream ss;
	ss << this->className()
	   << "["
	   <<this->_a1->getName()
	   <<"-"<<this->_a2->getName()
	   <<"-center["
	   <<this->_a3central->getName()
	   <<"]-"
	   <<this->_a4->getName()
	   <<"]";
	return ss.str();
    }


    

}; /* chem */
