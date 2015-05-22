#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "symbolTable.h"
#include "atom.h"
#include "properTorsion.h"
#include "core/wrappers.h"
namespace chem
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(chem,ProperTorsion_O);


    ProperTorsion_sp ProperTorsion_O::create(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4, core::Lisp_sp lisp)
    {_G();
	GC_ALLOCATE(ProperTorsion_O, obj );
	obj->setup(a1,a2,a3,a4);
	return obj;
    }


    void ProperTorsion_O::setup(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4)
    {_OF();
	this->_a1 = a1;
	this->_a2 = a2;
	this->_a3 = a3;
	this->_a4 = a4;
    }



    void ProperTorsion_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ProperTorsion_O>()
	;
    }

    void ProperTorsion_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ProperTorsion,"","",_lisp)
	;
#endif
    }



#ifdef XML_ARCHIVE
    void ProperTorsion_O::archiveBase(core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif


    void ProperTorsion_O::initialize()
    {_OF();
        this->Base::initialize();
    }



    string ProperTorsion_O::__repr__() const
    {
	stringstream ss;
	ss << this->className() << "["<<this->_a1->getName()<<"-"<<this->_a2->getName()<<"-"<<this->_a3->getName()<<"-"<<this->_a4->getName()<<"]";
	return ss.str();
    }

    

}; /* chem */
