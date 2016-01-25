#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/atom.h>
#include <cando/chem/properTorsion.h>
#include <clasp/core/wrappers.h>
namespace chem
{

// ----------------------------------------------------------------------
//



    ProperTorsion_sp ProperTorsion_O::create(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4)
    {
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
