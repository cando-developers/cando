#define	DEBUG_LEVEL_FULL

#include <cando/chem/virtualAtom.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <clasp/core/cons.h>
#include <cando/chem/residue.h>
#include <cando/chem/calculatePosition.h>
#include <clasp/core/wrappers.h>




namespace chem {



    VirtualAtom_sp VirtualAtom_O::create(MatterName name, CalculatePosition_sp proc)
{_G();
    LOG(BF("status") );
    GC_ALLOCATE(VirtualAtom_O, va );
    LOG(BF("status") );
    va->setName(name);
    LOG(BF("status") );
    va->_CalculatePosition = proc;
    LOG(BF("status") );
    return va;
}

//
// Copy constructor
//
VirtualAtom_O::VirtualAtom_O(const VirtualAtom_O& ss) :Atom_O(ss)
{
    this->_CalculatePosition = ss._CalculatePosition;
//    this->_Environment = ss._Environment;
}

//
//	copy
//
//	Allocate a copy of this atom and return it.
//	Keep track of the new atom so that we can redirect pointers to the copy.
//
Matter_sp VirtualAtom_O::copy()
{_G();
    LOG(BF("Copying atom @%X") % this );
    GC_COPY(VirtualAtom_O,aNew, *this); // VirtualAtom_sp aNew = RP_Copy<VirtualAtom_O>(this);
    this->copyAtom = aNew;
    LOG(BF("    copy atom== %s") % aNew->description().c_str() );
    return(aNew);
}



void	VirtualAtom_O::initialize()
{
    this->Base::initialize();
    this->_CalculatePosition = _Nil<CalculatePosition_O>();
}

    bool VirtualAtom_O::equal(core::T_sp obj) const
{
    if ( this->eq(obj) ) return true;
    if ( !obj.isA<VirtualAtom_O>() )
    {
	VirtualAtom_sp other = obj.as<VirtualAtom_O>();
	if ( other->getName() != this->getName() ) return false;
	if ( other->numberOfBonds() != this->numberOfBonds() ) return false;
	return true;
    }
    return false;
}




#ifdef XML_ARCHIVE
    void	VirtualAtom_O::archive(core::ArchiveP node)
{
    this->Base::archive(node);
    node->attributeIfNotNil("calculatePosition",this->_CalculatePosition);
//    node->attribute("environment",this->_Environment);
}
#endif



void	VirtualAtom_O::calculatePosition()
{_G();
    ASSERTNOTNULL(this->_CalculatePosition);
//    ASSERTNOTNULL(this->_Environment);
    Matter_wp wmres = this->containedBy();
    ASSERTNOTNULL(wmres);
    Matter_sp mres = wmres.lock();
    ASSERT(mres.notnilp());
    Residue_sp res = mres.as<Residue_O>();
    Vector3 vpos = this->_CalculatePosition->calculatePosition(res);
    this->setPosition(vpos);
}






    void VirtualAtom_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<VirtualAtom_O>()
	;
}

    void VirtualAtom_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON //[
	PYTHON_CLASS(ChemPkg,VirtualAtom,"","",_lisp)
    ;
#endif //]
}



    EXPOSE_CLASS(chem,VirtualAtom_O);
};




