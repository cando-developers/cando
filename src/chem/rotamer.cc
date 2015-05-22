#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "rotamer.h"
#include "core/wrappers.h"




namespace chem
{



#ifdef XML_ARCHIVE
void	RotamerAtom::archive(core::ArchiveP node)
{_G();
    node->attribute("d",this->_Distance);
    node->attribute("a",this->_Angle);
    node->attribute("t",this->_Torsion);
    LOG(BF("Serialized RotamerAtom@%X = %lf, %lf, %lf") % this % this->_Distance % this->_Angle % this->_Torsion );
};

#endif







// ----------------------------------------------------------------------
//



    EXPOSE_CLASS(chem,Rotamer_O);

    void Rotamer_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Rotamer_O>()
	;
    }

    void Rotamer_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Rotamer,"","",_lisp)
	;
#endif
    }


#ifdef XML_ARCHIVE
    void Rotamer_O::archiveBase(core::ArchiveP node)
    {
	if ( node->saving() ) {
	    uint sz = this->_RotamerAtoms.size();
	    node->attribute("size",sz);
	    int i;
	    core::ArchiveP vnode;
	    RotamerAtoms::iterator vi;
	    for ( i=0,vi=this->_RotamerAtoms.begin(); vi!=this->_RotamerAtoms.end(); i++,vi++ ) {
		vnode = node->createChildNode("V");
		vnode->attribute("i",i);
		(*vi).archive(vnode);
	    }
	} else {
	    uint sz=0, i=0;
	    node->attribute("size",sz);
	    this->_RotamerAtoms.resize(sz);
	    core::VectorArchiveNodes::iterator vi;
	    i = 0;
	    for ( vi=node->begin_Children(); vi!=node->end_Children(); vi++ ) {
		(*vi)->setRecognized(true);
		(*vi)->attribute("i",i);
		this->_RotamerAtoms[i].archive(*vi);
		i++;
	    }
	}
    }
#endif


    void Rotamer_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    

}; /* chem */
