#define	DEBUG_LEVEL_FULL       
//
// (C) 2004 Christian E. Schafmeister
//

/*
 *	ffTypes.cc
 *
 *	Maintains databases and structures to store types and type assignement
 *	rules.
 */
#include "ffTypesDb.h"
#include "loop.h"
#include "core/wrappers.h"





namespace chem {






void	FFTypesDb_O::initialize()
{
    this->Base::initialize();
    this->_TypeAssignmentRules.clear();
}



void    FFTypesDb_O::assignTypes(chem::Matter_sp matter)
{ _G();
chem::Loop    				lAtoms;
chem::Atom_sp  				atom;
chem::Matter_sp				c;
gctools::Vec0<chem::ChemInfo_sp>::iterator  it;
gctools::Vec0<chem::Matter_sp>::iterator   	mi, ri, ai;
chem::Molecule_sp                     	mol;
chem::Residue_sp				res;
string 		                        name;
    lAtoms.loopTopGoal(matter,ATOMS);
    LOG(BF("defined loop") );
    while ( lAtoms.advanceLoopAndProcess() ) {
        LOG(BF("Getting container") );
        c = lAtoms.getMatter();
//        rawGet = c.get();
//        castGet = dynamic_cast<chem::Atom_O*>(c.get());
//        LOG(BF("rawGet = %X") % rawGet );
//        LOG(BF("castGet = %X") % castGet );
//        LOG(BF("getting first atom in loop") );
        atom = downcast<chem::Atom_O>(c);
        LOG(BF("Got atom") );
        LOG(BF("atom name: %s") % atom->getName().c_str() );
        LOG(BF("Assigning type for atom: %s") % atom->description().c_str()  );
	{_BLOCK_TRACE("Testing every type rule");
	    for ( it=this->_TypeAssignmentRules.begin();
		    it!=this->_TypeAssignmentRules.end(); it++ ) 
	    {_BLOCK_TRACEF(BF("Testing rule code(%s)") % (*it)->getCode().c_str() );
//		LOG(BF("as xml: %s") % ((*it)->asXmlString().c_str() ) );
		if ( (*it)->matches(atom) ) 
		{
		    LOG(BF("Rule MATCH!!!") );
		    break;
		}
		LOG(BF("Rule does not match, keep going") );
	    }
	} 
    }
}



    void	FFTypesDb_O::archiveBase(core::ArchiveP node)
    {	
	this->FFBaseDb_O::archiveBase(node);
	node->attributeVector0(KW("type-rules"), this->_TypeAssignmentRules);
    }


chem::ChemInfo_sp	FFTypesDb_O::getRule(uint index)
{_OF();
    ASSERT_lessThan(index,this->_TypeAssignmentRules.size());
    return this->_TypeAssignmentRules[index];
}





    void FFTypesDb_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<FFTypesDb_O>()
	.def("assignTypes",&FFTypesDb_O::assignTypes)
    ;
}

    void FFTypesDb_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,FFTypesDb,"","",_lisp)
	.def("assignTypes",&FFTypesDb_O::assignTypes)
    ;
#endif
}






    EXPOSE_CLASS(chem, FFTypesDb_O);
};




