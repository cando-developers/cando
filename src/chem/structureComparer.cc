       
#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	Maintain a database of conformational searches
 *
 */

#include <cando/chem/structureComparer.h>
#include <cando/chem/loop.h>
#include <cando/chem/superposeEngine.h>
#include <clasp/core/wrappers.h>



namespace chem {

    REGISTER_CLASS(chem, StructureComparer_O );



void	StructureComparer_O::initialize()
{
    this->Base::initialize();
    this->_AllAtoms.clear();
    this->_SuperposeAtoms.clear();
}




void	StructureComparer_O::setFixedMatter(Matter_sp matter)
{_G();
Loop	loop;
Atom_sp	a;
    this->_Matter = matter;
    this->_AllAtoms.clear();
    loop.loopTopGoal(matter,ATOMS);
    while ( loop.advanceLoopAndProcess() ) {
	a = loop.getAtom();
	this->_AllAtoms.insert(a);
    }
    this->_SuperposeAtoms.clear();
}


void	StructureComparer_O::initializeFixedCoordinates()
{_G();
int		ia;
this->_FixedCoordinates = geom::CoordinateArray_O::create(this->_SuperposeAtoms.size());
gctools::SmallOrderedSet<Atom_sp>::iterator	oa;
    for ( oa=this->_SuperposeAtoms.begin(),ia=0; oa!=this->_SuperposeAtoms.end(); oa++,ia++ )
    {
        this->_FixedCoordinates->setElement(ia,(*oa)->getPosition());
    }
}


    void	StructureComparer_O::superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms)
{_G();
    this->_SuperposeAtoms = atoms;
    this->initializeFixedCoordinates();
};

void	StructureComparer_O::superposeAllHeavyAtoms()
{_G();
    gctools::SmallOrderedSet<Atom_sp>::iterator	oa;
    this->_SuperposeAtoms.clear();
    for ( oa=this->_AllAtoms.begin(); oa!=this->_AllAtoms.end(); oa++ )
    {
        if ( (*oa)->getElement() != element_H )
	{
	    this->_SuperposeAtoms.insert(*oa);
	}
    }
    this->initializeFixedCoordinates();
}





double StructureComparer_O::calculateRmsWithMatter(Matter_sp matter)
{_G();
    gctools::SmallOrderedSet<Atom_sp>::iterator			ai;
SuperposeEngine_sp				superposer;
Matrix					transform;
int					i;
geom::CoordinateArray_sp			moveableCoordinates;
    if ( matter != this->_Matter )
    {
        SIMPLE_ERROR(BF("The Matter passed must be the same as the one defined for the StructureComparer"));
    }

    if ( this->_SuperposeAtoms.size() < 3 )
    {
	SIMPLE_ERROR(BF("You must have defined at least three atoms to superpose"));
    }
    LOG(BF("Number of superpose atoms = %d") % this->_SuperposeAtoms.size()  );


	//
    	// Now check if the structure is new or not.
	// First assemble the superposable coordinates of this conformation
	//
    moveableCoordinates= geom::CoordinateArray_O::create(this->_SuperposeAtoms.size());
    for ( ai=this->_SuperposeAtoms.begin(), i=0;
		ai!=this->_SuperposeAtoms.end(); ai++, i++)
    {
	LOG(BF("Extracting coordinate for superpose atom(%s)") % (*ai)->getName().c_str()  );
	moveableCoordinates->setElement(i,(*ai)->getPosition());
    }
    LOG(BF("There are %d superposable atoms") % this->_SuperposeAtoms.size()  );

	//
	// Now compare this conformation to every one in the database
	//
    superposer = SuperposeEngine_O::create();
    superposer->setFixedAllPoints(this->_FixedCoordinates);
    superposer->setMoveableAllPoints(moveableCoordinates);
    double rms;
    this->_TransformToFixed = superposer->superpose();
    rms = superposer->rootMeanSquareDifference();
    return rms;
}



#ifdef XML_ARCHIVE
    void	StructureComparer_O::archive(core::ArchiveP node)
{_G();
    node->attribute( "Matter", this->_Matter );
    node->archiveOrderedSet( "AllAtoms", this->_AllAtoms );
    node->archiveOrderedSet( "SuperposeAtoms", this->_SuperposeAtoms );
    node->attribute("FixedCoordinates",this->_FixedCoordinates);
}
#endif



};
