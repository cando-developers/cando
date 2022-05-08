/*
    File: structureComparer.cc
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
       
#define	DEBUG_LEVEL_NONE


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




void	StructureComparer_O::initialize()
{
    this->Base::initialize();
    this->_AllAtoms.clear();
    this->_SuperposeAtoms.clear();
}




void	StructureComparer_O::setFixedMatter(Matter_sp matter)
{
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
{
int		ia;
this->_FixedCoordinates = geom::SimpleVectorCoordinate_O::make(this->_SuperposeAtoms.size());
gctools::SmallOrderedSet<Atom_sp>::iterator	oa;
    for ( oa=this->_SuperposeAtoms.begin(),ia=0; oa!=this->_SuperposeAtoms.end(); oa++,ia++ )
    {
      (*this->_FixedCoordinates)[ia] = (*oa)->getPosition();
    }
}


    void	StructureComparer_O::superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms)
{
    this->_SuperposeAtoms = atoms;
    this->initializeFixedCoordinates();
};

void	StructureComparer_O::superposeAllHeavyAtoms()
{
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
{
    gctools::SmallOrderedSet<Atom_sp>::iterator			ai;
SuperposeEngine_sp				superposer;
Matrix					transform;
int					i;
geom::SimpleVectorCoordinate_sp			moveableCoordinates;
    if ( matter != this->_Matter )
    {
        SIMPLE_ERROR(("The Matter passed must be the same as the one defined for the StructureComparer"));
    }

    if ( this->_SuperposeAtoms.size() < 3 )
    {
	SIMPLE_ERROR(("You must have defined at least three atoms to superpose"));
    }
    LOG("Number of superpose atoms = %d" , this->_SuperposeAtoms.size()  );


	//
    	// Now check if the structure is new or not.
	// First assemble the superposable coordinates of this conformation
	//
    moveableCoordinates= geom::SimpleVectorCoordinate_O::make(this->_SuperposeAtoms.size());
    for ( ai=this->_SuperposeAtoms.begin(), i=0;
		ai!=this->_SuperposeAtoms.end(); ai++, i++)
    {
	LOG("Extracting coordinate for superpose atom(%s)" , (*ai)->getName().c_str()  );
	(*moveableCoordinates)[i] = (*ai)->getPosition();
    }
    LOG("There are %d superposable atoms" , this->_SuperposeAtoms.size()  );

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
{
    node->attribute( "Matter", this->_Matter );
    node->archiveOrderedSet( "AllAtoms", this->_AllAtoms );
    node->archiveOrderedSet( "SuperposeAtoms", this->_SuperposeAtoms );
    node->attribute("FixedCoordinates",this->_FixedCoordinates);
}
#endif



};
