/*
    File: superposableConformationCollection.cc
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

#include <cando/chem/superposableConformationCollection.h>
//#include "core/archiveNode.h"
#include <clasp/core/lispStream.h>
#include <cando/chem/superposeEngine.h>
#include <cando/geom/color.h>
#include <clasp/core/wrappers.h>



namespace chem {
//
// Constructor
//


void	SuperposableConformationCollection_O::initialize()
{
  this->Base::initialize();
  this->_SuperposeAtomIndices = core::ComplexVector_byte32_t_O::make_vector(16,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
  this->_RmsCutOff = 0.1;
}

void SuperposableConformationCollection_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,rmsCutOff),this->_RmsCutOff);
  node->field(INTERN_(kw,atomIndices),this->_SuperposeAtomIndices);
  this->Base::fields(node);
}


//
// Destructor
//



CL_DEFMETHOD void SuperposableConformationCollection_O::clearSuperposeAtoms()
{
  ASSERTNOTNULL(this->_SuperposeAtomIndices);
  this->_SuperposeAtomIndices->fillPointerSet(0);
}



CL_DEFMETHOD void SuperposableConformationCollection_O::addSuperposeAtom(Atom_sp a)
{
  ConformationCollection_O::atomIterator	ai;
  int					idx;
  ASSERTNOTNULL(this->_SuperposeAtomIndices);
  for ( idx=0,ai=this->begin_AllAtoms(); ai!=this->end_AllAtoms(); idx++,ai++ )
  {
    if ( (*ai) == a )
    {
      this->_SuperposeAtomIndices->vectorPushExtend(idx);
      return;
    }
  }
  SIMPLE_ERROR(("You cannot add superpose atom: "+a->description()+" because it isn't in the ConformationCollection"));
}


CL_DEFMETHOD void	SuperposableConformationCollection_O::superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms)
{
  gctools::SmallOrderedSet<Atom_sp>::iterator	ai;
  this->clearSuperposeAtoms();
  for ( ai=atoms.begin(); ai!=atoms.end(); ai++ )
  {
    this->addSuperposeAtom(*ai);
  }
}

CL_DEFMETHOD void	SuperposableConformationCollection_O::superposeAllHeavyAtoms()
{
  ConformationCollection_O::atomIterator	ai;
  int					idx;
  ASSERTNOTNULL(this->_SuperposeAtomIndices);
  this->clearSuperposeAtoms();
  for ( idx=0,ai=this->begin_AllAtoms(); ai!=this->end_AllAtoms(); idx++,ai++ )
  {
    if ( (*ai)->getElement() != element_H )
    {
      this->_SuperposeAtomIndices->vectorPushExtend(idx);
    }
  }
}


CL_DEFMETHOD gctools::SmallOrderedSet<Atom_sp>	SuperposableConformationCollection_O::getSuperposeAtoms()
{
  gctools::SmallOrderedSet<Atom_sp>	result;
  Atom_sp			a;
  ASSERTNOTNULL(this->_SuperposeAtomIndices);
  for ( size_t si=0, siEnd(this->_SuperposeAtomIndices->length()); si<siEnd; si++ )
  {
    a = this->_AllAtoms[(*this->_SuperposeAtomIndices)[si]];
    result.insert(a);
  }
  return result;
}


int	SuperposableConformationCollection_O::numberOfSuperposeAtoms()
{
  return this->_SuperposeAtomIndices->length();
}



CL_DOCSTRING(R"dx(Create an entry if the conformation is new according to the rms deviation from existing conformations.  If a new entry is created - return it - otherwise return NIL.)dx");
CL_DEFMETHOD core::T_sp	SuperposableConformationCollection_O::createEntryIfConformationIsNew(Matter_sp matter)
{
  geom::SimpleVectorCoordinate_sp			newConf;
  SuperposeEngine_sp				superposer;
  Matrix					transform;
  ConformationCollectionEntry_sp		entry;
  ASSERT(matter==this->_Matter);
  LOG("Number of superpose atoms = %d" , this->_SuperposeAtomIndices->size()  );
  ASSERT(this->_SuperposeAtomIndices->length() >= 4 );


	//
    	// Now check if the structure is new or not.
	// First assemble the superposable coordinates of this conformation
	//
  newConf = this->_SimpleVectorCoordinate(matter);

  {

    LOG("Fixed points at the start of superposer: %s" , (newConf->asXmlString().c_str() ) );

    {

		//
		// Now compare this conformation to every one in the database
		//
      superposer = SuperposeEngine_O::create();
      superposer->setMoveablePoints(this->_SuperposeAtomIndices,newConf);
      double rms;
      gctools::Vec0<ConformationCollectionEntry_sp>::iterator	ci;
      geom::SimpleVectorCoordinate_sp				entryCoordinates;

		//
		// Loop through the low-high entry range and compare the structures
		// If this one is the same as any of those return that we've seen this structure before
		//
      for ( ci=this->_Entries.begin(); ci!=this->_Entries.end(); ci++ )
      {
        entryCoordinates = (*ci)->getAllCoordinates();
        LOG("Entry coordinates before superpose:%s" , (entryCoordinates->asXmlString()) );
        superposer->setFixedPoints(this->_SuperposeAtomIndices, entryCoordinates);
        LOG("About to superpose" );
        transform = superposer->superpose();
        rms = superposer->rootMeanSquareDifference();
#ifdef DEBUG_ON
        int					i=0;
        LOG("The rms difference(%lf) with structure(%d) " "compared to _RmsCutOff(%lf)" , rms , i , this->_RmsCutOff  );
#endif
        if ( rms < this->_RmsCutOff )
        {
          LOG("Found an identical minimum with rms(%lf)" , rms );
          LOG("Entry will not be added" );
          LOG("Fixed points at the end of superposer: %s" , (newConf->asXmlString().c_str() ) );
          return nil<core::T_O>();
        }
        LOG("EntryCoordinates points after superpose:%s" , (entryCoordinates->asXmlString().c_str() ) );
      }
      LOG("Fixed points at the end of superposer: %s" , (newConf->asXmlString().c_str() ) );
    }
  }

    	//
	// Ok, this is a new structure, so insert it into the list
	//
  entry = this->createEntry();
  geom__in_place_transform(newConf,transform);
  entry->setAllCoordinates(newConf);
  this->_Entries.push_back(entry);
  return entry;
}

#ifdef RENDER
geom::Render_sp	SuperposableConformationCollection_O::rendered(core::List_sp opts)
{
  IMPLEMENT_ME();
#if 0
  geom::FrameList_sp	frames;
  geom::DisplayList_sp	parts, dlSuperpose,dlAll;
  GrPickableMatter_sp		dlMatter;
  geom::Color_sp				color;
  Matter_sp			matter;
  ConformationCollection_O::entryIterator	si;
  ConformationCollectionEntry_sp	entry;
  geom::SimpleVectorCoordinate_sp		superposeCoords;
  geom::SimpleVectorCoordinate_O::iterator	ci;
  SuperposableConformationCollection_O::superposeAtomIndexIterator	ii;
  frames = geom::FrameList_O::create();
  matter = this->getMatter();
  for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
  {
    (*si)->writeCoordinatesToMatter(matter);
    dlAll = geom::DisplayList_O::create();
    dlMatter = GrPickableMatter_O::create();
    dlMatter->setName(_lisp->internKeyword("minimized"));
    dlMatter->setFromMatter(matter);
    dlAll->add(dlMatter);
    dlSuperpose = geom::DisplayList_O::create();
    dlSuperpose->setName(_lisp->internKeyword("superposeAtoms"));
    dlSuperpose->add(O_GrColor::systemColor(_kw_yellow));
    RPGrSphereList sphereList = O_GrSphereList::create();
    sphereList->setRadius(0.2);
    for ( ii=this->_SuperposeAtomIndices->begin(); ii!=this->_SuperposeAtomIndices->end(); ii++ )
    {
      sphereList->addVertex(this->_AllAtoms[*ii]->getPosition());
    }
    dlSuperpose->add(sphereList);
    dlAll->add(dlSuperpose);
    frames->add(dlAll);
  }
  return frames;
#endif
}
#endif



CL_DEFMETHOD
void	SuperposableConformationCollection_O::setEntryCoordinatesAsFixedWithinSuperposeEngine(
                                                                                              ConformationCollectionEntry_sp entry,
                                                                                              SuperposeEngine_sp superposer )
{
  superposer->setFixedPoints(this->_SuperposeAtomIndices,
                             entry->getAllCoordinates());
}

CL_DEFMETHOD
void	SuperposableConformationCollection_O::setEntryCoordinatesAsMoveableWithinSuperposeEngine(
                                                                                                 ConformationCollectionEntry_sp entry,
                                                                                                 SuperposeEngine_sp superposer )
{
  superposer->setMoveablePoints(this->_SuperposeAtomIndices, entry->getAllCoordinates());
}




};
