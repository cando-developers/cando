       
#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	Maintain a database of conformational searches
 *
 */

#include <cando/chem/structureList.h>
#include <cando/chem/loop.h>
#include <cando/chem/superposeEngine.h>
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include <clasp/core/hashTableEq.h>
#include <cando/geom/color.h>
#include <clasp/core/wrappers.h>




namespace chem {


Structure_Old_ListEntry_sp	Structure_Old_ListEntry_O::create(Structure_Old_List_sp s)
{
    IMPLEMENT_ME();
}



void	Structure_Old_ListEntry_O::initialize()
{_OF();
    this->Base::initialize();
    LOG(BF("Created StructureListEntry") );
    this->_Members = 0;
    this->_AllCoordinates = geom::CoordinateArray_O::create();
    this->_SuperposeCoordinates = geom::CoordinateArray_O::create();
    this->_Data = core::HashTableEq_O::create_default();
}



#ifdef XML_ARCHIVE
    void	Structure_Old_ListEntry_O::archive(core::ArchiveP node)
{
    node->archiveWeakPointer( "structureList", this->_WeakStructureList );
    node->attribute("allCoords", this->_AllCoordinates );
    node->attribute("superposeCoords", this->_SuperposeCoordinates );
    node->attribute("data",this->_Data);
    if ( node->loading() )
    {
        LOG(BF("Serialized object loaded from file") );
	LOG(BF("allCoords =%s") % (this->_AllCoordinates->asXmlString().c_str() ) );
	LOG(BF("superposeCoords =%s") % (this->_SuperposeCoordinates->asXmlString().c_str() ) );
    } else
    {
        LOG(BF("Serialized object output to file") );
	LOG(BF("allCoords =%s") % (this->_AllCoordinates->asXmlString().c_str() ) );
	LOG(BF("superposeCoords =%s") % (this->_SuperposeCoordinates->asXmlString().c_str() ) );
    }
}
#endif

void	Structure_Old_ListEntry_O::setAllCoordinates(geom::CoordinateArray_sp ac)
{
    ASSERTNOTNULL(ac);
    LOG(BF("setAllCoordinates:%s") % (ac->asXmlString().c_str() ) );
    LOG(BF("The address of the geom::CoordinateArray_sp is in o") );
    this->_AllCoordinates = ac;
}

void	Structure_Old_ListEntry_O::setSuperposeCoordinates(geom::CoordinateArray_sp ac)
{
    LOG(BF("setSuperposeCoordinates:%s") % (ac->asXmlString().c_str() ) );
    this->_SuperposeCoordinates = ac;
}



double	Structure_Old_ListEntry_O::rmsCompareToEntry(Structure_Old_ListEntry_sp other)
{
   IMPLEMENT_ME();
}





void	Structure_Old_List_O::initialize()
{
    this->Base::initialize();
    this->_RmsCutOff = 0.05;
    this->_Matter = _Nil<Matter_O>();
    this->_SuperposeAtoms.clear();
    this->_AllAtoms.clear();
    this->_Data = core::HashTableEq_O::create_default();
    this->clearEntries();
}


void	Structure_Old_ListEntry_O::writeCoordinatesToMatter(Matter_sp agg)
{
Structure_Old_List_sp			sl;
vector<Vector3>::iterator	ci;
gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
    sl = this->getStructureList();
    for ( ai=sl->begin_AllAtoms(),ci=this->_AllCoordinates->begin(); ai!=sl->end_AllAtoms(); ai++, ci++ )
    {
        (*ai)->setPosition(*ci);
    }
}


    void	Structure_Old_ListEntry_O::translateAllCoordinates( const Vector3& offset)
{
Structure_Old_List_sp			sl;
vector<Vector3>::iterator	ci;
    for ( ci=this->_AllCoordinates->begin(); ci!=this->_AllCoordinates->end(); ci++ )
    {
        *ci = (*ci).add(offset);
    }
    for ( ci=this->_SuperposeCoordinates->begin(); ci!=this->_SuperposeCoordinates->end(); ci++ )
    {
        *ci = (*ci).add(offset);
    }
}



void	Structure_Old_ListEntry_O::extractCoordinatesFromMatter(Matter_sp agg)
{
Structure_Old_List_sp			sl;
vector<Vector3>::iterator	ci;
gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
#ifdef	DEBUG_StructureListEntry
    this->_Status->addMessage("extractCoordinatesFromMatter");
#endif
    sl = this->getStructureList();
    for ( ai=sl->begin_AllAtoms(),ci=this->_AllCoordinates->begin(); ai!=sl->end_AllAtoms(); ai++, ci++ )
    {
        *ci = (*ai)->getPosition();
    }
    for ( ai=sl->begin_SuperposeAtoms(),ci=this->_SuperposeCoordinates->begin(); ai!=sl->end_SuperposeAtoms(); ai++, ci++ )
    {
        *ci = (*ai)->getPosition();
    }
}

void	Structure_Old_List_O::clearEntries()
{
    this->_Entries.clear();
}




void	Structure_Old_List_O::setMatter(Matter_sp matter)
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

/*!
 * Return True if the structure of the matter
 * has not been seen in the database.
 */
Structure_Old_ListEntry_sp Structure_Old_List_O::createStructureListEntryIfConformationIsNew(Matter_sp matter)
{
    gctools::SmallOrderedSet<Atom_sp>::iterator			ai;
    gctools::SmallOrderedSet<Atom_sp>::iterator			lai;
geom::CoordinateArray_sp			newConf;
SuperposeEngine_sp				superposer;
Matrix					transform;
uint					i;
Structure_Old_ListEntry_sp			entry;
    if ( matter != this->_Matter )
    {
        SIMPLE_ERROR(BF("The Matter passed must be the same as the one defined for the StructureList"));
    }

    if ( this->_SuperposeAtoms.size() < 4 )
    {
	SIMPLE_ERROR(BF("You must have defined at least three atoms to superpose"));
    }
    LOG(BF("Number of superpose atoms = %d") % this->_SuperposeAtoms.size()  );


	//
    	// Now check if the structure is new or not.
	// First assemble the superposable coordinates of this conformation
	//
    newConf = geom::CoordinateArray_O::create(this->_SuperposeAtoms.size());
    for ( ai=this->_SuperposeAtoms.begin(), i=0;
		ai!=this->_SuperposeAtoms.end(); ai++, i++)
    {
	LOG(BF("Extracting coordinate for superpose atom(%s)") % (*ai)->getName().c_str()  );
	newConf->setElement(i,(*ai)->getPosition());
    }
    LOG(BF("There are %d superposable atoms") % this->_SuperposeAtoms.size()  );

    {_BLOCK_TRACE("Comparing structure to known structures using superposer");

	LOG(BF("Fixed points at the start of superposer: %s") % (newConf->asXmlString().c_str() ) );

	{_BLOCK_TRACE("Creating and evaluating superposer");

		//
		// Now compare this conformation to every one in the database
		//
	    superposer = SuperposeEngine_O::create();
	    superposer->setFixedAllPoints(newConf);
	    double rms;
	    geom::CoordinateArray_sp				moveable;

		//
		// Loop through the low-high entry range and compare the structures
		// If this one is the same as any of those return that we've seen this structure before
		//
	    for ( auto ci=this->_Entries.begin(); ci!=this->_Entries.end(); ci++ )
	    {
	        moveable = (*ci)->getSuperposeCoordinates();
		LOG(BF("Moveable points before superpose:%s") % (moveable->asXmlString().c_str() ) );
		superposer->setMoveableAllPoints(moveable);
		transform = superposer->superpose();
		rms = superposer->rootMeanSquareDifference();
		LOG(BF("The rms difference(%lf) with structure(%d) " "compared to _RmsCutOff(%lf)") % rms % i % this->_RmsCutOff  );
		if ( rms < this->_RmsCutOff )
		{
		    LOG(BF("Found an identical minimum with rms(%lf)") % rms );
		    LOG(BF("Entry will not be added") );
		    (*ci)->setMembers((*ci)->getMembers()+1);
		    Structure_Old_ListEntry_sp zilch;
		    zilch = _Nil<Structure_Old_ListEntry_O>();
	    	    LOG(BF("Fixed points at the end of superposer: %s") % (newConf->asXmlString().c_str() ) );
		    return zilch;
		}
		LOG(BF("Moveable points after superpose:%s") % (moveable->asXmlString().c_str() ) );
	    }
	    LOG(BF("Fixed points at the end of superposer: %s") % (newConf->asXmlString().c_str() ) );
	}
    }

    	//
	// Ok, this is a new structure, so insert it into the list
	//
    entry = Structure_Old_ListEntry_O::create(this->sharedThis<Structure_Old_List_O>());
    entry->setMembers(1);
		//
		// Put all of the atom coordinates into a coordinate array
		// and put that into the database
		//
    newConf = geom::CoordinateArray_O::create(this->_AllAtoms.size());
    for ( lai=this->_AllAtoms.begin(),i=0;
		lai!=this->_AllAtoms.end(); lai++, i++)
    {
	newConf->setElement(i,(*lai)->getPosition());
    }
    entry->setAllCoordinates(newConf);
		//
		// Put all of the superpose coordinates into a coordinate array
		// and put that into the database
		//
    newConf = geom::CoordinateArray_O::create(this->_SuperposeAtoms.size());
    for ( ai=this->_SuperposeAtoms.begin(),i=0;
		ai!=this->_SuperposeAtoms.end(); ai++, i++)
    {
	newConf->setElement(i,(*ai)->getPosition());
    }
    LOG(BF("Writing superpose coordinates to entry:%s") % (newConf->asXmlString().c_str() ) );
    entry->setSuperposeCoordinates(newConf);
    return entry;
}

uint Structure_Old_List_O::addEntry(Structure_Old_ListEntry_sp entry)
{
    this->_Entries.push_back(entry);
    return 1;
}



void	Structure_Old_List_O::addAllEntries(Structure_Old_List_sp orig)
{_OF();
Matter_sp			matter;
Structure_Old_List_O::entryIterator	si;
Structure_Old_ListEntry_sp		entry;
    matter = orig->getMatter();
    ASSERTP(matter==this->getMatter(),"The Matter objects must match");
    for ( si=orig->begin_Entries(); si!=orig->end_Entries(); si++ )
    {
        (*si)->writeCoordinatesToMatter(matter);
	entry = this->createStructureListEntryIfConformationIsNew(matter);
	ANN(entry);
	if ( entry.notnilp() )
	{
            this->addEntry(entry);
	}
    }
}







    void	Structure_Old_List_O::superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms)
{
    this->_SuperposeAtoms = atoms;
};


void    Structure_Old_List_O::addSuperposeAtom(Atom_sp a)
{
    this->_SuperposeAtoms.insert(a);
}



void	Structure_Old_List_O::superposeAllHeavyAtoms()
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
}




#ifdef XML_ARCHIVE
void	Structure_Old_List_O::archiveBase(core::ArchiveP node)
{
    node->attribute( "Matter", this->_Matter );
    node->attributeIfNotNil( "Data", this->_Data );
    node->archiveVector0( "Entries", this->_Entries );
    node->archiveOrderedSet( "AllAtoms", this->_AllAtoms );
    node->archiveOrderedSet( "SuperposeAtoms", this->_SuperposeAtoms );
    node->attribute( "RmsCutOff", this->_RmsCutOff );
}
#endif

#ifdef XML_ARCHIVE
void	Structure_Old_List_O::archive(core::ArchiveP node)
{
    this->archiveBase(node);
}
#endif


#if 0
Structure_Old_List_sp	Structure_Old_List_O::open_StructureList(const string& fn)
{
Structure_Old_List_sp	res;
core::XmlLoadArchive_sp	xml;

xml = core::XmlLoadArchive_O::create();
    xml->parse(fn);
    res = downcast<Structure_Old_List_O>(xml->get("structureList"));
    return res;
}
#endif


#ifdef XML_ARCHIVE
void	Structure_Old_List_O::saveAs(const string& fn)
{
    core::XmlSaveArchive_sp	xml;

    xml = core::XmlSaveArchive_O::create();
    xml->put("structureList",this->sharedThis<Structure_Old_List_O>());
    xml->saveAs(fn);
}
#endif

#ifdef RENDER
geom::Render_sp	Structure_Old_List_O::rendered(core::List_sp opts)
{
    IMPLEMENT_ME();
#if 0
geom::FrameList_sp	frames;
geom::DisplayList_sp	parts, dlSuperpose;
GrPickableMatter_sp		dlMatter;
RPGrColor			grColor;
RPGrSphereList			grSpheres;
Matter_sp			matter;
Structure_Old_List_O::entryIterator	si;
Structure_Old_ListEntry_sp		entry;
geom::CoordinateArray_sp		superposeCoords;
geom::CoordinateArray_O::iterator	ci;
frames = geom::FrameList_O::create();
    matter = this->getMatter();
    if ( this->_Entries.size() == 0 )
    {
	dlMatter = GrPickableMatter_O::create();
	dlMatter->setName(_lisp->internKeyword("molecule"));
	dlMatter->setFromMatter(matter);
	frames->add(dlMatter);
    } else
    {
	for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
	{
	    parts = geom::DisplayList_O::create();
	    (*si)->writeCoordinatesToMatter(matter);
	    dlMatter = GrPickableMatter_O::create();
	    dlMatter->setName(_lisp->internKeyword("molecule"));
	    dlMatter->setFromMatter(matter);
	    parts->append(dlMatter);
	    dlSuperpose = geom::DisplayList_O::create();
	    dlSuperpose->setName(_lisp->internKeyword("superpose"));
	    grColor = O_GrColor::systemColor(_kw_green);
	    dlSuperpose->append(grColor);
	    superposeCoords = (*si)->getSuperposeCoordinates();
	    grSpheres = O_GrSphereList::create();
	    grSpheres->setRadius(0.2);
	    for ( ci=superposeCoords->begin(); ci!=superposeCoords->end(); ci++ )
	    {
		grSpheres->addVertex(*ci);
	    }
	    dlSuperpose->append(grSpheres);
	    parts->append(dlSuperpose);
	    frames->add(parts);
	}
    }
    return frames;
#endif
}
#endif





Structure_Old_ListEntry_sp	Structure_Old_List_O::getEntry(uint i)
{_OF();
    ASSERT_lessThan(i,this->_Entries.size());
    Structure_Old_ListEntry_sp e;
    e = this->_Entries[i];
    return e;
}

};

