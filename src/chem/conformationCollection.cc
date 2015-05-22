       
#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	Maintain a database of conformational searches
 *
 */

#include <cando/chem/conformationCollection.h>
#include <cando/chem/loop.h>
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include <cando/chem/superposeEngine.h>
#include <cando/chem/conformationExplorer.h>
#include <clasp/core/wrappers.h>




namespace chem {


    REGISTER_CLASS(chem, ConformationCollectionEntry_O );
    REGISTER_CLASS(chem, ConformationCollection_O );


void	ConformationCollectionEntry_O::initialize()
{
    this->Base::initialize();
    this->_Data = core::HashTableEq_O::create_default();
    this->_AllCoordinates = candoBase::CoordinateArray_O::create();
}



#ifdef XML_ARCHIVE
    void	ConformationCollectionEntry_O::archive(core::ArchiveP node)
{_G();
    node->archiveWeakPointer( "ConformationCollection",
   		 this->_WeakConformationCollection );
    node->attribute("AllCoordinates", this->_AllCoordinates );
    node->attribute("Data",this->_Data);
    if ( node->loading() )
    {
        LOG(BF("Serialized object loaded from file") );
	LOG(BF("allCoords =%s") % (this->_AllCoordinates->asXmlString().c_str() ) );
    } else
    {
        LOG(BF("Serialized object output to file") );
	LOG(BF("allCoords =%s") % (this->_AllCoordinates->asXmlString().c_str() ) );
    }
}
#endif




void	ConformationCollectionEntry_O::setConformationCollection(ConformationCollection_sp s)
{_G();
candoBase::CoordinateArray_sp	ca;
    this->_WeakConformationCollection = s;
    ca = candoBase::CoordinateArray_O::create(s->numberOfAllAtoms());
    this->_AllCoordinates = ca;
}


void	ConformationCollectionEntry_O::setAllCoordinates(candoBase::CoordinateArray_sp ac)
{_G();
    ASSERTNOTNULL(ac);
    LOG(BF("setAllCoordinates:%s") % (ac->asXmlString().c_str() ) );
    LOG(BF("The address of the candoBase::CoordinateArray_sp is in o") );
    this->_AllCoordinates = ac;
}



void	ConformationCollectionEntry_O::writeCoordinatesToMatter(Matter_sp agg)
{_G();
ConformationCollection_sp			sl;
vector<Vector3>::iterator	ci;
gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
    sl = this->getConformationCollection();
    for ( ai=sl->begin_AllAtoms(),ci=this->_AllCoordinates->begin(); ai!=sl->end_AllAtoms(); ai++, ci++ )
    {
        (*ai)->setPosition(*ci);
    }
}


void	ConformationCollectionEntry_O::extractCoordinatesFromMatter(Matter_sp matter)
{_G();
ConformationCollection_sp	sl;
vector<Vector3>::iterator	ci;
gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
#ifdef	DEBUG_ConformationCollectionEntry
    this->_Status->addMessage("extractCoordinatesFromMatter");
#endif
    LOG(BF("About to get ConformationCollection") );
    sl = this->getConformationCollection();
    LOG(BF("Got ConformationCollection") );
    ASSERTNOTNULL(sl);
    candoBase::CoordinateArray_sp coords = sl->_extractCoordinateArray(matter);
    ASSERTNOTNULL(coords);
    this->_AllCoordinates = coords;
}


void	ConformationCollection_O::initialize()
{
    this->Base::initialize();
    this->_Matter = _Nil<Matter_O>();
    this->_AllAtoms.clear();
    this->_Data = core::HashTableEq_O::create_default();
    this->clearEntries();
}


void	ConformationCollection_O::clearEntries()
{_G();
    this->_Entries.clear();
}



candoBase::CoordinateArray_sp ConformationCollection_O::_extractCoordinateArray(Matter_sp agg)
{_G();
    vector<Vector3>::iterator	ci;
    gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
#ifdef	DEBUG_ConformationCollectionEntry
    this->_Status->addMessage("extractCoordinatesFromMatter");
#endif
    LOG(BF("About to get ConformationCollection") );
    LOG(BF("About to iterate over atoms and get positions") );
    candoBase::CoordinateArray_sp coords = candoBase::CoordinateArray_O::create(this->numberOfAllAtoms());
    ASSERTNOTNULL(coords);
    for ( ai=this->begin_AllAtoms(),ci=coords->begin(); ai!=this->end_AllAtoms(); ai++, ci++ )
    {
        LOG(BF("Getting position of %s") % (*ai)->description().c_str()  );
        *ci = (*ai)->getPosition();
    }
    LOG(BF("Done") );
    return coords;
}


void	ConformationCollection_O::setMatter(Matter_sp matter)
{_G();
Loop	loop;
Atom_sp	a;
    this->_Matter = matter;
    this->_AllAtoms.clear();
    loop.loopTopGoal(matter,ATOMS);
    while ( loop.advanceLoopAndProcess() )
    {
	a = loop.getAtom();
	this->_AllAtoms.insert(a);
    }
}




uint	ConformationCollection_O::addEntry(ConformationCollectionEntry_sp entry)
{_G();
    this->_Entries.append(entry);
    return 1;
}




void	ConformationCollection_O::addConformationCollection(ConformationCollection_sp orig)
{_G();
Matter_sp			matter;
ConformationCollection_O::entryIterator	si;
ConformationCollectionEntry_sp		entry;
    matter = orig->getMatter();
    ASSERTP(matter==this->getMatter(),"The Matter objects must match");
    for ( si=orig->begin_Entries(); si!=orig->end_Entries(); si++ )
    {
        (*si)->writeCoordinatesToMatter(matter);
	entry = this->createEntry();
	entry->extractCoordinatesFromMatter(matter);
        this->addEntry(entry);
    }
}



void	ConformationCollection_O::addConformationExplorerLatestConformations(ConformationExplorer_sp orig)
{_G();
Matter_sp			matter;
ConformationExplorer_O::entryIterator	si;
ConformationExplorerEntryStage_sp		entryStage;
ConformationCollectionEntry_sp			entry;
    matter = orig->getMatter();
    ASSERTP(matter==this->getMatter(),"The Matter objects must match");
    for ( si=orig->begin_Entries(); si!=orig->end_Entries(); si++ )
    {
    	entryStage = (*si)->getLastCompleteEntryStage();
        entryStage->writeCoordinatesToMatter(matter);
	entry = this->createEntry();
	entry->extractCoordinatesFromMatter(matter);
        this->addEntry(entry);
    }
}


void	ConformationCollection_O::addConformationExplorerSelectedStageConformations(ConformationExplorer_sp orig)
{_G();
Matter_sp			matter;
ConformationExplorer_O::entryIterator	si;
ConformationExplorerEntryStage_sp		entryStage;
ConformationCollectionEntry_sp			entry;
    matter = orig->getMatter();
    ASSERTP(matter==this->getMatter(),"The Matter objects must match");
    for ( si=orig->begin_Entries(); si!=orig->end_Entries(); si++ )
    {
        ASSERTP((*si)->hasSelectedStage(),"Entry missing selected stage");
    	entryStage = (*si)->getSelectedStage();
        entryStage->writeCoordinatesToMatter(matter);
	entry = this->createEntry();
	entry->extractCoordinatesFromMatter(matter);
        this->addEntry(entry);
    }
}









#ifdef XML_ARCHIVE
    void	ConformationCollection_O::archiveBase(core::ArchiveP node)
{_G();
    node->attribute( "Matter", this->_Matter );
    node->archiveOrderedSet( "AllAtoms", this->_AllAtoms );
    node->archiveVector0( "Entries", this->_Entries );
    node->attributeIfNotNil( "Data", this->_Data );
}
#endif

#ifdef XML_ARCHIVE
    void	ConformationCollection_O::archive(core::ArchiveP node)
{_G();
    this->archiveBase(node);
}
#endif




#ifdef RENDER
candoBase::Render_sp	ConformationCollection_O::rendered(core::Cons_sp opts)
{_G();
candoBase::FrameList_sp	frames;
candoBase::DisplayList_sp	parts, dlSuperpose;
GrPickableMatter_sp		dlMatter;
candoBase::Color_sp				color;
Matter_sp			matter;
ConformationCollection_O::entryIterator	si;
ConformationCollectionEntry_sp	entry;
candoBase::CoordinateArray_sp		superposeCoords;
candoBase::CoordinateArray_O::iterator	ci;
frames = candoBase::FrameList_O::create();
    matter = this->getMatter();
    for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
    {
        (*si)->writeCoordinatesToMatter(matter);
	dlMatter = GrPickableMatter_O::create();
	FIX_ME(); dlMatter->setName(_lisp->internKeyword("gaussianMinimized"));
	dlMatter->setFromMatter(matter);
	frames->append(dlMatter);
    }
    return frames;
}
#endif








#ifdef XML_ARCHIVE
void	ConformationCollection_O::saveAs(const string& fn)
{_G();
    core::XmlSaveArchive_sp	xml;

    xml = core::XmlSaveArchive_O::create();
    xml->put("ConformationCollection",this->sharedThis<ConformationCollection_O>());
    xml->saveAs(fn);
}
#endif

uint	ConformationCollection_O::getEntryIndex(ConformationCollectionEntry_sp entry)
{_G();
uint		i;
entryIterator	ei;
    for ( i=0,ei=this->_Entries.begin(); ei!=this->_Entries.end(); ei++,i++ )
    {
        if ( (*ei) == entry )
	{
	    return i;
	}
    }
    SIMPLE_ERROR(BF("Could not find entry in ConformationCollection"));
}




ConformationCollectionEntry_sp	ConformationCollection_O::getEntry(uint i)
{_G();
    ASSERT_lessThan(i,this->_Entries.size());
    ConformationCollectionEntry_sp e;
    e = this->_Entries[i];
    return e;
}

ConformationCollectionEntry_sp	ConformationCollection_O::firstEntry()
{_G();
    ASSERT(this->_Entries.size() > 0 );
    return this->_Entries[0];
}


ConformationCollectionEntry_sp	ConformationCollection_O::createEntry()
{_G();
    GC_ALLOCATE(ConformationCollectionEntry_O, entry );
    entry->setConformationCollection(this->sharedThis<ConformationCollection_O>());
    return entry;
}




};

