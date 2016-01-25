       
#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	Maintain a database of conformational searches
 *
 */

#include "core/common.h"
#include "core/numbers.h"
#include "confSearchDatabase.h"
#include "loop.h"
#include "superposeEngine.h"
#include "core/wrappers.h"



namespace chem {




void	ConformationalSearchDatabase_O::initialize()
{
    this->Base::initialize();
    this->_EnergyCutOff = 5.0; // kCal/mole
    this->setRmsCutOff(0.1);
    this->enableDisplayProgress();
}


void	ConformationalSearchDatabase_O::clearEntries()
{
    this->SuperposableConformationCollection_O::clearEntries();
}




/*!
 * Return True if the structure is added to the database, False otherwise
 */
uint ConformationalSearchDatabase_O::addEntry(ConformationCollectionEntry_sp entry)
{
chem::ConformationCollectionEntry_sp	firstEntry;
double				lowestEnergy, newEnergy;

    if ( this->numberOfEntries() == 0 )
    {
        this->SuperposableConformationCollection_O::addEntry(entry);
	return success;
    }

    firstEntry = this->getEntry(0);
    lowestEnergy = firstEntry->getData()->getDouble("energy");
    newEnergy = entry->getData()->getDouble("energy");
    if ( newEnergy-lowestEnergy > this->_EnergyCutOff )
    {
        return energyTooHigh;
    }
    for ( SuperposableConformationCollection_O::entryIterator ei=this->begin_Entries(); ei!=this->end_Entries(); ei++ )
    {
        if ( (*ei)->getData()->getDouble("energy") > newEnergy )
	{
	    this->_Entries.insert(ei,entry);
	    return success;
	}
    }
    this->_Entries.append(entry);
    return success;
}








#ifdef XML_ARCHIVE
    void	ConformationalSearchDatabase_O::archive(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute( "EnergyCutOff", this->_EnergyCutOff );
    node->attribute( "DisplayProgress", this->_DisplayProgress );
}
#endif


/*!
 * Orient every structure from 1 on to the first one
 */
void	ConformationalSearchDatabase_O::finishSearch()
{
    geom::CoordinateArray_sp			newConf;
double					energy;
chem::SuperposeEngine_sp			superposer;
ConformationCollectionEntry_sp		firstEntry, entry;
gctools::Vec0<ConformationCollectionEntry_sp>::iterator	ei;
Matrix				transform;
uint					i;

    if ( this->numberOfEntries() <= 0 ) return;

	// Get a superposer that compares structures to the first one
    firstEntry = this->firstEntry();
    superposer = SuperposeEngine_O::create();
    this->setEntryCoordinatesAsFixedWithinSuperposeEngine(firstEntry,superposer);

	//
	// Get the first structure
	//
    firstEntry->writeCoordinatesToMatter(this->_Matter);
    firstEntry->getData()->setDouble("energyWrtFirst",0.0);
    firstEntry->getData()->setDouble("rmsWrtFirst",0.0);
    energy = firstEntry->getData()->getDouble("energy");
    LOG(BF("Entry 0  energy = %lf") % energy  );
    	//
	//
    	// Now check if the energy is within the limits
	//
    i=0;
    double rms;
    gctools::Vec0<ConformationCollectionEntry_sp>::iterator	closestIterator;
    uint highIndex = this->_Entries.size()-1;
    uint lastUsefulEntry = highIndex;
    double lowestEnergy = (*(this->_Entries.begin()))->getData()->getDouble("energy");
    LOG(BF("Looking for entry with energy closest to energy(%lf)") % energy );
    for ( i=1,ei=this->begin_Entries()+1; ei!=this->end_Entries(); i++,ei++ )
    {
        energy = (*ei)->getData()->getDouble("energy");
	LOG(BF("Looking at existing entry(%3d) with energy(%lf)") % i % energy  );
	double deltaEnergy = energy-lowestEnergy;
		// If delta > this->_EnergyCutOff then throw this and following
		// entries out
	if ( deltaEnergy > this->_EnergyCutOff )
	{
	    lastUsefulEntry = i-1;
	    break;
	}
        this->setEntryCoordinatesAsMoveableWithinSuperposeEngine(*ei,superposer);
	transform = superposer->superpose();
	rms = superposer->rootMeanSquareDifference();
	(*ei)->getData()->setDouble("energyWrtFirst",deltaEnergy);
	(*ei)->getData()->setDouble("rmsWrtFirst",rms);
        (*ei)->writeCoordinatesToMatter(this->_Matter);
	this->_Matter->applyTransformToAtoms(transform);
	(*ei)->extractCoordinatesFromMatter(this->_Matter);
	if ( this->_DisplayProgress)
	{
	    _lisp->print(BF("CONFSEARCH: Entry %d  energy(%lf) rms(%lf)") % i% energy % rms );
	}
    }

    if ( lastUsefulEntry < this->_Entries.size()-1 )
    {
	if ( this->_DisplayProgress )
	{
	    _lisp->print(BF( "CONFSEARCH: Trimming search entry table to %d elements") % (lastUsefulEntry+1) );
	}
	LOG(BF("CONFSEARCH: Trimming search entry table to %d elements") % (lastUsefulEntry+1 ) );
	this->_Entries.resize(lastUsefulEntry+1);
    }
}









};

