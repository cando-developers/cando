/*
    File: confSearchDatabase.cc
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
    this->_EnergyCutOff = 5.0; // kcal/mole
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
    geom::SimpleVectorCoordinate_sp			newConf;
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
    LOG("Entry 0  energy = %lf" , energy  );
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
    LOG("Looking for entry with energy closest to energy(%lf)" , energy );
    for ( i=1,ei=this->begin_Entries()+1; ei!=this->end_Entries(); i++,ei++ )
    {
        energy = (*ei)->getData()->getDouble("energy");
	LOG("Looking at existing entry(%3d) with energy(%lf)" , i , energy  );
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
          core::clasp_write_string(fmt::format("CONFSEARCH: Entry {}  energy({}) rms({})\n" , i , energy , rms ));
	}
    }

    if ( lastUsefulEntry < this->_Entries.size()-1 )
    {
	if ( this->_DisplayProgress )
	{
          core::clasp_write_string(fmt::format( "CONFSEARCH: Trimming search entry table to {} elements\n" , (lastUsefulEntry+1) ));
	}
	LOG("CONFSEARCH: Trimming search entry table to %d elements" , (lastUsefulEntry+1 ) );
	this->_Entries.resize(lastUsefulEntry+1);
    }
}









};
