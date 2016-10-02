/*
    File: conformationExplorer.cc
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
#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	Maintain a database of conformational searches
 *
 */
#include <clasp/core/common.h>
#include <clasp/core/str.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/conformationExplorer.h>
#include <cando/chem/loop.h>
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include <clasp/core/hashTableEq.h>
#include <cando/geom/color.h>
#include <cando/chem/superposeEngine.h>
#include <clasp/core/sort.h>
#include <clasp/core/binder.h>
#include <clasp/core/intArray.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <cando/geom/coordinateArray.h>
#include <clasp/core/numbers.h>
#include <clasp/core/wrappers.h>



namespace chem
{


    void ConformationExplorerEntryStage_O::initialize()
    {
	this->Base::initialize();
	this->_Binder = core::HashTableEq_O::create_default();
	this->_FinalCoordinates = geom::CoordinateArray_O::create();
	this->_WeakConformationExplorerEntry = _Nil<ConformationExplorerEntry_O>();
	this->_Complete = true;
	this->_EnergyKCal = 0.0;
	this->_ExternalInterfaceName = "";
	this->_Model = "";
	this->_StageName = "undefStage";
    }


    ConformationExplorerEntry_sp	ConformationExplorerEntryStage_O::getConformationExplorerEntry()
    {
	ConformationExplorerEntry_sp	entry;
	ASSERTNOTNULL(this->_WeakConformationExplorerEntry);
	entry = this->_WeakConformationExplorerEntry;
	return entry;
    }




CL_LISPIFY_NAME("setSuperposableCoordinatesAsFixedWithinSuperposeEngine");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsFixedWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setFixedPoints(explorer->_getSuperposeAtomIndices(), this->_FinalCoordinates );
    }


CL_LISPIFY_NAME("setSuperposableCoordinatesAsMoveableWithinSuperposeEngine");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsMoveableWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setMoveablePoints(explorer->_getSuperposeAtomIndices(), this->_FinalCoordinates );
    }



CL_LISPIFY_NAME("setFinalCoordinatesAsFixedWithinSuperposeEngine");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::setFinalCoordinatesAsFixedWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setFixedAllPoints(this->_FinalCoordinates);
    }


CL_LISPIFY_NAME("setFinalCoordinatesAsMoveableWithinSuperposeEngine");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::setFinalCoordinatesAsMoveableWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setMoveableAllPoints(this->_FinalCoordinates);
    }



    ConformationExplorer_sp	ConformationExplorerEntryStage_O::getConformationExplorer()
    {
	ConformationExplorerEntry_sp	entry;
	ASSERTNOTNULL(this->_WeakConformationExplorerEntry);
	entry = this->getConformationExplorerEntry();
	return entry->getConformationExplorer();
    }


    void ConformationExplorerEntryStage_O::alignCoordinatesTo(ConformationExplorerEntryStage_sp entryStage )
    {
	SuperposeEngine_sp superposer = SuperposeEngine_O::create();
	entryStage->setSuperposableCoordinatesAsFixedWithinSuperposeEngine(superposer);
	this->setSuperposableCoordinatesAsMoveableWithinSuperposeEngine(superposer);
	Matrix transform = superposer->superpose();
	this->_FinalCoordinates->transform(transform);
    }





#if 0
    geom::Render_sp	ConformationExplorerEntryStage_O::rendered(core::List_sp opts)
    {
	ConformationExplorer_sp		explorer;
	geom::DisplayList_sp			dlAll;
	GrPickableMatter_sp		dlMatter;
	Matter_sp			matter;
	geom::CoordinateArray_sp		superposeCoords;
	core::IntArray_O::iterator		ii;
	explorer = this->getConformationExplorer();
	dlAll = geom::DisplayList_O::create();
	FIX_ME(); dlAll->setName(_lisp->internKeyword(this->_StageName));
	{
	    LOG(BF("Writing coordinates to matter"));
	    matter = explorer->getMatter();
	    this->writeCoordinatesToMatter(matter);
	    dlMatter = GrPickableMatter_O::create();
	    dlMatter->setFromMatter(matter);
	    dlMatter->setName(_lisp->internKeyword("molFrame"));
	    dlAll->append(dlMatter);
#ifdef DEBUG_ON
	    core::SymbolVector_sp sv = dlMatter->getFullName();
	    LOG(BF("At time of creation dlMatter fullName=%s") % _rep_(sv) );
#endif
	}
	{
	    LOG(BF("Rendering superposed atoms as spheres" ));
	    // Render the superposed atom positions
	    geom::DisplayList_sp prims = geom::DisplayList_O::create();
	    FIX_ME(); prims->setName(_lisp->internKeyword("superposeAtoms"));
	    geom::GrSpheres_sp sphereList = geom::GrSpheres_O::create();
	    for ( ii=explorer->superposeAtomIndices_begin(); ii!=explorer->superposeAtomIndices_end(); ii++ )
	    {
		sphereList->appendSphere(explorer->_getAtomAtIndex(*ii)->getPosition(),0.2,geom::Color_O::yellow(_lisp));
	    }
	    prims->append(sphereList);
	    dlAll->append(prims);
	}
	geom::RenderInformation_sp info;
	stringstream sinfo;
	sinfo << "=====================================" << std::endl;
	sinfo << "Stage name: "<< this->_StageName << std::endl;
	sinfo << "Based on stage: "<< this->_BasedOnStageName << std::endl;
	sinfo << "Is complete: " << this->_Complete << std::endl;
	sinfo << "Energy (kCal/mole): " << fixed << setw(10) << setprecision(1) 
	      << this->_EnergyKCal << std::endl;
	sinfo << "Generated by package: " << this->_ExternalInterfaceName << std::endl;
	sinfo << "Model: " << this->_Model << std::endl;
	SYMBOL_EXPORT_SC_(ChemPkg,convergence);
	if ( this->getBinder()->contains(_sym_convergence))
	{
	    sinfo << this->getBinder()->lookup(_sym_convergence).as<core::Str_O>()->get();
	}
	info = geom::RenderInformation_O::create(sinfo.str());
	dlAll->append(info);
	geom::Render_sp dlDict = this->getBinder()->rendered(opts);
	dlAll->append(dlDict);
	return dlAll;
    }
#endif





#ifdef XML_ARCHIVE
    void	ConformationExplorerEntryStage_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("StageName",this->_StageName);
	node->attributeIfNotDefault<string>("BasedOnStageName",this->_BasedOnStageName,"");
	node->attributeIfNotDefault<string>("ExternalInterfaceName",this->_ExternalInterfaceName,"");
	node->attributeIfNotDefault<string>("Model",this->_Model,"");
	node->attributeIfNotDefault("EnergyKCal",this->_EnergyKCal,0.0);
	node->attributeIfNotDefault("Complete",this->_Complete,true);
	node->archiveWeakPointer( "ConformationExplorerEntry",
				  this->_WeakConformationExplorerEntry );
	node->attribute("AllCoordinates", this->_FinalCoordinates );
	node->attribute("Binder",this->_Binder);
    }
#endif




    void	ConformationExplorerEntryStage_O::setConformationExplorerEntry(ConformationExplorerEntry_sp s)
    {
	geom::CoordinateArray_sp	ca;
	ConformationExplorer_sp	explorer;
	this->_WeakConformationExplorerEntry = s;
	explorer = this->getConformationExplorer();
	ca = geom::CoordinateArray_O::create(explorer->numberOfAllAtoms());
	this->_FinalCoordinates = ca;
    }


    void	ConformationExplorerEntryStage_O::setFinalCoordinates(geom::CoordinateArray_sp ac)
    {
	ASSERTNOTNULL(ac);
	LOG(BF("setFinalCoordinates:%s") % (ac->asXmlString().c_str() ) );
	LOG(BF("The address of the geom::CoordinateArray_sp is in o") );
	this->_FinalCoordinates = ac;
    }



    void	ConformationExplorerEntryStage_O::setCoordinateForAtom(Atom_sp atom, Vector3 pos )
    {
	int	idx;
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	idx = explorer->indexOfAtom(atom);
	ASSERTNOTNULL(this->_FinalCoordinates);
	ASSERT(this->_FinalCoordinates->size() == explorer->numberOfAllAtoms());
	this->_FinalCoordinates->setElement(idx,pos);
    }


CL_LISPIFY_NAME("writeCoordinatesToMatter");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::writeCoordinatesToMatter(Matter_sp agg)
    {
	ConformationExplorer_sp			sl;
	vector<Vector3>::iterator	ci;
        gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
	sl = this->getConformationExplorer();
	for ( ai=sl->begin_AllAtoms(),ci=this->_FinalCoordinates->begin(); ai!=sl->end_AllAtoms(); ai++, ci++ )
	{
	    (*ai)->setPosition(*ci);
	}
    }


CL_LISPIFY_NAME("extractCoordinatesFromMatter");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::extractCoordinatesFromMatter(Matter_sp matter)
    {
	ConformationExplorer_sp	sl;
	vector<Vector3>::iterator	ci;
	LOG(BF("About to get ConformationExplorer") );
	sl = this->getConformationExplorer();
	LOG(BF("Got ConformationExplorer") );
	ASSERTNOTNULL(sl);
	geom::CoordinateArray_sp coords = sl->_extractCoordinateArray(matter);
	ASSERTNOTNULL(coords);
	this->_FinalCoordinates = coords;
    }




















    void	ConformationExplorerEntry_O::initialize()
    {
	this->Base::initialize();
	this->_Binder = core::HashTableEq_O::create_default();
	this->_UniqueEntryIndex = -1;
    }


#ifdef XML_ARCHIVE
    void	ConformationExplorerEntry_O::archiveBase(core::ArchiveP node)
    {
	node->archiveWeakPointer( "ConformationExplorer",
				  this->_WeakConformationExplorer );
	node->archiveVector0("Stages", this->_Stages);
	node->attributeIfNotNil("SelectedStage",this->_SelectedStage);
	node->attribute("Binder",this->_Binder);
	node->attributeIfNotDefault("uidx",this->_UniqueEntryIndex,-1);
    }
#endif

CL_LISPIFY_NAME("getSelectedStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getSelectedStage()
    {
	if ( !this->hasSelectedStage() )
	{
	    return _Nil<ConformationExplorerEntryStage_O>();
	}
	return this->_SelectedStage;
    }

CL_LISPIFY_NAME("setSelectedStage");
CL_DEFMETHOD     void	ConformationExplorerEntry_O::setSelectedStage(ConformationExplorerEntryStage_sp stage)
    {
	ConformationExplorerEntry_sp	stageEntry;
	stageEntry = stage->getConformationExplorerEntry();
	ASSERTP(stageEntry.get()==this,"StageEntry != this");
	this->_SelectedStage = stage;
    }



CL_LISPIFY_NAME("createEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp	ConformationExplorerEntry_O::createEntryStage(string const& name)
    {
	ASSERTF(!this->hasEntryStageWithName(name),
		BF("Stage with key[%s] already exists!") % name );
	GC_ALLOCATE(ConformationExplorerEntryStage_O, stage );
	stage->setConformationExplorerEntry(this->sharedThis<ConformationExplorerEntry_O>());
	stage->setStageName(name);
	this->_Stages.push_back(stage);
	return stage;
    }


    void	ConformationExplorerEntry_O::setConformationExplorer(ConformationExplorer_sp s)
    {
	geom::CoordinateArray_sp	ca;
	this->_WeakConformationExplorer = s;
    }




CL_LISPIFY_NAME("lastEntryStageName");
CL_DEFMETHOD     string	ConformationExplorerEntry_O::lastEntryStageName()
    {
	if ( this->_Stages.size() == 0 ) return "";
	return this->_Stages[this->_Stages.size()-1]->getStageName();
    }



CL_LISPIFY_NAME("getLastEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getLastEntryStage()
    {
	ASSERTP(this->_Stages.size()>0,"There must be at least one complete stage");
	int idx = this->_Stages.size()-1;
	return this->_Stages[idx];
    }


CL_LISPIFY_NAME("alreadyHasLastCompleteStage");
CL_DEFMETHOD     bool	ConformationExplorerEntry_O::alreadyHasLastCompleteStage(const string& stageName )
    {
	ConformationExplorerEntryStage_sp	lastStage;
	lastStage = this->getLastCompleteEntryStage();
	return( lastStage->getStageName() == stageName );
    }


CL_LISPIFY_NAME("getOrCreateLastIncompleteEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getOrCreateLastIncompleteEntryStage(const string& stageName)
    {
	ConformationExplorerEntryStage_sp	lastStage, newStage;
	lastStage = this->getLastEntryStage();
	if ( lastStage->isComplete() )
	{
	    if ( lastStage->getStageName() != stageName )
	    {
		newStage = this->createEntryStage(stageName);
		newStage->setComplete(false);
	    } else
	    {
		newStage = _Nil<ConformationExplorerEntryStage_O>();
		return newStage;
	    }
	} else
	{
	    newStage = lastStage;
	    newStage->setStageName(stageName);
	    newStage->setComplete(false);
	}
	return newStage;
    }


CL_LISPIFY_NAME("getLastCompleteEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getLastCompleteEntryStage()
    {
	ASSERTP(this->_Stages.size()>0,"There must be at least one complete stage");
	int idx = this->_Stages.size()-1;
	if ( this->_Stages[idx]->isComplete() )
	{
	    return this->_Stages[idx];
	}
	ASSERTP(idx!=0, "The first stage must always be complete");
	idx--;
	if ( this->_Stages[idx]->isComplete() )
	{
	    return this->_Stages[idx];
	}
	SIMPLE_ERROR(BF("Only the last stage can be incomplete"));
    }







CL_LISPIFY_NAME("hasEntryStageWithName");
CL_DEFMETHOD     bool	ConformationExplorerEntry_O::hasEntryStageWithName(const string& key)
    {
	stageIterator	si;
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{
	    if ( (*si)->getStageName() == key ) return true;
	}
	return false;
    }


CL_LISPIFY_NAME("getEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp	ConformationExplorerEntry_O::getEntryStage(const string& key)
    {
	stageIterator	si;
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{
	    if ( (*si)->getStageName() == key )
	    {
		return (*si);
	    }
	}
	SIMPLE_ERROR(BF("Could not find key: "+key));
    }



#if 0
    geom::Render_sp	ConformationExplorerEntry_O::rendered(core::List_sp opts)
    {
	geom::DisplayList_sp	dlAll;
	geom::Render_sp		dlStage;
	stageIterator		si;
	dlAll = geom::DisplayList_O::create();
	FIX_ME(); dlAll->setName(_lisp->internKeyword("stages"));
	LOG(BF("Rendering %d stages")%this->_Stages.size());
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{_BLOCK_TRACEF(BF("Rendering stage: %s")%(*si)->getStageName());
	    if ( (*si)->isComplete() )
	    {
		LOG(BF("It was complete"));
		dlStage = (*si)->rendered(opts);
		dlAll->append(dlStage);
	    } else
	    {
		LOG(BF("It was NOT complete"));
	    }
	}
	return dlAll;
    }
#endif


    adapt::StringSet_sp ConformationExplorerEntry_O::getEntryStageNames()
    {
	adapt::StringSet_sp stageNames = adapt::StringSet_O::create();
	stageIterator it;
	for ( it=this->begin_Stages(); it!=this->end_Stages(); it++ )
	{
	    stageNames->insert((*it)->getStageName());
	}
	return stageNames;
    }

    void	ConformationExplorer_O::initialize()
    {
	this->Base::initialize();
	this->_Matter = _Nil<Matter_O>();
	this->_AllAtoms.clear();
	this->_Binder = core::HashTableEq_O::create_default();
	this->clearEntries();
	this->_SuperposeAtomIndices = core::IntArray_O::create();
    }


CL_LISPIFY_NAME("clearEntries");
CL_DEFMETHOD     void ConformationExplorer_O::clearEntries()
    {
	this->_Entries.clear();
    }




    Atom_sp	ConformationExplorer_O::_getAtomAtIndex(unsigned i)
    {
	ASSERT_lessThan(i,this->_AllAtoms.size());
	return this->_AllAtoms[i];
    }



    geom::CoordinateArray_sp ConformationExplorer_O::_extractCoordinateArray(Matter_sp agg)
    {
	vector<Vector3>::iterator	ci;
        gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
#ifdef	DEBUG_ConformationExplorerEntry
	this->_Status->addMessage("extractCoordinatesFromMatter");
#endif
	LOG(BF("About to get ConformationExplorer"));
	LOG(BF("About to iterate over atoms and get positions"));
	uint numberOfAtoms = agg->numberOfAtoms();
	geom::CoordinateArray_sp coords = geom::CoordinateArray_O::create(numberOfAtoms);
	ASSERTNOTNULL(coords);
	uint aidx = 0;
	for ( ai=this->begin_AllAtoms(),ci=coords->begin(); ai!=this->end_AllAtoms(); ai++, ci++ )
	{
	    LOG(BF("Getting position of ATOM(%d) %s") % aidx % (*ai)->description().c_str() );
	    if ( !(*ai)->isContainedBy(agg) )
	    {
		SIMPLE_ERROR(BF("You asked a ConformationExplorer to extract coordinates from the Matter(%s) but the ConformationExplorer atom(%s) is not contained within that Matter")%agg->description()%(*ai)->description());
	    } else
	    {
#ifdef DEBUG_ON
		LOG(BF("You asked a ConformationExplorer to extract coordinates from the Matter(%s) and the ConformationExplorer atom(%s) IS contained within that Matter")%agg->description()%(*ai)->description());
#endif
	    }
	    *ci = (*ai)->getPosition();
	    LOG(BF("Position = %s") % (*ci).asString() );
	    aidx++;
	}
	LOG(BF("Done"));
	return coords;
    }


CL_LISPIFY_NAME("setMatter");
CL_DEFMETHOD     void	ConformationExplorer_O::setMatter(Matter_sp matter)
    {
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




CL_LISPIFY_NAME("appendEntry");
CL_DEFMETHOD     void	ConformationExplorer_O::appendEntry(ConformationExplorerEntry_sp entry)
    {
	entry->setConformationExplorer(this->sharedThis<ConformationExplorer_O>());
	this->_Entries.push_back(entry);
    }


CL_LISPIFY_NAME("entriesAsList");
CL_DEFMETHOD     core::List_sp ConformationExplorer_O::entriesAsList()
    {
	core::Cons_sp first = core::Cons_O::create(_Nil<core::T_O>(), _Nil<core::T_O>());
	core::Cons_sp cur = first;
	ConformationExplorer_O::entryIterator si;
	for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
	{
	    core::Cons_sp newCons = core::Cons_O::create(*si,_Nil<core::T_O>());
	    cur->setCdr(newCons);
	    cur = newCons;
	}
	return first->cdr();
    }









#ifdef XML_ARCHIVE
    void	ConformationExplorer_O::archiveBase(core::ArchiveP node)
    {
	if ( !node->loading() )
	{
	    ASSERTNOTNULL(this->_Matter);
	}
	node->attribute( "Matter", this->_Matter );
	node->archiveVector0( "Entries", this->_Entries );
	node->archiveOrderedSet( "AllAtoms", this->_AllAtoms );
	node->attribute( "SuperposeAtomIndices", this->_SuperposeAtomIndices );
	node->attributeIfNotNil( "Binder", this->_Binder );
    }
#endif



#if 0
    geom::Render_sp	ConformationExplorer_O::rendered(core::List_sp opts)
    {
	geom::FrameList_sp	frames;
	geom::Render_sp		dlEntry;
	_LINE();
	ConformationExplorer_O::entryIterator	si;
	ConformationExplorerEntry_sp	entry;
	geom::CoordinateArray_sp		superposeCoords;
	geom::CoordinateArray_O::iterator	ci;
	_LINE();
	geom::DisplayList_sp dlAll = geom::DisplayList_O::create();
	_LINE();
	frames = geom::FrameList_O::create();
	FIX_ME(); frames->setName(_lisp->internKeyword("conformation"));
	_LINE();
	if ( this->numberOfEntries() == 0 )
	{
	    LOG(BF("FLAG! There are no entries"));
	    ASSERTNOTNULL(this->_Matter);
	    if ( this->_Matter.notnilp() )
	    {
		LOG(BF("Rendering ConformationExplorer::matter"));
		dlEntry = this->_Matter->rendered(opts);
		ASSERTNOTNULL(dlEntry);
		frames->append(dlEntry);
	    }
	} else
	{
	    LOG(BF("FLAG! There are %d entries") % this->_Entries.size());
	    for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
	    {_BLOCK_TRACEF(BF("Rendering entry"));
		dlEntry = (*si)->rendered(opts);
		LOG(BF("The name of the Frame entry is[%s]") % _rep_(dlEntry->getName()) );
		frames->append(dlEntry);
	    }
	}
	_LINE();
	dlAll->append(frames);
	_LINE();
	dlAll->append(this->_Binder->rendered(opts));
	_LINE();
	return dlAll;
    }
#endif



CL_LISPIFY_NAME("clearSuperposeAtoms");
CL_DEFMETHOD     void	ConformationExplorer_O::clearSuperposeAtoms()
    {
	ASSERTNOTNULL(this->_SuperposeAtomIndices);
	this->_SuperposeAtomIndices->clear();
    }



CL_LISPIFY_NAME("addSuperposeAtom");
CL_DEFMETHOD     void	ConformationExplorer_O::addSuperposeAtom(Atom_sp a)
    {
	atomIterator		ai;
	int			idx;
	ASSERTNOTNULL(this->_SuperposeAtomIndices);
	for ( idx=0,ai=this->begin_AllAtoms(); ai!=this->end_AllAtoms(); idx++,ai++ )
	{
	    if ( (*ai) == a )
	    {
		this->_SuperposeAtomIndices->append(idx);
		return;
	    }
	}
	SIMPLE_ERROR(BF("You cannot add superpose atom: "+a->description()+" because it isn't in the ConformationCollection"));
    }


    void	ConformationExplorer_O::superposeAtomsFromSet(gctools::SmallOrderedSet<Atom_sp> atoms)
    {
        gctools::SmallOrderedSet<Atom_sp>::iterator	ai;
	this->clearSuperposeAtoms();
	for ( ai=atoms.begin(); ai!=atoms.end(); ai++ )
	{
	    this->addSuperposeAtom(*ai);
	}
    }

CL_LISPIFY_NAME("superposeAllHeavyAtoms");
CL_DEFMETHOD     void	ConformationExplorer_O::superposeAllHeavyAtoms()
    {
	atomIterator		ai;
	int			idx;
	ASSERTNOTNULL(this->_SuperposeAtomIndices);
	this->clearSuperposeAtoms();
	for ( idx=0,ai=this->begin_AllAtoms(); ai!=this->end_AllAtoms(); idx++,ai++ )
	{
	    if ( (*ai)->getElement() != element_H )
	    {
		this->_SuperposeAtomIndices->append(idx);
	    }
	}
    }


    gctools::SmallOrderedSet<Atom_sp>	ConformationExplorer_O::getSuperposeAtoms()
    {
        gctools::SmallOrderedSet<Atom_sp>	result;
	core::IntArray_O::iterator	si;
	Atom_sp			a;
	ASSERTNOTNULL(this->_SuperposeAtomIndices);
	for ( si=this->_SuperposeAtomIndices->begin();si!=this->_SuperposeAtomIndices->end(); si++ )
	{
	    a = this->_AllAtoms[*si];
	    result.insert(a);
	}
	return result;
    }


    uint	ConformationExplorer_O::numberOfSuperposeAtoms()
    {
	return this->_SuperposeAtomIndices->size();
    }


    struct OrderStageAndKeyValue
    {
	string _StageName;
	core::Symbol_sp _KeyName;

	bool operator()(ConformationExplorerEntry_sp x,ConformationExplorerEntry_sp y )
	{
	    ConformationExplorerEntryStage_sp xStage, yStage;
	    xStage = x->getEntryStage(this->_StageName);
	    yStage = y->getEntryStage(this->_StageName);
	    core::T_sp xValue = xStage->getBinder()->gethash(this->_KeyName);
	    core::T_sp yValue = yStage->getBinder()->gethash(this->_KeyName);
	    if ( xValue.as<core::General_O>()->operator<(yValue) ) return true;
	    return false;
	}
    };




CL_LISPIFY_NAME("ConformationExplorer-sort");
CL_DEFMETHOD     void	ConformationExplorer_O::sort(const string& stageName, core::Symbol_sp keyName )
    {
	OrderStageAndKeyValue comparer;
	comparer._StageName = stageName;
	comparer._KeyName = keyName;
	sort::quickSort(this->_Entries.begin(), this->_Entries.end(), comparer);
    }







CL_LISPIFY_NAME("saveAs");
CL_DEFMETHOD     void	ConformationExplorer_O::saveAs(const string& fn)
    {
	IMPLEMENT_ME();
#ifdef XML_ARCHIVE
	core::XmlSaveArchive_sp	xml;
	xml = core::XmlSaveArchive_O::create();
	xml->put("ConformationExplorer",this->sharedThis<ConformationExplorer_O>());
	xml->saveAs(fn);
#endif
    }

CL_LISPIFY_NAME("getEntryIndex");
CL_DEFMETHOD     unsigned ConformationExplorer_O::getEntryIndex(ConformationExplorerEntry_sp entry)
    {
	int		i;
	entryIterator	ei;
	for ( i=0,ei=this->_Entries.begin(); ei!=this->_Entries.end(); ei++,i++ )
	{
	    if ( (*ei) == entry )
	    {
		return i;
	    }
	}
	SIMPLE_ERROR(BF("Could not find entry in ConformationExplorer"));
    }




CL_LISPIFY_NAME("getEntry");
CL_DEFMETHOD     ConformationExplorerEntry_sp	ConformationExplorer_O::getEntry(unsigned i)
    {
	ASSERT_lessThan(i,this->_Entries.size());
	ConformationExplorerEntry_sp e;
	e = this->_Entries[i];
	return e;
    }

CL_LISPIFY_NAME("firstEntry");
CL_DEFMETHOD     ConformationExplorerEntry_sp	ConformationExplorer_O::firstEntry()
    {
	ASSERT(this->_Entries.size() > 0 );
	return this->_Entries[0];
    }





CL_LISPIFY_NAME("createEntry");
CL_DEFMETHOD     ConformationExplorerEntry_sp	ConformationExplorer_O::createEntry()
    {
	GC_ALLOCATE(ConformationExplorerEntry_O, entry );
	entry->setConformationExplorer(this->sharedThis<ConformationExplorer_O>());
	entry->_UniqueEntryIndex = this->numberOfEntries();
	return entry;
    }





    bool	ConformationExplorer_O::hasStageNameInAllEntries(const string& stageKey)
    {
	entryIterator		ei;
	bool			inAll;
	inAll = false;
	for ( ei=this->begin_Entries(); ei!=this->end_Entries(); ei++ )
	{
	    if ( !(*ei)->hasEntryStageWithName(stageKey) )
	    {
		return false;
	    }
	}
	return true;
    }




    bool ConformationExplorer_O::findMostSimilarConformationEntryStageWithStageName(
	Matter_sp 	matter,
	const string&	stageName,
	double&		bestRms,
	ConformationExplorerEntryStage_sp& bestStage,
	uint&		bestEntryIndex )
    {
	ConformationExplorerEntryStage_sp	stage;
	double					rms;
	bool					gotBest;
	SuperposeEngine_sp			superposer;
	geom::CoordinateArray_sp			matterConf;
	bestStage = _Nil<ConformationExplorerEntryStage_O>();
	bestRms = 9.9e99;
	gotBest = false;
	uint numSuperposeAtoms = this->numberOfSuperposeAtoms();
	if ( numSuperposeAtoms == 0 )
	{
	    SIMPLE_ERROR(BF("There must be at least one superpose atom"));
	} else if ( numSuperposeAtoms == 1 || numSuperposeAtoms == 2 )
	{
	    uint entryIndex = 0;
	    entryIterator ei;
	    for ( ei=this->begin_Entries(); ei!=this->end_Entries(); ei++, entryIndex++ )
	    {
		if ( (*ei)->hasEntryStageWithName(stageName) )
		{
		    stage = (*ei)->getEntryStage(stageName);
		    bestRms = 0.0;
		    bestStage = stage;
		    bestEntryIndex = entryIndex;
		    gotBest = true;
		}
	    }
	} else if ( numSuperposeAtoms >= 3 )
	{
	    superposer = SuperposeEngine_O::create();
	    matterConf = this->_extractCoordinateArray(matter);
	    superposer->setFixedPoints(this->_SuperposeAtomIndices,matterConf);
	    uint entryIndex = 0;
	    entryIterator ei;
	    for ( ei=this->begin_Entries(); ei!=this->end_Entries(); ei++, entryIndex++ )
	    {
		if ( (*ei)->hasEntryStageWithName(stageName) )
		{
		    stage = (*ei)->getEntryStage(stageName);
		    stage->setSuperposableCoordinatesAsMoveableWithinSuperposeEngine(superposer);
		    superposer->superpose();
		    rms = superposer->rootMeanSquareDifference();
		    if ( rms < bestRms )
		    {
			bestRms = rms;
			bestStage = stage;
			bestEntryIndex = entryIndex;
			gotBest = true;
		    }
#if 0
		} else
		{
		    return false;
#endif
		}
	    }
	}
	return gotBest;
    }

CL_LISPIFY_NAME("getClosestMatchingEntry");
CL_DEFMETHOD     ConformationExplorerEntry_sp ConformationExplorer_O::getClosestMatchingEntry(ConformationExplorerMatch_sp match)
    {
	uint idx = match->getClosestMatchIndex();
	return this->getEntry(idx);
    }

CL_LISPIFY_NAME("getClosestMatchingEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp ConformationExplorer_O::getClosestMatchingEntryStage(ConformationExplorerMatch_sp match)
    {
	uint idx = match->getClosestMatchIndex();
	ConformationExplorerEntry_sp e = this->getEntry(idx);
	ConformationExplorerEntryStage_sp s = e->getEntryStage(match->getMatchStageName());
	return s;
    }



    string ConformationExplorer_O::__str__()
    {
	stringstream ss;
	ss << this->Base::__str__();
	ss << "#ConformationExplorerEntries: " << this->_Entries.size() << std::endl;
	adapt::StringSet_sp stageNames = adapt::StringSet_O::create();
        gctools::Vec0<ConformationExplorerEntry_sp>::const_iterator ei;
	for ( ei=this->_Entries.begin(); ei!=this->_Entries.end(); ei++ )
	{
	    adapt::StringSet_sp oneEntryStageNames = (*ei)->getEntryStageNames();
	    stageNames->insertStringSet(oneEntryStageNames);
	}
	ss << "EntryStageNames: %s\n", stageNames->asString().c_str();
	return ss.str();
    }

CL_LISPIFY_NAME("findClosestMatchingConformation");
CL_DEFMETHOD     void ConformationExplorer_O::findClosestMatchingConformation(
	ConformationExplorerMatch_sp match,
	Matter_sp matter, const string& stageKey,
	double rms )
    {
	double					bestRms;
	ConformationExplorerEntryStage_sp	bestStage;
	uint					bestIndex;
	match->setMatches(false);
	match->setMatchStageName(stageKey);
	bool gotBest = this->findMostSimilarConformationEntryStageWithStageName(matter,stageKey,bestRms,bestStage,bestIndex);
	if ( !gotBest )
	{
	    LOG(BF("There were no entries to match") );
	    return;
	}
	if ( bestRms < rms )
	{
	    LOG(BF("The best match had an rms(%lf) and we are looking for one better than(%lf)") % bestRms % rms  );
	    match->setMatches(true);
	}
	match->setClosestMatchRms(bestRms);
	match->setClosestMatchIndex(bestIndex);
	return;
    }



    uint	ConformationExplorer_O::indexOfAtom(Atom_sp atom)
    {
	return this->_AllAtoms.indexOf(atom);
    }



CL_LISPIFY_NAME("alignAllConformationsToTheFirstForStage");
CL_DEFMETHOD     void	ConformationExplorer_O::alignAllConformationsToTheFirstForStage(const string& stageName )
    {
	if ( !this->hasStageNameInAllEntries(stageName) )
	{
	    SIMPLE_ERROR(BF("You have to have the stageName(%s) in all entries")%stageName);
	}
	if ( this->_Entries.size() == 0 )
	{
	    SIMPLE_ERROR(BF("You have to have at least one entry"));
	}
	if ( this->_Entries.size() < 2 ) return;
	ConformationExplorer_O::entryIterator first = this->_Entries.begin();
	ConformationExplorerEntryStage_sp firstEntryStage = (*first)->getEntryStage(stageName);
	ConformationExplorer_O::entryIterator rest = first;
	rest++;
	for ( ; rest!=this->end_Entries(); rest++ )
	{
	    ConformationExplorerEntryStage_sp entryStage = (*rest)->getEntryStage(stageName);
	    (entryStage)->alignCoordinatesTo(firstEntryStage);
	}
    }


    class	OrderByEnergyAscending
    {
    public:
	string	_StageName;
    public:
	bool operator()(ConformationExplorerEntry_sp x, ConformationExplorerEntry_sp y )
	{
	    ConformationExplorerEntryStage_sp xStage = x->getEntryStage(this->_StageName);
	    ConformationExplorerEntryStage_sp yStage = y->getEntryStage(this->_StageName);
	    return xStage->getEnergyKCal()<yStage->getEnergyKCal();
	}
    };


CL_LISPIFY_NAME("sortByEnergyAscendingForStage");
CL_DEFMETHOD     void	ConformationExplorer_O::sortByEnergyAscendingForStage(const string& stageName )
    {
	OrderByEnergyAscending order;
	order._StageName = stageName;
	sort::quickSort(this->_Entries.begin(),this->_Entries.end(),order);
    }

    void	ConformationExplorerMatch_O::initialize()
    {
	this->_Matches = false;
    }


















};
