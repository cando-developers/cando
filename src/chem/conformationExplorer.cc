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
    {_G();
	ConformationExplorerEntry_sp	entry;
	ASSERTNOTNULL(this->_WeakConformationExplorerEntry);
	entry = this->_WeakConformationExplorerEntry;
	return entry;
    }




    void	ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsFixedWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {_G();
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setFixedPoints(explorer->_getSuperposeAtomIndices(), this->_FinalCoordinates );
    }


    void	ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsMoveableWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {_G();
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setMoveablePoints(explorer->_getSuperposeAtomIndices(), this->_FinalCoordinates );
    }



    void	ConformationExplorerEntryStage_O::setFinalCoordinatesAsFixedWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {_G();
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setFixedAllPoints(this->_FinalCoordinates);
    }


    void	ConformationExplorerEntryStage_O::setFinalCoordinatesAsMoveableWithinSuperposeEngine( SuperposeEngine_sp superposer )
    {_G();
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	superposer->setMoveableAllPoints(this->_FinalCoordinates);
    }



    ConformationExplorer_sp	ConformationExplorerEntryStage_O::getConformationExplorer()
    {_G();
	ConformationExplorerEntry_sp	entry;
	ASSERTNOTNULL(this->_WeakConformationExplorerEntry);
	entry = this->getConformationExplorerEntry();
	return entry->getConformationExplorer();
    }


    void ConformationExplorerEntryStage_O::alignCoordinatesTo(ConformationExplorerEntryStage_sp entryStage )
    {_G();
	SuperposeEngine_sp superposer = SuperposeEngine_O::create();
	entryStage->setSuperposableCoordinatesAsFixedWithinSuperposeEngine(superposer);
	this->setSuperposableCoordinatesAsMoveableWithinSuperposeEngine(superposer);
	Matrix transform = superposer->superpose();
	this->_FinalCoordinates->transform(transform);
    }





#if 0
    geom::Render_sp	ConformationExplorerEntryStage_O::rendered(core::List_sp opts)
    {_G();
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
	    LOG(BF("At time of creation dlMatter fullName=%s") % sv->__repr__() );
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
	info = geom::RenderInformation_O::create(sinfo.str(),_lisp);
	dlAll->append(info);
	geom::Render_sp dlDict = this->getBinder()->rendered(opts);
	dlAll->append(dlDict);
	return dlAll;
    }
#endif





#ifdef XML_ARCHIVE
    void	ConformationExplorerEntryStage_O::archiveBase(core::ArchiveP node)
    {_G();
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
    {_G();
	geom::CoordinateArray_sp	ca;
	ConformationExplorer_sp	explorer;
	this->_WeakConformationExplorerEntry = s;
	explorer = this->getConformationExplorer();
	ca = geom::CoordinateArray_O::create(explorer->numberOfAllAtoms());
	this->_FinalCoordinates = ca;
    }


    void	ConformationExplorerEntryStage_O::setFinalCoordinates(geom::CoordinateArray_sp ac)
    {_G();
	ASSERTNOTNULL(ac);
	LOG(BF("setFinalCoordinates:%s") % (ac->asXmlString().c_str() ) );
	LOG(BF("The address of the geom::CoordinateArray_sp is in o") );
	this->_FinalCoordinates = ac;
    }



    void	ConformationExplorerEntryStage_O::setCoordinateForAtom(Atom_sp atom, Vector3 pos )
    {_G();
	int	idx;
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	idx = explorer->indexOfAtom(atom);
	ASSERTNOTNULL(this->_FinalCoordinates);
	ASSERT(this->_FinalCoordinates->size() == explorer->numberOfAllAtoms());
	this->_FinalCoordinates->setElement(idx,pos);
    }


    void	ConformationExplorerEntryStage_O::writeCoordinatesToMatter(Matter_sp agg)
    {_G();
	ConformationExplorer_sp			sl;
	vector<Vector3>::iterator	ci;
        gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
	sl = this->getConformationExplorer();
	for ( ai=sl->begin_AllAtoms(),ci=this->_FinalCoordinates->begin(); ai!=sl->end_AllAtoms(); ai++, ci++ )
	{
	    (*ai)->setPosition(*ci);
	}
    }


    void	ConformationExplorerEntryStage_O::extractCoordinatesFromMatter(Matter_sp matter)
    {_G();
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
    {_G();
	node->archiveWeakPointer( "ConformationExplorer",
				  this->_WeakConformationExplorer );
	node->archiveVector0("Stages", this->_Stages);
	node->attributeIfNotNil("SelectedStage",this->_SelectedStage);
	node->attribute("Binder",this->_Binder);
	node->attributeIfNotDefault("uidx",this->_UniqueEntryIndex,-1);
    }
#endif

    ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getSelectedStage()
    {_G();
	if ( !this->hasSelectedStage() )
	{
	    return _Nil<ConformationExplorerEntryStage_O>();
	}
	return this->_SelectedStage;
    }

    void	ConformationExplorerEntry_O::setSelectedStage(ConformationExplorerEntryStage_sp stage)
    {_G();
	ConformationExplorerEntry_sp	stageEntry;
	stageEntry = stage->getConformationExplorerEntry();
	ASSERTP(stageEntry.get()==this,"StageEntry != this");
	this->_SelectedStage = stage;
    }



    ConformationExplorerEntryStage_sp	ConformationExplorerEntry_O::createEntryStage(string const& name)
    {_G();
	ASSERTF(!this->hasEntryStageWithName(name),
		BF("Stage with key[%s] already exists!") % name );
	GC_ALLOCATE(ConformationExplorerEntryStage_O, stage );
	stage->setConformationExplorerEntry(this->sharedThis<ConformationExplorerEntry_O>());
	stage->setStageName(name);
	this->_Stages.push_back(stage);
	return stage;
    }


    void	ConformationExplorerEntry_O::setConformationExplorer(ConformationExplorer_sp s)
    {_G();
	geom::CoordinateArray_sp	ca;
	this->_WeakConformationExplorer = s;
    }




    string	ConformationExplorerEntry_O::lastEntryStageName()
    {_G();
	if ( this->_Stages.size() == 0 ) return "";
	return this->_Stages[this->_Stages.size()-1]->getStageName();
    }



    ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getLastEntryStage()
    {_G();
	ASSERTP(this->_Stages.size()>0,"There must be at least one complete stage");
	int idx = this->_Stages.size()-1;
	return this->_Stages[idx];
    }


    bool	ConformationExplorerEntry_O::alreadyHasLastCompleteStage(const string& stageName )
    {_G();
	ConformationExplorerEntryStage_sp	lastStage;
	lastStage = this->getLastCompleteEntryStage();
	return( lastStage->getStageName() == stageName );
    }


    ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getOrCreateLastIncompleteEntryStage(const string& stageName)
    {_G();
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


    ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getLastCompleteEntryStage()
    {_G();
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







    bool	ConformationExplorerEntry_O::hasEntryStageWithName(const string& key)
    {_G();
	stageIterator	si;
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{
	    if ( (*si)->getStageName() == key ) return true;
	}
	return false;
    }


    ConformationExplorerEntryStage_sp	ConformationExplorerEntry_O::getEntryStage(const string& key)
    {_G();
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
    {_G();
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


    core::StringSet_sp ConformationExplorerEntry_O::getEntryStageNames()
    {_G();
	core::StringSet_sp stageNames = core::StringSet_O::create();
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


    void ConformationExplorer_O::clearEntries()
    {_G();
	this->_Entries.clear();
    }




    Atom_sp	ConformationExplorer_O::_getAtomAtIndex(unsigned i)
    {_G();
	ASSERT_lessThan(i,this->_AllAtoms.size());
	return this->_AllAtoms[i];
    }



    geom::CoordinateArray_sp ConformationExplorer_O::_extractCoordinateArray(Matter_sp agg)
    {_G();
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


    void	ConformationExplorer_O::setMatter(Matter_sp matter)
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




    void	ConformationExplorer_O::appendEntry(ConformationExplorerEntry_sp entry)
    {_G();
	entry->setConformationExplorer(this->sharedThis<ConformationExplorer_O>());
	this->_Entries.push_back(entry);
    }


    core::List_sp ConformationExplorer_O::entriesAsCons()
    {_G();
	core::Cons_sp first = core::T_O::create(_Nil<core::T_O>(), _Nil<core::T_O>());
	core::Cons_sp cur = first;
	ConformationExplorer_O::entryIterator si;
	for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
	{
	    core::Cons_sp newCons = core::T_O::create(*si,_Nil<core::T_O>());
	    cur->setCdr(newCons);
	    cur = newCons;
	}
	return first->cdr();
    }









#ifdef XML_ARCHIVE
    void	ConformationExplorer_O::archiveBase(core::ArchiveP node)
    {_G();
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
    {_G();
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
		LOG(BF("The name of the Frame entry is[%s]") % dlEntry->getName()->__repr__() );
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



    void	ConformationExplorer_O::clearSuperposeAtoms()
    {_G();
	ASSERTNOTNULL(this->_SuperposeAtomIndices);
	this->_SuperposeAtomIndices->clear();
    }



    void	ConformationExplorer_O::addSuperposeAtom(Atom_sp a)
    {_G();
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
    {_G();
        gctools::SmallOrderedSet<Atom_sp>::iterator	ai;
	this->clearSuperposeAtoms();
	for ( ai=atoms.begin(); ai!=atoms.end(); ai++ )
	{
	    this->addSuperposeAtom(*ai);
	}
    }

    void	ConformationExplorer_O::superposeAllHeavyAtoms()
    {_G();
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
    {_G();
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
    {_G();
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
	    if ( xValue->operator<(yValue) ) return true;
	    return false;
	}
    };




    void	ConformationExplorer_O::sort(const string& stageName, core::Symbol_sp keyName )
    {_G();
	OrderStageAndKeyValue comparer;
	comparer._StageName = stageName;
	comparer._KeyName = keyName;
	core::sort::quickSort(this->_Entries.begin(), this->_Entries.end(), comparer);
    }







    void	ConformationExplorer_O::saveAs(const string& fn)
    {_G();
	IMPLEMENT_ME();
#ifdef XML_ARCHIVE
	core::XmlSaveArchive_sp	xml;
	xml = core::XmlSaveArchive_O::create();
	xml->put("ConformationExplorer",this->sharedThis<ConformationExplorer_O>());
	xml->saveAs(fn);
#endif
    }

    unsigned ConformationExplorer_O::getEntryIndex(ConformationExplorerEntry_sp entry)
    {_G();
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




    ConformationExplorerEntry_sp	ConformationExplorer_O::getEntry(unsigned i)
    {_G();
	ASSERT_lessThan(i,this->_Entries.size());
	ConformationExplorerEntry_sp e;
	e = this->_Entries[i];
	return e;
    }

    ConformationExplorerEntry_sp	ConformationExplorer_O::firstEntry()
    {_G();
	ASSERT(this->_Entries.size() > 0 );
	return this->_Entries[0];
    }





    ConformationExplorerEntry_sp	ConformationExplorer_O::createEntry()
    {_G();
	GC_ALLOCATE(ConformationExplorerEntry_O, entry );
	entry->setConformationExplorer(this->sharedThis<ConformationExplorer_O>());
	entry->_UniqueEntryIndex = this->numberOfEntries();
	return entry;
    }





    bool	ConformationExplorer_O::hasStageNameInAllEntries(const string& stageKey)
    {_G();
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
    {_G();
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

    ConformationExplorerEntry_sp ConformationExplorer_O::getClosestMatchingEntry(ConformationExplorerMatch_sp match)
    {_G();
	uint idx = match->getClosestMatchIndex();
	return this->getEntry(idx);
    }

    ConformationExplorerEntryStage_sp ConformationExplorer_O::getClosestMatchingEntryStage(ConformationExplorerMatch_sp match)
    {_G();
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
	core::StringSet_sp stageNames = core::StringSet_O::create();
        gctools::Vec0<ConformationExplorerEntry_sp>::const_iterator ei;
	for ( ei=this->_Entries.begin(); ei!=this->_Entries.end(); ei++ )
	{
	    core::StringSet_sp oneEntryStageNames = (*ei)->getEntryStageNames();
	    stageNames->insertStringSet(oneEntryStageNames);
	}
	ss << "EntryStageNames: %s\n", stageNames->asString().c_str();
	return ss.str();
    }

    void ConformationExplorer_O::findClosestMatchingConformation(
	ConformationExplorerMatch_sp match,
	Matter_sp matter, const string& stageKey,
	double rms )
    {_G();
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
    {_G();
	return this->_AllAtoms.indexOf(atom);
    }



    void	ConformationExplorer_O::alignAllConformationsToTheFirstForStage(const string& stageName )
    {_G();
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


    void	ConformationExplorer_O::sortByEnergyAscendingForStage(const string& stageName )
    {_G();
	OrderByEnergyAscending order;
	order._StageName = stageName;
	core::sort::quickSort(this->_Entries.begin(),this->_Entries.end(),order);
    }

    void	ConformationExplorerMatch_O::initialize()
    {
	this->_Matches = false;
    }




    void ConformationExplorerEntryStage_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ConformationExplorerEntryStage_O>()
	    .def("setBasedOnStageName",&ConformationExplorerEntryStage_O::setBasedOnStageName)
	    .def("getBasedOnStageName",&ConformationExplorerEntryStage_O::getBasedOnStageName)
	    .def("setStageName",&ConformationExplorerEntryStage_O::setStageName)
	    .def("getStageName",&ConformationExplorerEntryStage_O::getStageName)
	    .def("getFinalCoordinates",&ConformationExplorerEntryStage_O::getFinalCoordinates)
	    .def("setFinalCoordinatesAsFixedWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setFinalCoordinatesAsFixedWithinSuperposeEngine)
	    .def("setFinalCoordinatesAsMoveableWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setFinalCoordinatesAsMoveableWithinSuperposeEngine)
	    .def("setSuperposableCoordinatesAsFixedWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsFixedWithinSuperposeEngine)
	    .def("setSuperposableCoordinatesAsMoveableWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsMoveableWithinSuperposeEngine)
	    .def("setModel",&ConformationExplorerEntryStage_O::setModel)
	    .def("getModel",&ConformationExplorerEntryStage_O::getModel)
	    .def("setEnergyKCal",&ConformationExplorerEntryStage_O::setEnergyKCal)
	    .def("getEnergyKCal",&ConformationExplorerEntryStage_O::getEnergyKCal)
	    .def("setExternalInterfaceName",&ConformationExplorerEntryStage_O::setExternalInterfaceName)
	    .def("getExternalInterfaceName",&ConformationExplorerEntryStage_O::getExternalInterfaceName)
	    .def("writeCoordinatesToMatter",&ConformationExplorerEntryStage_O::writeCoordinatesToMatter)
	    .def("extractCoordinatesFromMatter",&ConformationExplorerEntryStage_O::extractCoordinatesFromMatter)
	    .def("getBinder",&ConformationExplorerEntryStage_O::getBinder)
	    .def("isComplete",&ConformationExplorerEntryStage_O::isComplete)
	    .def("setComplete",&ConformationExplorerEntryStage_O::setComplete)
	    ;
    }
    void ConformationExplorerEntryStage_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg, ConformationExplorerEntryStage,"","",_lisp)
	.def("setBasedOnStageName",&ConformationExplorerEntryStage_O::setBasedOnStageName)
	.def("getBasedOnStageName",&ConformationExplorerEntryStage_O::getBasedOnStageName)
	.def("setStageName",&ConformationExplorerEntryStage_O::setStageName)
	.def("getStageName",&ConformationExplorerEntryStage_O::getStageName)
	.def("getFinalCoordinates",&ConformationExplorerEntryStage_O::getFinalCoordinates)
	.def("setFinalCoordinatesAsFixedWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setFinalCoordinatesAsFixedWithinSuperposeEngine)
	.def("setFinalCoordinatesAsMoveableWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setFinalCoordinatesAsMoveableWithinSuperposeEngine)
	.def("setSuperposableCoordinatesAsFixedWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsFixedWithinSuperposeEngine)
	.def("setSuperposableCoordinatesAsMoveableWithinSuperposeEngine",&ConformationExplorerEntryStage_O::setSuperposableCoordinatesAsMoveableWithinSuperposeEngine)
	.def("setModel",&ConformationExplorerEntryStage_O::setModel)
	.def("getModel",&ConformationExplorerEntryStage_O::getModel)
	.def("setEnergyKCal",&ConformationExplorerEntryStage_O::setEnergyKCal)
	.def("getEnergyKCal",&ConformationExplorerEntryStage_O::getEnergyKCal)
	.def("setExternalInterfaceName",&ConformationExplorerEntryStage_O::setExternalInterfaceName)
	.def("getExternalInterfaceName",&ConformationExplorerEntryStage_O::getExternalInterfaceName)
    	.def("writeCoordinatesToMatter",&ConformationExplorerEntryStage_O::writeCoordinatesToMatter)
    	.def("extractCoordinatesFromMatter",&ConformationExplorerEntryStage_O::extractCoordinatesFromMatter)
	.def("getBinder",&ConformationExplorerEntryStage_O::getBinder)
	.def("isComplete",&ConformationExplorerEntryStage_O::isComplete)
	.def("setComplete",&ConformationExplorerEntryStage_O::setComplete)
	//	.def("coordinatesAsString",&ConformationExplorerEntryStage_O::coordinatesAsString)
	//	.def("parseStructureFromFileName",&ConformationExplorerEntryStage_O::parseStructureFromFileName)
	;
#endif
}


    void ConformationExplorerEntry_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<ConformationExplorerEntry_O>()
	.def("getUniqueEntryIndex",&ConformationExplorerEntry_O::getUniqueEntryIndex)
    	.def("getConformationExplorer",&ConformationExplorerEntry_O::getConformationExplorer)
    	.def("lastEntryStageName",&ConformationExplorerEntry_O::lastEntryStageName)
    	.def("hasSelectedStage",&ConformationExplorerEntry_O::hasSelectedStage)
    	.def("getSelectedStage",&ConformationExplorerEntry_O::getSelectedStage)
    	.def("setSelectedStage",&ConformationExplorerEntry_O::setSelectedStage)
    	.def("alreadyHasLastCompleteStage",&ConformationExplorerEntry_O::alreadyHasLastCompleteStage)
    	.def("getOrCreateLastIncompleteEntryStage",&ConformationExplorerEntry_O::getOrCreateLastIncompleteEntryStage)
    	.def("hasEntryStageWithName",&ConformationExplorerEntry_O::hasEntryStageWithName)
    	.def("getEntryStage",&ConformationExplorerEntry_O::getEntryStage)
	.def("createEntryStage",&ConformationExplorerEntry_O::createEntryStage)
	.def("getLastCompleteEntryStage",&ConformationExplorerEntry_O::getLastCompleteEntryStage)
	.def("getLastEntryStage",&ConformationExplorerEntry_O::getLastEntryStage)
	;
}
    void ConformationExplorerEntry_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef	USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg, ConformationExplorerEntry,"","",_lisp)
	.add_property("iterate_Stages",
		      boost::python::range(&ConformationExplorerEntry_O::begin_Stages,
					   &ConformationExplorerEntry_O::end_Stages))
    	.def("getConformationExplorer",&ConformationExplorerEntry_O::getConformationExplorer)
    	.def("lastEntryStageName",&ConformationExplorerEntry_O::lastEntryStageName)
    	.def("hasSelectedStage",&ConformationExplorerEntry_O::hasSelectedStage)
    	.def("getSelectedStage",&ConformationExplorerEntry_O::getSelectedStage)
    	.def("setSelectedStage",&ConformationExplorerEntry_O::setSelectedStage)
    	.def("alreadyHasLastCompleteStage",&ConformationExplorerEntry_O::alreadyHasLastCompleteStage)
    	.def("getOrCreateLastIncompleteEntryStage",&ConformationExplorerEntry_O::getOrCreateLastIncompleteEntryStage)
    	.def("hasEntryStageWithName",&ConformationExplorerEntry_O::hasEntryStageWithName)
    	.def("getEntryStage",&ConformationExplorerEntry_O::getEntryStage)
	.def("createEntryStage",&ConformationExplorerEntry_O::createEntryStage)
	.def("getLastCompleteEntryStage",&ConformationExplorerEntry_O::getLastCompleteEntryStage)
	.def("getLastEntryStage",&ConformationExplorerEntry_O::getLastEntryStage)
	;
#endif
}




void ConformationExplorer_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<ConformationExplorer_O>()
//	    .add_property("iterate_AllAtoms",
//		boost::python::range(&ConformationExplorer_O::begin_AllAtoms,
//			&ConformationExplorer_O::end_AllAtoms))
//	    .add_property("iterate_Entries",
//		boost::python::range(&ConformationExplorer_O::begin_Entries,
//			&ConformationExplorer_O::end_Entries))
	.def("entriesAsCons",&ConformationExplorer_O::entriesAsCons)
	.def("getBinder",&ConformationExplorer_O::getBinder)
	.def("setMatter",&ConformationExplorer_O::setMatter)
	.def("getMatter",&ConformationExplorer_O::getMatter)
	.def("firstEntry",&ConformationExplorer_O::firstEntry)
	.def("clearSuperposeAtoms",&ConformationExplorer_O::clearSuperposeAtoms)
	.def("superposeAllHeavyAtoms",&ConformationExplorer_O::superposeAllHeavyAtoms)
	.def("addSuperposeAtom",&ConformationExplorer_O::addSuperposeAtom)
	.def("clearEntries",&ConformationExplorer_O::clearEntries)
	.def("appendEntry",&ConformationExplorer_O::appendEntry)
	.def("getEntry",&ConformationExplorer_O::getEntry)
	.def("createEntry",&ConformationExplorer_O::createEntry)
	.def("getEntryIndex",&ConformationExplorer_O::getEntryIndex)
	.def("saveAs",&ConformationExplorer_O::saveAs,"","","",false)
	.def("numberOfEntries",&ConformationExplorer_O::numberOfEntries)
	.def("findClosestMatchingConformation", &ConformationExplorer_O::findClosestMatchingConformation)
	.def("sort", &ConformationExplorer_O::sort,"","","",false)
	.def("getClosestMatchingEntry", &ConformationExplorer_O::getClosestMatchingEntry)
	.def("getClosestMatchingEntryStage", &ConformationExplorer_O::getClosestMatchingEntryStage)
	.def("alignAllConformationsToTheFirstForStage",&ConformationExplorer_O::alignAllConformationsToTheFirstForStage)
	.def("sortByEnergyAscendingForStage", &ConformationExplorer_O::sortByEnergyAscendingForStage)
	;
}
void ConformationExplorer_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef	USEBOOSTPYTHON
    PYTHON_CLASS(ChemPkg, ConformationExplorer,"","",_lisp)
	.add_property("iterate_AllAtoms",
		      boost::python::range(&ConformationExplorer_O::begin_AllAtoms,
					   &ConformationExplorer_O::end_AllAtoms))
	.add_property("iterate_Entries",
		      boost::python::range(&ConformationExplorer_O::begin_Entries,
					   &ConformationExplorer_O::end_Entries))
	.def("entriesAsCons",&ConformationExplorer_O::entriesAsCons)
	.def("getBinder",&ConformationExplorer_O::getBinder)
	.def("setMatter",&ConformationExplorer_O::setMatter)
	.def("getMatter",&ConformationExplorer_O::getMatter)
	.def("firstEntry",&ConformationExplorer_O::firstEntry)
	.def("clearSuperposeAtoms",&ConformationExplorer_O::clearSuperposeAtoms)
	.def("superposeAllHeavyAtoms",&ConformationExplorer_O::superposeAllHeavyAtoms)
	.def("addSuperposeAtom",&ConformationExplorer_O::addSuperposeAtom)
	.def("clearEntries",&ConformationExplorer_O::clearEntries)
	.def("appendEntry",&ConformationExplorer_O::appendEntry)
	.def("getEntry",&ConformationExplorer_O::getEntry)
	.def("createEntry",&ConformationExplorer_O::createEntry)
	.def("getEntryIndex",&ConformationExplorer_O::getEntryIndex)
//	.def("saveAs",&ConformationExplorer_O::saveAs,"","","",false)
    //	    .def("addConformationExplorer",&ConformationExplorer_O::addConformationExplorer)
    //	    .def("addStructureList",&ConformationExplorer_O::addStructureList)
	.def("numberOfEntries",&ConformationExplorer_O::numberOfEntries)
	.def("findClosestMatchingConformation", &ConformationExplorer_O::findClosestMatchingConformation)
//	.def("sort", &ConformationExplorer_O::sort,"","","",false)
	.def("getClosestMatchingEntry", &ConformationExplorer_O::getClosestMatchingEntry)
	.def("getClosestMatchingEntryStage", &ConformationExplorer_O::getClosestMatchingEntryStage)
	.def("sortByEnergyAscendingForStage", &ConformationExplorer_O::sortByEnergyAscendingForStage)
	;
#endif
}


void ConformationExplorerMatch_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<ConformationExplorerMatch_O>()
	.def("getMatches",&ConformationExplorerMatch_O::getMatches)
	.def("getClosestMatchIndex",&ConformationExplorerMatch_O::getClosestMatchIndex)
	.def("getClosestMatchRms",&ConformationExplorerMatch_O::getClosestMatchRms)
	;
}
void ConformationExplorerMatch_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef	USEBOOSTPYTHON

    PYTHON_CLASS(ChemPkg, ConformationExplorerMatch,"","",_lisp)
	.def("getMatches",&ConformationExplorerMatch_O::getMatches)
	.def("getClosestMatchIndex",&ConformationExplorerMatch_O::getClosestMatchIndex)
	;
#endif
}

EXPOSE_CLASS(chem, ConformationExplorerEntryStage_O );
EXPOSE_CLASS(chem, ConformationExplorerEntry_O );
EXPOSE_CLASS(chem, ConformationExplorer_O );
EXPOSE_CLASS(chem, ConformationExplorerMatch_O );

};




