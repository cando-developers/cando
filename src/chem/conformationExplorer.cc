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
#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	Maintain a database of conformational searches
 *
 */
#include <clasp/core/common.h>
#include <clasp/core/array.h>
#include <cando/adapt/stringSet.h>
#include <cando/geom/coordinateArray.h>
#include <cando/chem/conformationExplorer.h>
#include <cando/chem/loop.h>
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include <clasp/core/hashTableEq.h>
#include <cando/geom/color.h>
#include <cando/chem/superposeEngine.h>
#include <clasp/core/sort.h>
#include <clasp/core/array.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/nVector.h>
#include <clasp/core/numbers.h>
#include <clasp/core/wrappers.h>



namespace chem
{


    void ConformationExplorerEntryStage_O::initialize()
    {
	this->Base::initialize();
	this->_Binder = core::HashTableEq_O::create_default();
	this->_FinalCoordinates = geom::SimpleVectorCoordinate_O::create();
	this->_WeakConformationExplorerEntry = nil<ConformationExplorerEntry_O>();
	this->_Complete = true;
	this->_EnergyKCal = 0.0;
	this->_ExternalInterfaceName = unbound<core::T_O>();
	this->_Model = unbound<core::T_O>();
	this->_StageName = unbound<core::T_O>();
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
	geom__in_place_transform(this->_FinalCoordinates,transform);
    }





#if 0
    geom::Render_sp	ConformationExplorerEntryStage_O::rendered(core::List_sp opts)
    {
	ConformationExplorer_sp		explorer;
	geom::DisplayList_sp			dlAll;
	GrPickableMatter_sp		dlMatter;
	Matter_sp			matter;
	geom::SimpleVectorCoordinate_sp		superposeCoords;
	core::IntArray_O::iterator		ii;
	explorer = this->getConformationExplorer();
	dlAll = geom::DisplayList_O::create();
	FIX_ME(); dlAll->setName(_lisp->internKeyword(this->_StageName));
	{
	    LOG("Writing coordinates to matter");
	    matter = explorer->getMatter();
	    this->writeCoordinatesToMatter(matter);
	    dlMatter = GrPickableMatter_O::create();
	    dlMatter->setFromMatter(matter);
	    dlMatter->setName(_lisp->internKeyword("molFrame"));
	    dlAll->append(dlMatter);
#ifdef DEBUG_ON
	    core::SymbolVector_sp sv = dlMatter->getFullName();
	    LOG("At time of creation dlMatter fullName={}" , _rep_(sv) );
#endif
	}
	{
	    LOG("Rendering superposed atoms as spheres" );
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
	sinfo << "Energy (kcal/mole): " << fixed << setw(10) << setprecision(1) 
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
	geom::SimpleVectorCoordinate_sp	ca;
	ConformationExplorer_sp	explorer;
	this->_WeakConformationExplorerEntry = s;
	explorer = this->getConformationExplorer();
	ca = geom::SimpleVectorCoordinate_O::make(explorer->numberOfAllAtoms());
	this->_FinalCoordinates = ca;
    }


    void	ConformationExplorerEntryStage_O::setFinalCoordinates(geom::SimpleVectorCoordinate_sp ac)
    {
	ASSERTNOTNULL(ac);
	LOG("setFinalCoordinates:{}" , (ac->asXmlString().c_str() ) );
	LOG("The address of the geom::SimpleVectorCoordinate_sp is in o" );
	this->_FinalCoordinates = ac;
    }



    void	ConformationExplorerEntryStage_O::setCoordinateForAtom(Atom_sp atom, Vector3 pos )
    {
	int	idx;
	ConformationExplorer_sp	explorer;
	explorer = this->getConformationExplorer();
	idx = explorer->indexOfAtom(atom);
	ASSERTNOTNULL(this->_FinalCoordinates);
	ASSERT((*this->_FinalCoordinates).length() == explorer->numberOfAllAtoms());
	(*this->_FinalCoordinates)[idx] = pos;
    }


CL_LISPIFY_NAME("writeCoordinatesToMatter");
CL_DEFMETHOD     void	ConformationExplorerEntryStage_O::writeCoordinatesToMatter(Matter_sp agg)
    {
	ConformationExplorer_sp			sl;
        geom::SimpleVectorCoordinate_O::iterator	ci;
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
	LOG("About to get ConformationExplorer" );
	sl = this->getConformationExplorer();
	LOG("Got ConformationExplorer" );
	ASSERTNOTNULL(sl);
	geom::SimpleVectorCoordinate_sp coords = sl->_SimpleVectorCoordinate(matter);
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
	    return nil<ConformationExplorerEntryStage_O>();
	}
	return this->_SelectedStage;
    }

CL_LISPIFY_NAME("setSelectedStage");
CL_DEFMETHOD     void	ConformationExplorerEntry_O::setSelectedStage(ConformationExplorerEntryStage_sp stage)
    {
	ConformationExplorerEntry_sp	stageEntry;
	stageEntry = stage->getConformationExplorerEntry();
	this->_SelectedStage = stage;
    }



CL_LISPIFY_NAME("createEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp	ConformationExplorerEntry_O::createEntryStage(core::T_sp name)
    {
	ASSERTF(!this->hasEntryStageWithName(name), ("Stage with key[{}] already exists!") , name );
	auto  stage  = gctools::GC<ConformationExplorerEntryStage_O>::allocate_with_default_constructor();
	stage->setConformationExplorerEntry(this->sharedThis<ConformationExplorerEntry_O>());
	stage->setStageName(name);
	this->_Stages.push_back(stage);
	return stage;
    }


    void	ConformationExplorerEntry_O::setConformationExplorer(ConformationExplorer_sp s)
    {
	geom::SimpleVectorCoordinate_sp	ca;
	this->_WeakConformationExplorer = s;
    }




CL_LISPIFY_NAME("lastEntryStageName");
CL_DEFMETHOD     core::T_sp	ConformationExplorerEntry_O::lastEntryStageName()
    {
      if ( this->_Stages.size() == 0 ) return nil<core::T_O>();
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
CL_DEFMETHOD     bool	ConformationExplorerEntry_O::alreadyHasLastCompleteStage(core::T_sp stageName )
    {
	ConformationExplorerEntryStage_sp	lastStage;
	lastStage = this->getLastCompleteEntryStage();
	return( lastStage->getStageName() == stageName );
    }


CL_LISPIFY_NAME("getOrCreateLastIncompleteEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp ConformationExplorerEntry_O::getOrCreateLastIncompleteEntryStage(core::T_sp stageName)
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
		newStage = nil<ConformationExplorerEntryStage_O>();
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
	SIMPLE_ERROR("Only the last stage can be incomplete");
    }







CL_LISPIFY_NAME("hasEntryStageWithName");
CL_DEFMETHOD     bool	ConformationExplorerEntry_O::hasEntryStageWithName(core::T_sp key)
    {
	stageIterator	si;
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{
	    if ( (*si)->getStageName() == key ) return true;
	}
	return false;
    }


CL_LISPIFY_NAME("getEntryStage");
CL_DEFMETHOD     ConformationExplorerEntryStage_sp	ConformationExplorerEntry_O::getEntryStage(core::T_sp key)
    {
	stageIterator	si;
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{
	    if ( (*si)->getStageName() == key )
	    {
		return (*si);
	    }
	}
	SIMPLE_ERROR("Could not find key: {}" , core::_rep_(key));
    }



#if 0
    geom::Render_sp	ConformationExplorerEntry_O::rendered(core::List_sp opts)
    {
	geom::DisplayList_sp	dlAll;
	geom::Render_sp		dlStage;
	stageIterator		si;
	dlAll = geom::DisplayList_O::create();
	FIX_ME(); dlAll->setName(_lisp->internKeyword("stages"));
	LOG("Rendering {} stages"%this->_Stages.size());
	for ( si=this->_Stages.begin(); si!=this->_Stages.end(); si++ )
	{
	    if ( (*si)->isComplete() )
	    {
		LOG("It was complete");
		dlStage = (*si)->rendered(opts);
		dlAll->append(dlStage);
	    } else
	    {
		LOG("It was NOT complete");
	    }
	}
	return dlAll;
    }
#endif


core::HashTable_sp ConformationExplorerEntry_O::getEntryStageNames()
{
  core::HashTableEq_sp stageNames = core::HashTableEq_O::create_default();
  for ( auto it=this->begin_Stages(); it!=this->end_Stages(); it++ ) {
    stageNames->setf_gethash((*it)->getStageName(), _lisp->_true());
  }
  return stageNames;
}

    void	ConformationExplorer_O::initialize()
    {
	this->Base::initialize();
	this->_Matter = nil<Matter_O>();
	this->_AllAtoms.clear();
	this->_Binder = core::HashTableEq_O::create_default();
	this->clearEntries();
	this->_SuperposeAtomIndices = core::ComplexVector_byte32_t_O::make_vector(0,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
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



    geom::SimpleVectorCoordinate_sp ConformationExplorer_O::_SimpleVectorCoordinate(Matter_sp agg)
    {
      geom::SimpleVectorCoordinate_O::iterator	ci;
        gctools::SmallOrderedSet<Atom_sp>::iterator		ai;
#ifdef	DEBUG_ConformationExplorerEntry
	this->_Status->addMessage("extractCoordinatesFromMatter");
#endif
	LOG("About to get ConformationExplorer");
	LOG("About to iterate over atoms and get positions");
	uint numberOfAtoms = agg->numberOfAtoms();
	geom::SimpleVectorCoordinate_sp coords = geom::SimpleVectorCoordinate_O::make(numberOfAtoms);
	ASSERTNOTNULL(coords);
	uint aidx = 0;
	for ( ai=this->begin_AllAtoms(),ci=coords->begin(); ai!=this->end_AllAtoms(); ai++, ci++ )
	{
	    LOG("Getting position of ATOM({}) {}" , aidx , (*ai)->description().c_str() );
	    *ci = (*ai)->getPosition();
	    LOG("Position = {}" , (*ci).asString() );
	    aidx++;
	}
	LOG("Done");
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
	core::Cons_sp first = core::Cons_O::create(nil<core::T_O>(), nil<core::T_O>());
	core::Cons_sp cur = first;
	ConformationExplorer_O::entryIterator si;
	for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
	{
	    core::Cons_sp newCons = core::Cons_O::create(*si,nil<core::T_O>());
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
	geom::SimpleVectorCoordinate_sp		superposeCoords;
	geom::SimpleVectorCoordinate_O::iterator	ci;
	_LINE();
	geom::DisplayList_sp dlAll = geom::DisplayList_O::create();
	_LINE();
	frames = geom::FrameList_O::create();
	FIX_ME(); frames->setName(_lisp->internKeyword("conformation"));
	_LINE();
	if ( this->numberOfEntries() == 0 )
	{
	    LOG("FLAG! There are no entries");
	    ASSERTNOTNULL(this->_Matter);
	    if ( this->_Matter.notnilp() )
	    {
		LOG("Rendering ConformationExplorer::matter");
		dlEntry = this->_Matter->rendered(opts);
		ASSERTNOTNULL(dlEntry);
		frames->append(dlEntry);
	    }
	} else
	{
	    LOG("FLAG! There are {} entries" , this->_Entries.size());
	    for ( si=this->begin_Entries(); si!=this->end_Entries(); si++ )
	    {
		dlEntry = (*si)->rendered(opts);
		LOG("The name of the Frame entry is[{}]" , _rep_(dlEntry->getName()) );
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
	this->_SuperposeAtomIndices = core::ComplexVector_byte32_t_O::make_vector(0,0,core::make_fixnum(0),nil<core::T_O>(),false,core::make_fixnum(0));
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
              this->_SuperposeAtomIndices->vectorPushExtend(idx);
		return;
	    }
	}
	SIMPLE_ERROR("You cannot add superpose atom: {} because it isn't in the ConformationCollection", a->description());
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
		this->_SuperposeAtomIndices->vectorPushExtend(idx);
	    }
	}
    }


    gctools::SmallOrderedSet<Atom_sp>	ConformationExplorer_O::getSuperposeAtoms()
    {
        gctools::SmallOrderedSet<Atom_sp>	result;
	Atom_sp			a;
	ASSERTNOTNULL(this->_SuperposeAtomIndices);
	for ( size_t si=0; si<this->_SuperposeAtomIndices->length(); si++ )
	{
          a = this->_AllAtoms[(*this->_SuperposeAtomIndices)[si]];
          result.insert(a);
	}
	return result;
    }


    uint	ConformationExplorer_O::numberOfSuperposeAtoms()
    {
	return this->_SuperposeAtomIndices->length();
    }


    struct OrderStageAndKeyValue
    {
      core::T_sp _StageName;
      core::Symbol_sp _KeyName;
      
      bool operator()(ConformationExplorerEntry_sp x,ConformationExplorerEntry_sp y )
      {
        ConformationExplorerEntryStage_sp xStage, yStage;
        xStage = x->getEntryStage(this->_StageName);
        yStage = y->getEntryStage(this->_StageName);
        core::T_sp xValue = xStage->getBinder()->gethash(this->_KeyName);
        core::T_sp yValue = yStage->getBinder()->gethash(this->_KeyName);
        if ( core::clasp_number_compare(xValue.as<core::Number_O>(),yValue.as<core::Number_O>())) return true;
        return false;
      }
    };




CL_LISPIFY_NAME("ConformationExplorer-sort");
CL_DEFMETHOD     void	ConformationExplorer_O::sort(core::T_sp stageName, core::Symbol_sp keyName )
    {
	OrderStageAndKeyValue comparer;
	comparer._StageName = stageName;
	comparer._KeyName = keyName;
	sort::quickSortVec0(this->_Entries, 0, this->_Entries.size(), comparer);
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
	SIMPLE_ERROR("Could not find entry in ConformationExplorer");
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
      auto  entry  = gctools::GC<ConformationExplorerEntry_O>::allocate_with_default_constructor();
	entry->setConformationExplorer(this->sharedThis<ConformationExplorer_O>());
	entry->_UniqueEntryIndex = this->numberOfEntries();
	return entry;
    }





bool	ConformationExplorer_O::hasStageNameInAllEntries(core::T_sp stageKey)
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
	core::T_sp	stageName,
	double&		bestRms,
	ConformationExplorerEntryStage_sp& bestStage,
	uint&		bestEntryIndex )
    {
	ConformationExplorerEntryStage_sp	stage;
	double					rms;
	bool					gotBest;
	SuperposeEngine_sp			superposer;
	geom::SimpleVectorCoordinate_sp			matterConf;
	bestStage = nil<ConformationExplorerEntryStage_O>();
	bestRms = 9.9e99;
	gotBest = false;
	uint numSuperposeAtoms = this->numberOfSuperposeAtoms();
	if ( numSuperposeAtoms == 0 )
	{
	    SIMPLE_ERROR("There must be at least one superpose atom");
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
	    matterConf = this->_SimpleVectorCoordinate(matter);
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

CL_LISPIFY_NAME("findClosestMatchingConformation");
CL_DEFMETHOD     void ConformationExplorer_O::findClosestMatchingConformation(
	ConformationExplorerMatch_sp match,
	Matter_sp matter,
        core::T_sp stageKey,
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
	    LOG("There were no entries to match" );
	    return;
	}
	if ( bestRms < rms )
	{
	    LOG("The best match had an rms({}) and we are looking for one better than({})" , bestRms , rms  );
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
CL_DEFMETHOD     void	ConformationExplorer_O::alignAllConformationsToTheFirstForStage(core::T_sp stageName )
    {
	if ( !this->hasStageNameInAllEntries(stageName) )
	{
          SIMPLE_ERROR("You have to have the stageName({}) in all entries" , stageName);
	}
	if ( this->_Entries.size() == 0 )
	{
	    SIMPLE_ERROR("You have to have at least one entry");
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
      core::T_sp	_StageName;
    public:
	bool operator()(ConformationExplorerEntry_sp x, ConformationExplorerEntry_sp y )
	{
	    ConformationExplorerEntryStage_sp xStage = x->getEntryStage(this->_StageName);
	    ConformationExplorerEntryStage_sp yStage = y->getEntryStage(this->_StageName);
	    return xStage->getEnergyKCal()<yStage->getEnergyKCal();
	}
    };


CL_LISPIFY_NAME("sortByEnergyAscendingForStage");
CL_DEFMETHOD     void	ConformationExplorer_O::sortByEnergyAscendingForStage(core::T_sp stageName )
    {
	OrderByEnergyAscending order;
	order._StageName = stageName;
	sort::quickSortVec0(this->_Entries, 0, this->_Entries.size(), order);
    }

    void	ConformationExplorerMatch_O::initialize()
    {
	this->_Matches = false;
    }


















};
