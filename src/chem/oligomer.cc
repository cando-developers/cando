#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//

#include <cando/main/foundation.h>
#include <cando/chem/oligomer.h>
#include <clasp/core/str.h>
#include <cando/chem/molecule.h>
#include <clasp/core/numerics.h>
#include <cando/chem/candoScript.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/binder.h>
#include <cando/chem/topology.h>
#include <clasp/core/environment.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/monomer.h>
#include <cando/chem/coupling.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/plug.h>
#include <cando/chem/specificContext.h>
#include <clasp/core/wrappers.h>

/*! \page oligomerBuilding  How to build oligomers
 * Building oligomers requires the following steps.
 * 	-# Define the oligomer
 * 		- Create the _Monomers, add them to the oligomer and couple them together.
 * 	-# Construct the alchemists Oligomer_O::createAlchemists (you may need to dismissAlchemists first)
 * 	-# Either accept the default build plan or construct one
 *	 	-# For each monomer define the domainId using setDomainId
 * 		-# (optional) Fix the coordinates of the fragmentAlchemists that have defined coordinates
 * 			-# Clear all Fixed alchemists using Oligomer_O::clearAllFragmentAlchemistsCoordinatesAreFixed
 * 			-# Use setAllFragemntAlchemistsCoordinatesAreFixed(True)
 * 			-# Define the coordinates for the residue of the monomer using applyCoordinatesFromResidue(res)
 * 			-# Create a BuildPlan Dumb_BuildPlan::Dumb_BuildPlan(_kinemat_)
 * 			-# Define a BuildPlan using Dumb_BuildPlan::recordPlan()
 * 			-# Tell the Conformation that you want to use that plan using Conformation_O::applyBuildPlan
 * 	-# Build the coordinates by calling:
 * 		- Conformation_O::build
 * 		- or Conformation_O::buildCore followed by Conformation_O::buildFinish
 * 		- or Oligomer_O::buildCoordinates()
 */




namespace chem {





void Oligomer_O::initialize()
{_OF();
    this->Base::initialize();
    this->empty();
    this->_Verbose = false;
    this->_HasError = false;
    this->_Name = _Nil<NameType>();
    this->_ErrorMessage = core::clasp_make_string_output_stream();
}


#if 0
void	Oligomer_O::signal_monomerContentsChanged()
{_G();
    this->notify(Oligomer_monomerContentsChanged);
}


void	Oligomer_O::signalConnectivityChanged()
{_G();
//    this->notify(Oligomer_connectivityChanged);
}
#endif


    void Oligomer_O::catchSignal(core::Symbol_sp signal, core::Model_sp source, core::List_sp data)
{_G();
    MultiMonomer_sp		monomer;
//    this->core::Model_O::catchSignal(signal,source,data);
    IMPLEMENT_ME();
#if 0
    switch ( signal )
    {
	case MultiMonomer_contentsChanged:
	    monomer = downcast<MultiMonomer_O>(source);
	    this->checkForErrors();
	    this->_WeakLastMultiMonomerChanged = monomer;
	    break;
    }
#endif
}

CL_LISPIFY_NAME("hasLastMultiMonomerChanged");
CL_DEFMETHOD bool	Oligomer_O::hasLastMultiMonomerChanged()
{_OF();
    ASSERTNOTNULL(this->_WeakLastMultiMonomerChanged);
    return this->_WeakLastMultiMonomerChanged.notnilp();
}
CL_LISPIFY_NAME("getLastMultiMonomerChanged");
CL_DEFMETHOD MultiMonomer_sp Oligomer_O::getLastMultiMonomerChanged()
{_OF();
    ASSERTNOTNULL(this->_WeakLastMultiMonomerChanged);
    return this->_WeakLastMultiMonomerChanged;
};


void	Oligomer_O::propagateSignal(core::Symbol_sp signal)
{_G();
	// Do nothing
}




CL_LISPIFY_NAME("empty");
CL_DEFMETHOD void	Oligomer_O::empty()
{
    this->_Monomers.clear();
    this->_Couplings.clear();
}

#if 0 //[
void Oligomer_O::setCandoDatabase(CandoDatabase_sp bdb)
{
    this->_WeakCandoDatabase = bdb->sharedThis<CandoDatabase_O>();
};
#endif //]

#ifdef XML_ARCHIVE
    void	Oligomer_O::archiveBase(core::ArchiveP node)
{_G();
//    node->archiveWeakPointer("candoDatabase",this->_WeakCandoDatabase);
    this->Base::archiveBase(node);
    node->attribute("name",this->_Name);
    node->archiveVector0("monomers",this->_Monomers);
    node->archiveVector0("couplings",this->_Couplings);
    node->attribute("hasError",this->_HasError);
    node->attribute("verbose",this->_Verbose);
    node->attribute("errorMessage",this->_ErrorMessage);
}
#endif


CL_LISPIFY_NAME("monomersAsCons");
CL_DEFMETHOD core::List_sp Oligomer_O::monomersAsCons()
{_G();
    core::List_sp cons = _Nil<core::T_O>();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	cons = core::Cons_O::create(*mi,cons);
    }
    return cons;
}

    void	Oligomer_O::expandMonomerListToNeighbors(gctools::SmallOrderedSet<Monomer_sp>& monomers)
{_G();
    gctools::SmallOrderedSet<Monomer_sp>				expanded;
    gctools::SmallOrderedSet<Monomer_sp>::iterator		mi;
Monomer_sp				neighbor;
Monomer_O::Couplings::const_iterator	ci;
Coupling_sp				coupling;
    expanded.clear();
    LOG(BF("Starting set of monomers:") );
    for ( mi=monomers.begin(); mi!=monomers.end(); mi++ )
    {
	LOG(BF("    %s") % (*mi)->description().c_str()  );
        for ( ci = (*mi)->_Couplings.begin();
		ci != (*mi)->_Couplings.end(); ci++ )
	{
	    ASSERTNOTNULL(ci->second);
	    coupling = ci->second;
	    neighbor = coupling->getOtherSideMonomer(*mi);
	    expanded.insert(neighbor);
	}
    }
    monomers.clear();
    LOG(BF("Expanded set of monomers:") );
    for ( gctools::SmallOrderedSet<Monomer_sp>::iterator si=expanded.begin(); si!=expanded.end(); si++ )
    {
	LOG(BF("    %s") % (*si)->description().c_str()  );
        monomers.insert(*si);
    };
}



    void	Oligomer_O::checkMonomersAndNeighborsForErrors(CandoDatabase_sp cdb, gctools::SmallOrderedSet<Monomer_sp>& monomers)
{_G();
    gctools::SmallOrderedSet<Monomer_sp>			expanded;
    ASSERTNOTNULL(cdb);
    expanded = monomers;
    this->expandMonomerListToNeighbors(expanded);
    this->checkMonomersAndNotNeighborsForErrors(cdb,expanded);
}


    void	Oligomer_O::checkMonomersAndNotNeighborsForErrors(CandoDatabase_sp cdb, gctools::SmallOrderedSet<Monomer_sp>& monomers)
{_G();
    gctools::SmallOrderedSet<Monomer_sp>::iterator	it;
    ASSERTNOTNULL(cdb);
    for ( it=monomers.begin(); it!=monomers.end(); it++ )
    {
        (*it)->checkForErrorsAndUnknownContexts(cdb);
    }
    this->checkForErrors();
}



/*! Add a monomer, don't change its id
 * Oligomer._State must be (alchemistsAreUndefined)
 */
CL_LISPIFY_NAME("addMonomer");
CL_DEFMETHOD uint Oligomer_O::addMonomer(Monomer_sp s)
{_G();
Monomer_sp	monomer;
uint		idx;
    monomer = s->sharedThis<Monomer_O>();
    LOG(BF("s@%p    monomer@%X") % s.get() % monomer.get() );
    LOG(BF("Adding monomer@%X with name (%s)") % monomer.get() % s->description().c_str() );
    idx = this->_Monomers.size();
    monomer->setSequenceNumber(idx);
    this->_Monomers.push_back(monomer);
    if ( monomer->getHasError() )
    {
        this->_HasError = true;
    }
//    monomer->connectListener(this->sharedThis<Oligomer_O>(),MultiMonomer_contentsChanged);
    monomer->setOligomer(this->sharedThis<Oligomer_O>());
    return idx;
}

CL_LISPIFY_NAME("hasMonomerWithId");
CL_DEFMETHOD bool Oligomer_O::hasMonomerWithId(core::Symbol_sp id)
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	if ( (*mi)->getId() == id ) return true;
    }
    return false;
}


CL_LISPIFY_NAME("getMonomerWithId");
CL_DEFMETHOD Monomer_sp Oligomer_O::getMonomerWithId(core::Symbol_sp id)
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	if ( (*mi)->getId() == id ) return *mi;
    }
    return _Nil<Monomer_O>();
}




/*! Remove a monomer.
Oligomer._State must be (alchemistsAreUndefined) */
CL_LISPIFY_NAME("removeMonomer");
CL_DEFMETHOD void	Oligomer_O::removeMonomer(Monomer_sp s)
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	if ( *mi == s ) 
	{
//	    (*mi)->disconnectAll(this->sharedThis<Oligomer_O>());
	    this->_Monomers.erase(mi);
	    break;
	}
    }
}


/*! Remove a coupling.
Oligomer._State must be (alchemistsAreUndefined) */
CL_LISPIFY_NAME("Oligomer-removeCoupling");
CL_DEFMETHOD void	Oligomer_O::removeCoupling(Coupling_sp s)
{
    gctools::Vec0<Coupling_sp>::iterator	mi;
    for ( mi=this->_Couplings.begin(); mi!=this->_Couplings.end(); mi++ ) {
	if ( *mi == s ) {
//	    (*mi)->disconnectAll(this->sharedThis<Oligomer_O>());
	    this->_Couplings.erase(mi);
	    break;
	}
    }
}



CL_LISPIFY_NAME("checkForErrors");
CL_DEFMETHOD bool	Oligomer_O::checkForErrors()
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    gctools::Vec0<Coupling_sp>::iterator	ci;
    this->_HasError = false;
    this->_ErrorMessage = core::clasp_make_string_output_stream();
    for ( mi=this->_Monomers.begin(); mi != this->_Monomers.end(); mi++ )
    {
	if ( (*mi)->getHasError() )
	{
	    this->addErrorMessage("Error in "+(*mi)->description() );
	}
    }
    for ( ci=this->_Couplings.begin();
		ci != this->_Couplings.end(); ci++ )
    {
	if ( (*ci)->getHasError() )
	{
	    this->_HasError = true;
	    this->addErrorMessage("Error in "+(*ci)->description() );
	}
    }
    return this->_HasError;
}



CL_LISPIFY_NAME("updateMultiMonomers");
CL_DEFMETHOD void	Oligomer_O::updateMultiMonomers()
{_G();
Oligomer_O::monomerIterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        (*mi)->updateGroupName();
    }
}



/*! Add a Monomer as a leaf to the indicated monomer
 */
Monomer_sp	Oligomer_O::addMultiMonomerAsLeaf( CandoDatabase_sp db,
						   Monomer_sp monFrom, core::Symbol_sp couplingName,
						   core::Symbol_sp representedEntityNameSetName, geom::Vector2 pos )
{_G();
DirectionalCoupling_sp 	coup;
MultiMonomer_sp	monTo;

    coup = DirectionalCoupling_O::create();
    coup->setName(couplingName);
    monTo = MultiMonomer_O::create();
    monTo->setGroupName(representedEntityNameSetName);
    monTo->setPosition2D(pos);
    coup->setOutMonomer(monTo);
    monTo->setInCoupling(coup);
    	// Now attach it
    monFrom->addOutCoupling(coup);
    coup->setInMonomer(monFrom);
    this->addMonomer(monTo);
    this->addCoupling(coup);
    monFrom->throwIfBadConnections();
    monTo->throwIfBadConnections();
    coup->throwIfBadConnections();
    return monTo;
}


/*! Remove a Monomer that is a leaf.
 * Returns the monomer that it was attached to.
 */
CL_LISPIFY_NAME("removeLeafMonomer");
CL_DEFMETHOD Monomer_sp Oligomer_O::removeLeafMonomer(Monomer_sp monomerToRemove)
{_G();
uint	couplings;
DirectionalCoupling_sp	couplingToRemove;
Monomer_sp	monomerToRemain;

    if ( !monomerToRemove->hasInCoupling() )
    {
        SIMPLE_ERROR(BF("Trying to remove a leaf monomer that has no in coupling!") %monomerToRemove->description());
    }
    couplings = monomerToRemove->numberOfCouplings();
    if ( couplings > 1 )
    {
        SIMPLE_ERROR(BF("Trying to remove a leaf monomer but it has out couplings!")%monomerToRemove->description());
    }
    couplingToRemove = monomerToRemove->getInCoupling().as<DirectionalCoupling_O>();
    monomerToRemain = couplingToRemove->getMonomer1();
    monomerToRemain->removeCoupling(couplingToRemove);
    this->removeMonomer(monomerToRemove);
    this->removeCoupling(couplingToRemove);
    monomerToRemain->throwIfBadConnections();
    return monomerToRemain;
}



/*! Add a coupling.
Oligomer._State must be (alchemistsAreUndefined) */
CL_LISPIFY_NAME("Oligomer-addCoupling");
CL_DEFMETHOD void	Oligomer_O::addCoupling(Coupling_sp c)
{_G();
    this->_Couplings.push_back(c);
//    c->connectListener(this->sharedThis<Oligomer_O>(),Coupling_nameChanged);
    c->setOligomer(this->sharedThis<Oligomer_O>());
}



/*! Return the monomer at the root of the tree.
Oligomer._State can be anything. */
CL_LISPIFY_NAME("rootMonomer");
CL_DEFMETHOD Monomer_sp	Oligomer_O::rootMonomer()
{_G();
Monomer_sp	sub;
DirectionalCoupling_sp	coupling;
	// _State can be anything
coupling = _Nil<DirectionalCoupling_O>();
    sub = this->_Monomers[0];
    do {
        LOG(BF("Looking at monomer: %s") % sub->description().c_str()  );
	if ( sub->hasInCoupling() )
	{
	    LOG(BF("    It has an in-coupling") );
	    coupling = sub->getInCoupling();
	    if ( coupling.nilp() ) 
	    {
	        LOG(BF("     But its Nil") );
		break;
	    }
	    sub = coupling->getInMonomer();
	} else
	{
	    coupling = _Nil<DirectionalCoupling_O>();
	}
	ANN(coupling);
    } while ( coupling.notnilp() );
    return sub;
}


CL_LISPIFY_NAME("throwIfBadConnections");
CL_DEFMETHOD void	Oligomer_O::throwIfBadConnections()
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    gctools::Vec0<Coupling_sp>::iterator	ci;
    for ( mi=this->_Monomers.begin(); mi != this->_Monomers.end(); mi++ )
    {
        (*mi)->throwIfBadConnections();
    }
    for ( ci=this->_Couplings.begin();
		ci != this->_Couplings.end(); ci++ )
    {
        (*ci)->throwIfBadConnections();
    }
}








CL_LISPIFY_NAME("setVerbose");
CL_DEFMETHOD void	Oligomer_O::setVerbose(bool v)
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    this->_Verbose = v;
    LOG(BF("Setting verbosity of monomers") );
    for ( mi = this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) {
        ASSERTNOTNULLP(*mi,"Monomer is undefined!!!");
        LOG(BF("Setting verbosity of monomer: %s") % (*mi)->description().c_str() );
	(*mi)->setVerbose(v);
    }
}





CL_LISPIFY_NAME("couple");
CL_DEFMETHOD DirectionalCoupling_sp	Oligomer_O::couple( Monomer_sp inMon, core::Symbol_sp name, Monomer_sp outMon )
{_OF();
DirectionalCoupling_sp	coupling;
bool		foundIn, foundOut;
gctools::Vec0<Monomer_sp>::iterator	mi;
    foundIn = false;
    foundOut = false;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	if ( *mi == inMon ) {
	    foundIn = true;
	}
	if ( *mi == outMon ) {
	    foundOut = true;
	}
    }
    if ( !foundIn ) {
	SIMPLE_ERROR(BF("Could not find first monomer in oligomer"));
    }
    if ( !foundOut ) {
	SIMPLE_ERROR(BF("Could not find second monomer in oligomer"));
    }
    core::Symbol_sp outMonPlugName = DirectionalCoupling_O::inPlugName(name);
    if ( outMon->hasCouplingWithPlugName(outMonPlugName) ) {
	SIMPLE_ERROR(BF("The second monomer already has an in coupling"));
    }
    core::Symbol_sp inMonPlugName = DirectionalCoupling_O::outPlugName(name);
    if ( inMon->hasCouplingWithPlugName(inMonPlugName) ) {
	SIMPLE_ERROR(BF("The first monomer[%s] already has an out coupling with the name: %s") % inMon->getName() % name );
    }
    coupling = DirectionalCoupling_O::create();
    coupling->setName(name);
    coupling->setInMonomer(inMon->sharedThis<Monomer_O>());
    coupling->setOutMonomer(outMon->sharedThis<Monomer_O>());
    inMon->addCoupling(inMonPlugName,coupling);
    outMon->addCoupling(outMonPlugName,coupling);
    this->addCoupling(coupling);
    LOG(BF("after coupling coupling: %s") % coupling->description().c_str()  );
    return coupling;
}


CL_LISPIFY_NAME("ringCouple");
CL_DEFMETHOD RingCoupling_sp	Oligomer_O::ringCouple( Monomer_sp mon1, Monomer_sp mon2 )
{_OF();
RingCoupling_sp	coupling;
bool		found1, found2;
gctools::Vec0<Monomer_sp>::iterator	mi;
    found1 = false;
    found2 = false;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	if ( *mi == mon1 ) {
	    found1 = true;
	}
	if ( *mi == mon2 ) {
	    found2 = true;
	}
    }
    if ( !found1 ) 
    {
	SIMPLE_ERROR(BF("Could not find monomer1 in oligomer"));
    }
    if ( !found2 ) 
    {
	SIMPLE_ERROR(BF("Could not find monomer2 in oligomer"));
    }
    RingClosingPlug_sp rcPlug1 = mon1->getMissingRingClosingPlug(mon2);
    if ( rcPlug1.nilp() )
    {
	SIMPLE_ERROR(BF("The %s can not make a ring closing coupling")
					 % mon1->description() );
    }
    RingClosingPlug_sp rcPlug2 = mon2->getMissingRingClosingPlug(mon1);
    if ( rcPlug2.nilp() )
    {
	SIMPLE_ERROR(BF("The %s can not make a ring closing coupling")
					 % mon2->description()
		  );
    }
    core::Symbol_sp mon1PlugName = rcPlug1->getName();
    core::Symbol_sp mon2PlugName = rcPlug2->getName();
    coupling = RingCoupling_O::create();
    coupling->setMonomer1(mon1->sharedThis<Monomer_O>());
    coupling->setPlug1(mon1PlugName);
    coupling->setMonomer2(mon2->sharedThis<Monomer_O>());
    coupling->setPlug2(mon2PlugName);
    mon1->addCoupling(mon1PlugName,coupling);
    mon2->addCoupling(mon2PlugName,coupling);
    this->addCoupling(coupling);
    return coupling;
}




CL_LISPIFY_NAME("getMolecule");
CL_DEFMETHOD Molecule_sp Oligomer_O::getMolecule()
{_G();
Topology_sp		topology;
int			residueNetCharge;
    if ( this->getHasError() ) {
	SIMPLE_ERROR(BF("The molecule can not be returned, the Oligomer has an error"));
    }
    Residue_sp			res;
    Molecule_sp 			mol;
    mol = Molecule_O::create();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	ASSERT((*mi).notnilp());
    	res = (*mi)->createResidue();
	topology = (*mi)->getTopology();
	residueNetCharge = topology->getResidueNetCharge();
	res->setNetCharge(residueNetCharge);
	(*mi)->setTempResidue(res);
	mol->addMatter(res);
    }
    gctools::Vec0<Coupling_sp>::iterator	ci;
    Monomer_sp	mon1, mon2;
    Residue_sp	res1, res2;
    for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ ) 
    {
	mon1 = (*ci)->getMonomer1();
	ASSERT(mon1.notnilp());
	mon2 = (*ci)->getMonomer2();
	ASSERT(mon2.notnilp());
	res1 = mon1->getTemporaryResidue();
	res2 = mon2->getTemporaryResidue();
    	(*ci)->doCoupling(res1,res2);
    };
    return mol;
}


CL_LISPIFY_NAME("getErrorMessage");
CL_DEFMETHOD string	Oligomer_O::getErrorMessage()
{_G();
stringstream			msg;
gctools::Vec0<Monomer_sp>::iterator	mi;
gctools::Vec0<Coupling_sp>::iterator	ci;
    msg.str("");
    msg << core::cl__get_output_stream_string(this->_ErrorMessage).as<core::Str_O>()->get();
    msg << "Monomer errors: " << std::endl;
    LOG(BF("Adding monomer errors") );
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) {
	if ( (*mi)->getHasError() ) {
	    msg << (*mi)->getStatusMessage() << std::endl;
	}
    }
    LOG(BF("Adding coupling errors") );
    for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ ) {
	if ( (*ci)->getHasError() ) {
	    msg << (*mi)->getStatusMessage() << std::endl;
	}
    }
    LOG(BF("Done adding errors") );
    return msg.str();
}


void Oligomer_O::_gatherMultiMonomers(gctools::Vec0<Monomer_sp>& multiMonomers)
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    multiMonomers.clear();
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	if ( (*mi)->numberOfPossibleMonomers()>1 )
	{
	    multiMonomers.push_back(*mi);
	}
    }
    if ( multiMonomers.size() == 0 )
    {
	LOG(BF("There are no Monomers with more than 1 possibles") );
    }
}



/*!
 * Perturb one MultiMonomer in the Oligomer.
 * If there are none then do nothing.
 */
void	Oligomer_O::perturbMonomerSequence()
{_G();
    gctools::Vec0<Monomer_sp>	multiMonomers;
uint				iperturb;
    this->_gatherMultiMonomers(multiMonomers);
    	// 
	// If there are no MultiMonomers then just return
	//
    if ( multiMonomers.size() == 0 ) 
    {
	//
	// For now throw an exception
	//
	SIMPLE_ERROR(BF("There are no multimonomers to perturb"));
    }
    iperturb = core::randomNumber01()*multiMonomers.size();
    LOG(BF("Randomizing MultiMonomer to %d out of %d") % iperturb % multiMonomers.size());
    ASSERT_lt(iperturb,multiMonomers.size());
    multiMonomers[iperturb]->randomizeMonomer();
}


void	Oligomer_O::gotoSequence(const Bignum& idx)
{_G();
vector<int>	digits, bases;
gctools::Vec0<Monomer_sp>::iterator	mi;
vector<int>::iterator	di;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        bases.push_back((*mi)->numberOfPossibleMonomers());
    }
    digits = core::bignumToMixedBaseDigits(idx,bases);
    for ( mi=this->_Monomers.begin(),di=digits.begin();
    		mi!=this->_Monomers.end(); mi++,di++ )
    {
        (*mi)->setMonomerIndex(*di);
    }
}


/*!
 * Set all MultiMonomers to their first monomer.
 */
void	Oligomer_O::firstSequence()
{_G();
    gctools::Vec0<Monomer_sp>	multiMonomers;
    this->gotoSequence(0);
}


/*! Convert the current sequence into an index
 */
Bignum Oligomer_O::currentSequenceIndex()
{_G();
vector<int>	digits, bases;
gctools::Vec0<Monomer_sp>::iterator	mi;
Bignum			index;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        bases.push_back((*mi)->numberOfPossibleMonomers());
	digits.push_back((*mi)->getMonomerIndex());
    }
    index = core::mixedBaseDigitsToBignum(bases,digits);
    return index;
}


/*!
 * Increment the monomer sequence.
 * Return -1 if there are too many sequences.
 */
CL_LISPIFY_NAME("numberOfSequences");
CL_DEFMETHOD Bignum Oligomer_O::numberOfSequences()
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
vector<int>			bases;
Bignum				numSeq;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        bases.push_back((*mi)->numberOfPossibleMonomers());
    }
    numSeq = core::numberOfIndicesForMixedBase(bases);
    return numSeq;
}

/*!
 * Increment the monomer sequence.
 * Return true if it incremented, return false if it overflowed
 */
bool Oligomer_O::incrementSequence()
{_G();
    gctools::Vec0<Monomer_sp>::iterator	mi;
long unsigned int numSeq;
    numSeq = 1;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        if ( (*mi)->incrementMonomerIndex() ) return true;
        numSeq *= (*mi)->numberOfPossibleMonomers();
    }
    return false;
}




CL_LISPIFY_NAME("getFirstMonomer");
CL_DEFMETHOD Monomer_sp	Oligomer_O::getFirstMonomer()
{
    if ( this->_Monomers.size() != 0 ) {
        return this->_Monomers[0];
    }
    return _Nil<Monomer_O>();
}



void Oligomer_O::_fillMonomerAsString(Monomer_sp mon, stringstream& seq)
{_G();
    seq <<"( ";
    seq << "monomer '" << mon->getId()->__repr__();
    if ( mon->getGroupName() == mon->getName() )
    {
	seq << " '" << mon->getGroupName()->symbolName();
    } else {
	seq << " (GroupPart :group '" << mon->getGroupName()->symbolName() <<" :part '" <<mon->getName()->symbolName() << ")";
    }
    adapt::SymbolSet_sp aliases = mon->getMonomerAliases();
    if ( aliases->size() > 0 )
    {
	seq << " :aliases '(" << aliases->asString() << " ) ";
    }
    seq <<" )";
}


void Oligomer_O::_fillSequenceAsStringForChildren(Monomer_sp rootMonomer, stringstream& seq)
{_G();
    gctools::Vec0<Coupling_sp>	outCouplings;
    gctools::Vec0<Coupling_sp>::iterator	oci;
    outCouplings = rootMonomer->getOutCouplings();
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
        Monomer_sp mon1 = (*oci)->getMonomer1();
        Monomer_sp mon2 = (*oci)->getMonomer2();
        seq << "( link '"<<mon1->getId()->symbolName();
	seq <<" '"<<(*oci)->getName()->symbolName() << " ";
	this->_fillMonomerAsString(mon2,seq);
	seq <<" )"<<std::endl;
    }
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
	if ((*oci).isA<RingCoupling_O>()) continue;
        Monomer_sp mon2 = (*oci)->getMonomer2();
        this->_fillSequenceAsStringForChildren(mon2,seq);
    }
}

CL_LISPIFY_NAME("sequenceAsString");
CL_DEFMETHOD string	Oligomer_O::sequenceAsString()
{_G();

    Monomer_sp mon2 = this->rootMonomer();
    stringstream seq;
    seq << "(Oligomer :name " << this->_Name->__repr__() << " :parts (list " << std::endl;
    this->_fillMonomerAsString(mon2,seq);
    seq << std::endl;
    this->_fillSequenceAsStringForChildren(mon2,seq);
    for (gctools::Vec0<Coupling_sp>::iterator ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ )
    {
	if ( (*ci).isA<RingCoupling_O>() )
	{
	    RingCoupling_sp rc = (*ci).as<RingCoupling_O>();
	    Monomer_sp mon1 = rc->getMonomer1();
	    Monomer_sp mon2 = rc->getMonomer2();
	    seq << "[ ringLink '";
	    seq << mon1->getId()->__repr__();
	    seq <<" '"<<rc->getPlug1() << " ";
	    seq << mon2->getId()->__repr__();
	    seq <<" '"<<rc->getPlug2() << " ";
	    seq <<" ]"<<std::endl;
	}
    }
    seq << "))"<<std::endl;
    return seq.str();
}



void Oligomer_O::_fillSequenceAsFileNameForChildren(Monomer_sp rootMonomer, stringstream& seq)
{_G();
    gctools::Vec0<Coupling_sp>	outCouplings;
    gctools::Vec0<Coupling_sp>::iterator	oci;
    outCouplings = rootMonomer->getOutCouplings();
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
//        Monomer_sp mon1 = (*oci)->getMonomer1();
        Monomer_sp mon2 = (*oci)->getMonomer2();
	seq <<"-" << mon2->getName()->symbolName();
    }
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
	if ( (*oci).isA<RingCoupling_O>() ) continue;
        Monomer_sp mon2 = (*oci)->getMonomer2();
        this->_fillSequenceAsFileNameForChildren(mon2,seq);
    }
}


adapt::ObjectSet_sp Oligomer_O::getAllAliases()
{_G();
    adapt::ObjectSet_sp aliases = adapt::ObjectSet_O::create();
    Oligomer_O::monomerIterator mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	adapt::ObjectSet_sp oneMonomerAliases = (*mi)->getAllAliases();
	aliases->addObjects(oneMonomerAliases);
    }
    return aliases;
}

adapt::SymbolSet_sp Oligomer_O::allMonomerAliases()
{_G();
    adapt::SymbolSet_sp aliases = adapt::SymbolSet_O::create();
    Oligomer_O::monomerIterator mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	adapt::SymbolSet_sp oneMonomerAliases = (*mi)->getMonomerAliases();
	aliases->insertSymbolSet(oneMonomerAliases);
    }
    return aliases;
}

CL_LISPIFY_NAME("sequenceAsFileName");
CL_DEFMETHOD string	Oligomer_O::sequenceAsFileName()
{_G();

    Monomer_sp mon2 = this->rootMonomer();
    stringstream seq;
    seq << mon2->getName()->symbolName();
    this->_fillSequenceAsFileNameForChildren(mon2,seq);
//    seq << "_" << this->getName()->symbolName();
    string name = core::stripCharacters(seq.str(),"()");
    return name.substr(0,255);
}

CL_LISPIFY_NAME("allSpecificMonomerContexts");
CL_DEFMETHOD SpecificContextSet_sp Oligomer_O::allSpecificMonomerContexts()
{_G();
    SpecificContextSet_sp allContexts = SpecificContextSet_O::create();
    Oligomer_O::monomerIterator mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	MonomerContext_sp monomerContext = (*mi)->getGeneralMonomerContext();
	SpecificContextSet_sp oneSpecificSet = monomerContext->getAllSpecificContexts();
	allContexts->merge(oneSpecificSet);
    }
    return allContexts;
}







#define ARGS_chem__oligomer_sequence "(olig)"
#define DECL_chem__oligomer_sequence ""
#define DOCS_chem__oligomer_sequence "oligomerSequence"
CL_DEFUN core::T_sp chem__oligomer_sequence(Oligomer_sp olig)
{_G();
    string seq = olig->sequenceAsString();
    return core::Str_O::create(seq);
}


void	Oligomer_O::_assembleFromParts(core::List_sp parts, CandoDatabase_sp bdb)
{_G();
    core::HashTableEq_sp monomerMap = core::HashTableEq_O::create_default();
    for ( auto p : parts ) {
	OligomerPart_Base_sp oligPart = p->car<OligomerPart_Base_O>();
	MultiMonomer_sp mon = oligPart->createMonomer(bdb);
	this->addMonomer(mon);
	monomerMap->setf_gethash(mon->getId(), mon);
	if ( oligPart.isA<OligomerPart_Link_O>() )
	{
          OligomerPart_Link_sp link = gc::As<OligomerPart_Link_sp>(oligPart);
	    core::Symbol_sp	mon1Id = link->_Monomer1Id;
	    core::Symbol_sp	mon2Id = link->_Monomer2->_MonomerId;
	    if ( !monomerMap->contains(mon1Id) ) SIMPLE_ERROR(BF("Unknown monomer id: %s")%mon1Id);
	    if ( !monomerMap->contains(mon2Id) ) SIMPLE_ERROR(BF("Unknown monomer id: %s")%mon2Id);
	    MultiMonomer_sp mon1 = monomerMap->gethash(mon1Id).as<MultiMonomer_O>();
	    MultiMonomer_sp mon2 = monomerMap->gethash(mon2Id).as<MultiMonomer_O>();
	    this->couple(mon1,link->_Coupling,mon2);
	}
    }
}

#if INIT_TO_FACTORIES

#define ARGS_Oligomer_O_make "(parts)"
#define DECL_Oligomer_O_make ""
#define DOCS_Oligomer_O_make "make Oligomer"
Oligomer_sp Oligomer_O::make(core::List_sp parts)
  {_G();
      GC_ALLOCATE(Oligomer_O, me );
    if ( parts.notnilp() )
      {
	CandoDatabase_sp bdb = getCandoDatabase();
	me->_assembleFromParts(parts,bdb);
      }
    return me;
  };

#else

core::T_sp Oligomer_O::__init__(core::Function_sp exec, core::List_sp args, core::Environment_sp env, core::Lisp_sp lisp)
{_OF();
    core::List_sp parts = translate::from_object<core::Cons_O>::convert(env->lookup(Pkg(),"parts"));
    if ( parts.notnilp() )
    {
	CandoDatabase_sp bdb = getCandoDatabase(lisp);
	this->_assembleFromParts(parts,bdb);
    }
    return _Nil<core::T_O>();
}

#endif




#define ARGS_chem__set_oligomer "(oligomerName parts)"
#define DECL_chem__set_oligomer ""
#define DOCS_chem__set_oligomer "setOligomer"
CL_DEFUN core::T_sp chem__set_oligomer(Oligomer_O::NameType::smart_ptr oligomerName, core::List_sp parts)
{_G();
    Oligomer_sp olig = Oligomer_O::create();
    olig->setName(oligomerName);
    CandoDatabase_sp bdb = getCandoDatabase();
    olig->_assembleFromParts(parts,bdb);
    core::Symbol_sp sym = oligomerName;
    sym->defparameter(olig);
//    lisp->globalEnvironment()->extend(sym,olig);
    return olig;
}



void Oligomer_O::exposeCando(core::Lisp_sp lisp)
{
    core::class_<Oligomer_O>()
//	.def_raw("core:__init__",&Oligomer_O::__init__,"(self &key parts)")
	.def("empty",&Oligomer_O::empty)
	.def("allSpecificMonomerContexts",&Oligomer_O::allSpecificMonomerContexts)
	.def("monomersAsCons",&Oligomer_O::monomersAsCons)
//	.def("addMultiMonomerAsLeaf",&Oligomer_O::addMultiMonomerAsLeaf)
	.def("removeLeafMonomer",&Oligomer_O::removeLeafMonomer)
	.def("hasLastMultiMonomerChanged",&Oligomer_O::hasLastMultiMonomerChanged)
	.def("getLastMultiMonomerChanged",&Oligomer_O::getLastMultiMonomerChanged)
	.def("sequenceAsFileName",&Oligomer_O::sequenceAsFileName)
	.def("setName",&Oligomer_O::setName,"","","",false)
	.def("getName",&Oligomer_O::getName,"","","",false)
//	.def("getCandoDatabase",&Oligomer_O::getCandoDatabase)
	.def("addMonomer", &Oligomer_O::addMonomer)
	.def("removeMonomer", &Oligomer_O::removeMonomer)
	.def("Oligomer-addCoupling", &Oligomer_O::addCoupling)
	.def("Oligomer-removeCoupling", &Oligomer_O::removeCoupling)
	.def("couple",&Oligomer_O::couple)
	.def("getMonomerWithId",&Oligomer_O::getMonomerWithId)
	.def("hasMonomerWithId",&Oligomer_O::hasMonomerWithId)
	.def("checkForErrors",&Oligomer_O::checkForErrors)
	.def("ringCouple",&Oligomer_O::ringCouple)
	.def("rootMonomer", &Oligomer_O::rootMonomer)
	.def("getMolecule", &Oligomer_O::getMolecule)
	.def("getHasError",&Oligomer_O::getHasError)
	.def("getErrorMessage",&Oligomer_O::getErrorMessage)
//	.def("getMonomers", &Oligomer_O::python_getMonomers)
//	.def("getMonomersInSequence",
//		&Oligomer_O::python_getMonomersInSequence)
//	.def("getCouplings", &Oligomer_O::python_getCouplings)
//	.def("getMonomerWithId",&Oligomer_O::getMonomerWithId)
	.def("getFirstMonomer", &Oligomer_O::getFirstMonomer)
	.def("getVerbose", &Oligomer_O::getVerbose)
	.def("throwIfBadConnections", &Oligomer_O::throwIfBadConnections)
	.def("sequenceAsString", &Oligomer_O::sequenceAsString)
	.def("setVerbose", &Oligomer_O::setVerbose)
	.def("updateMultiMonomers", &Oligomer_O::updateMultiMonomers)
	.def("numberOfSequences", &Oligomer_O::numberOfSequences)
//	.def("getCouplingWithSequenceNumber", &Oligomer_O::getCouplingWithSequenceNumber)
//	.def("getConformation",&Oligomer_O::getConformation)
//	.def("asXml",&Oligomer_O::asXml)
//	.def("assignResidueStorageIdsToMonomers",&Oligomer_O::assignResidueStorageIdsToMonomers)
//	.def("parseFromXml",&Oligomer_O::parseFromXml)
//	.def("mutateStereoisomer",&Oligomer_O::mutateStereoisomer)
//	.def("buildMolecule",&Oligomer_O::buildMolecule)
//	.def("unbuildMolecule",&Oligomer_O::unbuildMolecule)
//	.def("getAreAlchemistsReady",&Oligomer_O::getAreAlchemistsReady)
//	.def("clearAllFragmentAlchemistsCoordinatesAreFixed",
//		&Oligomer_O::clearAllFragmentAlchemistsCoordinatesAreFixed)
    ;
//    Chem_temp_Defun(set_oligomer);
//    Chem_temp_Defun(oligomer_sequence);
}




void Oligomer_O::exposePython(core::Lisp_sp lisp)
{_G();
#ifdef	USEBOOSTPYTHON //[
    PYTHON_CLASS(ChemPkg,Oligomer,"","",_lisp)
//	.def("defineFromMolecule", &Oligomer_O::defineFromMolecule )
	.add_property("iterate_Monomers",
			boost::python::range(&Oligomer_O::begin_Monomers,
				&Oligomer_O::end_Monomers))
	.add_property("iterate_Couplings",
			boost::python::range(&Oligomer_O::begin_Couplings,
				&Oligomer_O::end_Couplings))
	.def("empty",&Oligomer_O::empty)
	.def("sequenceAsFileName",&Oligomer_O::sequenceAsFileName)
	.def("monomersAsCons",&Oligomer_O::monomersAsCons)
//	.def("addMultiMonomerAsLeaf",&Oligomer_O::addMultiMonomerAsLeaf)
	.def("removeLeafMonomer",&Oligomer_O::removeLeafMonomer)
	.def("hasLastMultiMonomerChanged",&Oligomer_O::hasLastMultiMonomerChanged)
	.def("getLastMultiMonomerChanged",&Oligomer_O::getLastMultiMonomerChanged)
    //	.def("setName",&Oligomer_O::setName,"","","",false)
    //	.def("getName",&Oligomer_O::getName,"","","",false)
//	.def("getCandoDatabase",&Oligomer_O::getCandoDatabase)
	.def("addMonomer", &Oligomer_O::addMonomer)
	.def("removeMonomer", &Oligomer_O::removeMonomer)
	.def("addCoupling", &Oligomer_O::addCoupling)
	.def("getMonomerWithId",&Oligomer_O::getMonomerWithId)
	.def("hasMonomerWithId",&Oligomer_O::hasMonomerWithId)
	.def("removeCoupling", &Oligomer_O::removeCoupling)
	.def("couple",&Oligomer_O::couple)
	.def("checkForErrors",&Oligomer_O::checkForErrors)
	.def("ringCouple",&Oligomer_O::ringCouple)
	.def("rootMonomer", &Oligomer_O::rootMonomer)
	.def("getMolecule", &Oligomer_O::getMolecule)
	.def("getHasError",&Oligomer_O::getHasError)
	.def("getErrorMessage",&Oligomer_O::getErrorMessage)
//	.def("getMonomers", &Oligomer_O::python_getMonomers)
//	.def("getMonomersInSequence",
//		&Oligomer_O::python_getMonomersInSequence)
//	.def("getCouplings", &Oligomer_O::python_getCouplings)
//	.def("getMonomerWithId",&Oligomer_O::getMonomerWithId)
	.def("getFirstMonomer", &Oligomer_O::getFirstMonomer)
	.def("getVerbose", &Oligomer_O::getVerbose)
	.def("throwIfBadConnections", &Oligomer_O::throwIfBadConnections)
	.def("sequenceAsString", &Oligomer_O::sequenceAsString)
	.def("setVerbose", &Oligomer_O::setVerbose)
	.def("updateMultiMonomers", &Oligomer_O::updateMultiMonomers)
	.def("numberOfSequences", &Oligomer_O::numberOfSequences)
//	.def("getCouplingWithSequenceNumber", &Oligomer_O::getCouplingWithSequenceNumber)
//	.def("getConformation",&Oligomer_O::getConformation)
//	.def("asXml",&Oligomer_O::asXml)
//	.def("assignResidueStorageIdsToMonomers",&Oligomer_O::assignResidueStorageIdsToMonomers)
//	.def("parseFromXml",&Oligomer_O::parseFromXml)
//	.def("mutateStereoisomer",&Oligomer_O::mutateStereoisomer)
//	.def("buildMolecule",&Oligomer_O::buildMolecule)
//	.def("unbuildMolecule",&Oligomer_O::unbuildMolecule)
//	.def("getAreAlchemistsReady",&Oligomer_O::getAreAlchemistsReady)
//	.def("clearAllFragmentAlchemistsCoordinatesAreFixed",
//		&Oligomer_O::clearAllFragmentAlchemistsCoordinatesAreFixed)
	;

//    boost::python::def("create_Oligomer",&create_Oligomer );
#endif
}


EXPOSE_CLASS(chem,Oligomer_O);

}; // namespace chem



