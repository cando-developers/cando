/*
    File: oligomer.cc
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
//#define	DEBUG_LEVEL_FULL


//
// (C) 2004 Christian E. Schafmeister
//

#include <cando/chem/oligomer.h>
#include <clasp/core/array.h>
#include <clasp/core/hashTableEq.h>
#include <clasp/core/bformat.h>
#include <cando/chem/molecule.h>
#include <clasp/core/numerics.h>
#include <cando/chem/candoScript.h>
#include <cando/adapt/adapters.h>
#include <cando/chem/topology.h>
#include <clasp/core/evaluator.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/monomer.h>
#include <cando/chem/coupling.h>
#include <cando/adapt/adapters.h>
#include <clasp/core/lispStream.h>
#include <cando/chem/plug.h>
#include <cando/chem/constitution.h>
#include <cando/chem/specificContext.h>
#include <clasp/core/wrappers.h>

/*! \page oligomerBuilding  How to build oligomers
 * Building oligomers requires the following steps.
 * 	-# Define the oligomer
 * 		- Create the _Monomers, add them to the oligomer and couple them together.
 * 	-# Construct the alchemists Oliggomer_O::createAlchemists (you may need to dismissAlchemists first)
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


void Oligomer_O::fields(core::Record_sp node)
{
  node->field_if_not_empty( INTERN_(kw,monomers), this->_Monomers);
  node->field_if_not_empty( INTERN_(kw,couplings), this->_Couplings);
}


Oligomer_O::Oligomer_O(const Oligomer_O& original)
{
  core::HashTableEq_sp newMonomersFromOld = core::HashTableEq_O::create_default();
  gctools::Vec0<Monomer_sp> copiedMonomers;
  this->_Monomers.resize(original._Monomers.size());
  for ( size_t i=0, iEnd(this->_Monomers.size()); i<iEnd; ++i ) {
    Monomer_sp oldMonomer = original._Monomers[i];
    Monomer_O& ref = *oldMonomer;
    GC_COPY(Monomer_O,monomerCopy,ref);
    this->_Monomers[i] = monomerCopy;
    newMonomersFromOld->setf_gethash(oldMonomer,monomerCopy);
  }
  this->_Couplings.resize(original._Couplings.size());
  for ( size_t ci=0, ciEnd(this->_Couplings.size()); ci<ciEnd; ++ci ) {
    Coupling_sp oldCoupling = original._Couplings[ci];
    Monomer_sp oldMonomer1 = oldCoupling->getMonomer1();
    Monomer_sp oldMonomer2 = oldCoupling->getMonomer2();
    Monomer_sp newMonomer1 = gc::As<Monomer_sp>(newMonomersFromOld->gethash(oldMonomer1));
    Monomer_sp newMonomer2 = gc::As<Monomer_sp>(newMonomersFromOld->gethash(oldMonomer2));
    if ( gc::IsA<DirectionalCoupling_sp>(oldCoupling) ) {
      DirectionalCoupling_sp dc = gc::As_unsafe<DirectionalCoupling_sp>(oldCoupling);
      DirectionalCoupling_sp newCoupling = DirectionalCoupling_O::make(newMonomer1,dc->getSourcePlugName(),dc->getTargetPlugName(),newMonomer2);
      this->_Couplings[ci] = newCoupling;
    } else if (gc::IsA<RingCoupling_sp>(oldCoupling) ) {
      RingCoupling_sp rc = gc::As_unsafe<RingCoupling_sp>(oldCoupling);
      RingCoupling_sp newCoupling = RingCoupling_O::make(newMonomer1,rc->getPlug1(),newMonomer2,rc->getPlug2());
      this->_Couplings[ci] = newCoupling;
    }
  }
  this->reinitializeData();
}

void Oligomer_O::initialize()
{_OF();
    this->Base::initialize();
    this->empty();
}


/*! Update the Monomers so that they refer to the proper couplings */
void Oligomer_O::reinitializeData()
{
  for ( size_t i(0),iEnd(this->_Monomers.size()); i<iEnd; ++i ) {
    this->_Monomers[i]->eraseCouplings();
    this->_Monomers[i]->setSequenceNumber(i);
  }
  for ( size_t ci(0),ciEnd(this->_Couplings.size()); ci<ciEnd; ++ci ) {
    Coupling_sp coupling = this->_Couplings[ci];
    if ( gc::IsA<DirectionalCoupling_sp>(coupling) ) {
      DirectionalCoupling_sp dc = gc::As_unsafe<DirectionalCoupling_sp>(coupling);
      core::Symbol_sp inPlugName = DirectionalCoupling_O::inPlugName(dc->getName());
      core::Symbol_sp outPlugName = DirectionalCoupling_O::outPlugName(dc->getName());
      coupling->getMonomer1()->addCoupling(outPlugName,dc);
      coupling->getMonomer2()->addCoupling(inPlugName,dc);
    } else if (gc::IsA<RingCoupling_sp>(coupling) ) {
      RingCoupling_sp rc = gc::As_unsafe<RingCoupling_sp>(coupling);
      core::Symbol_sp plug1Name = rc->getPlug1();
      core::Symbol_sp plug2Name = rc->getPlug2();
      coupling->getMonomer1()->addCoupling(plug1Name,rc);
      coupling->getMonomer2()->addCoupling(plug2Name,rc);
    }
  }
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

CL_LISPIFY_NAME("monomersAsList");
CL_DEFMETHOD core::List_sp Oligomer_O::monomersAsList()
{
  ql::list result;
  gctools::Vec0<Monomer_sp>::iterator	mi;
  for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
  {
    result << *mi;
  }
  return result.cons();
}

CL_LISPIFY_NAME("couplings-as-list");
CL_DEFMETHOD core::List_sp Oligomer_O::couplingsAsList()
{
  ql::list result;
  core::List_sp cons = _Nil<core::T_O>();
  gctools::Vec0<Coupling_sp>::iterator	mi;
  for ( mi=this->_Couplings.begin(); mi!=this->_Couplings.end(); mi++ )
  {
    result << *mi;
  }
  return result.cons();
}



void	Oligomer_O::expandMonomerListToNeighbors(gctools::SmallOrderedSet<Monomer_sp>& monomers)
{
  gctools::SmallOrderedSet<Monomer_sp>				expanded;
  gctools::SmallOrderedSet<Monomer_sp>::iterator		mi;
  Monomer_sp				neighbor;
  Monomer_O::Couplings::iterator	        ci;
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
{
    gctools::SmallOrderedSet<Monomer_sp>			expanded;
    ASSERTNOTNULL(cdb);
    expanded = monomers;
    this->expandMonomerListToNeighbors(expanded);
    this->checkMonomersAndNotNeighborsForErrors(cdb,expanded);
}


    void	Oligomer_O::checkMonomersAndNotNeighborsForErrors(CandoDatabase_sp cdb, gctools::SmallOrderedSet<Monomer_sp>& monomers)
{
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
{
Monomer_sp	monomer;
uint		idx;
    monomer = s->sharedThis<Monomer_O>();
    LOG(BF("s@%p    monomer@%X") % s.get() % monomer.get() );
    LOG(BF("Adding monomer@%X with name (%s)") % monomer.get() % s->description().c_str() );
    idx = this->_Monomers.size();
    monomer->setSequenceNumber(idx);
    this->_Monomers.push_back(monomer);
//    monomer->connectListener(this->sharedThis<Oligomer_O>(),Monomer_contentsChanged);
    return idx;
}

CL_LISPIFY_NAME("hasMonomerWithId");
CL_DEFMETHOD bool Oligomer_O::hasMonomerWithId(core::Symbol_sp id)
{
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
    {
	if ( (*mi)->getId() == id ) return true;
    }
    return false;
}


CL_LISPIFY_NAME("getMonomerWithId");
CL_DEFMETHOD Monomer_sp Oligomer_O::getMonomerWithId(core::Symbol_sp id)
{
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
{
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
{
  IMPLEMENT_ME();
#if 0
    gctools::Vec0<Monomer_sp>::iterator	mi;
    gctools::Vec0<Coupling_sp>::iterator	ci;
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
#endif
}



CL_LISPIFY_NAME("updateMultiMonomers");
CL_DEFMETHOD void	Oligomer_O::updateMultiMonomers()
{
Oligomer_O::monomerIterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        (*mi)->updateGroupName();
    }
}



/*! Add a Monomer as a leaf to the indicated monomer
 */
Monomer_sp	Oligomer_O::addMultiMonomerAsLeaf( Monomer_sp monFrom, core::Symbol_sp couplingName,
						   core::Symbol_sp representedEntityNameSetName, geom::Vector2 pos )
{
DirectionalCoupling_sp 	coup;
Monomer_sp	monTo;

    coup = DirectionalCoupling_O::create();
    coup->setPlugNames(DirectionalCoupling_O::outPlugName(couplingName),DirectionalCoupling_O::inPlugName(couplingName));
    monTo = Monomer_O::create();
//    monTo->setGroupName(representedEntityNameSetName);
    coup->setTargetMonomer(monTo);
    monTo->setInCoupling(coup);
    	// Now attach it
    monFrom->addOutCoupling(coup);
    coup->setSourceMonomer(monFrom);
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
{
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
{
    this->_Couplings.push_back(c);
}



/*! Return the monomer at the root of the tree.
Oligomer._State can be anything. */
CL_LISPIFY_NAME("rootMonomer");
CL_DEFMETHOD Monomer_sp	Oligomer_O::rootMonomer() const
{
  Monomer_sp	sub;
  Coupling_sp	coupling;
  coupling = _Unbound<Coupling_O>();
  if (this->_Monomers.size()==0) {
    SIMPLE_ERROR(BF("There are no monomers, root-monomer can not be determined"));
  }
  sub = this->_Monomers[0];
  do {
    LOG(BF("Looking at monomer: %s") % sub->description().c_str()  );
    if ( sub->hasInCoupling()) {
      LOG(BF("    It has an in-coupling") );
      coupling = sub->getInCoupling();
      if ( !gc::IsA<DirectionalCoupling_sp>(coupling)) {
        LOG(BF("     But its not a directional") );
        break;
      }
      sub = gc::As_unsafe<DirectionalCoupling_sp>(coupling)->getSourceMonomer();
    } else break;
  } while (!coupling.unboundp());
  return sub;
}


CL_LISPIFY_NAME("throwIfBadConnections");
CL_DEFMETHOD void	Oligomer_O::throwIfBadConnections()
{
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




CL_LAMBDA("(oligomer !) source-monomer coupling-or-source-plug-name target-monomer &optional (target-plug-name nil)")
CL_LISPIFY_NAME("couple");
CL_DEFMETHOD DirectionalCoupling_sp	Oligomer_O::couple( Monomer_sp sourceMon, core::T_sp couplingOrSourcePlugName, Monomer_sp targetMon, core::Symbol_sp targetPlugName  )
{
  DirectionalCoupling_sp	coupling;
  bool		foundIn, foundOut;
  gctools::Vec0<Monomer_sp>::iterator	mi;
  core::Symbol_sp sourcePlugName;
  if (couplingOrSourcePlugName.consp()) {
    core::Cons_sp pair = gc::As_unsafe<core::Cons_sp>(couplingOrSourcePlugName);
    sourcePlugName = gc::As<core::Symbol_sp>(CONS_CAR(pair));
    if (!targetPlugName.nilp()) {
      SIMPLE_ERROR(BF("coupling-or-source-plug-name was a pair and so target-plug-name must be NIL - but it is %s") % _rep_(targetPlugName));
    }
    targetPlugName = gc::As<core::Symbol_sp>(CONS_CDR(pair));
  } else if (gc::IsA<core::Symbol_sp>(couplingOrSourcePlugName)) {
    core::Symbol_sp sym = gc::As_unsafe<core::Symbol_sp>(couplingOrSourcePlugName);
    if (!DirectionalCoupling_O::isOutPlugName(couplingOrSourcePlugName)&&!targetPlugName.nilp()) {
      SIMPLE_ERROR(BF("coupling-or-source-plug-name was a coupling name and so target-plug-name must be NIL - but it is %s") % _rep_(targetPlugName));
    }
    if (DirectionalCoupling_O::isCouplingName(sym)) {
      sourcePlugName = DirectionalCoupling_O::outPlugName(sym);
      targetPlugName = DirectionalCoupling_O::inPlugName(sym);
    } else {
      if (!DirectionalCoupling_O::isOutPlugName(sym)) {
        SIMPLE_ERROR(BF("coupling-or-source-plug-name must be a source-plug-name - instead it is %s") % _rep_(sym));
      }
      sourcePlugName = sym;
      if (!DirectionalCoupling_O::isInPlugName(targetPlugName)) {
        SIMPLE_ERROR(BF("target-plug-name must be a target-plug-name - instead it is %s") % _rep_(targetPlugName));
      }
    }
  } else {
    SIMPLE_ERROR(BF("coupling-or-source-plug-name must be a cons or a symbol - it is %s") % _rep_(couplingOrSourcePlugName));
  }
  LOG(BF("Creating a coupling from %s to %s") % _rep_(sourceMon) % _rep_(targetMon));
  foundIn = false;
  foundOut = false;
  for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) 
  {
    if ( *mi == sourceMon ) {
      foundIn = true;
    }
    if ( *mi == targetMon ) {
      foundOut = true;
    }
  }
  if ( !foundIn ) {
    SIMPLE_ERROR(BF("Could not find first monomer in oligomer"));
  }
  if ( !foundOut ) {
    SIMPLE_ERROR(BF("Could not find second monomer in oligomer"));
  }
#if 0
  // Turn this off for now until we know how the new coupling names work
  core::Symbol_sp targetMonPlugName = DirectionalCoupling_O::inPlugName(name);
  if ( targetMon->hasCouplingWithPlugName(targetMonPlugName) ) {
    SIMPLE_ERROR(BF("The second monomer already has an in coupling"));
  }
  core::Symbol_sp sourceMonPlugName = DirectionalCoupling_O::outPlugName(name);
  if ( sourceMon->hasCouplingWithPlugName(sourceMonPlugName) ) {
    SIMPLE_ERROR(BF("The first monomer[%s] already has an out coupling with the name: %s") % sourceMon->currentStereoisomerName() % name );
  }
#endif
  LOG(BF("in/out plug names: %s/%s") % _rep_(sourcePlugName) % _rep_(targetPlugName));
  coupling = DirectionalCoupling_O::make(sourceMon,sourcePlugName,targetPlugName,targetMon);
  sourceMon->addCoupling(sourcePlugName,coupling);
  targetMon->addCoupling(targetPlugName,coupling);
  this->addCoupling(coupling);
  LOG(BF("after coupling coupling: %s") % coupling->description()  );
  return coupling;
}


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
    RingClosingPlug_sp rcPlug1 = gc::As<OutPlug_sp>(mon1->getMissingRingClosingPlug(mon2));
    if ( rcPlug1.nilp() )
    {
	SIMPLE_ERROR(BF("The %s can not make a ring closing coupling")
					 % mon1->description() );
    }
    RingClosingPlug_sp rcPlug2 = gc::As<OutPlug_sp>(mon2->getMissingRingClosingPlug(mon1));
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

CL_DEFMETHOD RingCoupling_sp	Oligomer_O::ringCoupleWithPlugNames( Monomer_sp mon1, core::Symbol_sp plug1name, Monomer_sp mon2, core::Symbol_sp plug2name )
{
  RingCoupling_sp	coupling;
  bool		found1, found2;
  gctools::Vec0<Monomer_sp>::iterator	mi;
  found1 = false;
  found2 = false;
  for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) {
    if ( *mi == mon1 ) {
      found1 = true;
    }
    if ( *mi == mon2 ) {
      found2 = true;
    }
  }
  if ( !found1 ) 
    SIMPLE_ERROR(BF("Could not find monomer1 in oligomer"));
  if ( !found2 ) 
    SIMPLE_ERROR(BF("Could not find monomer2 in oligomer"));
  coupling = RingCoupling_O::make(mon1->sharedThis<Monomer_O>(),plug1name,mon2->sharedThis<Monomer_O>(),plug2name);
  mon1->addCoupling(plug1name,coupling);
  mon2->addCoupling(plug2name,coupling);
  this->addCoupling(coupling);
  return coupling;
}




CL_LISPIFY_NAME("getMolecule");
CL_DEFMETHOD Molecule_sp Oligomer_O::getMolecule()
{
Topology_sp		topology;
int			residueNetCharge;
    Residue_sp			res;
    Molecule_sp 			mol;
    mol = Molecule_O::create();
    core::write_bf_stream(BF("%s:%d Creating residues\n") % __FILE__ % __LINE__ );
    core::HashTableEq_sp monomersToResidues = core::HashTableEq_O::create_default();
    gctools::Vec0<Monomer_sp>::iterator	mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ ) {
      Topology_sp topology = (*mi)->currentTopology();
      Constitution_sp constitution = topology->_Constitution;
      core::Symbol_sp stereoisomerName = (*mi)->currentStereoisomerName(); // Should be stereoisomerName()
      if (chem__verbose(0)) {
        core::write_bf_stream(BF("Building residue for monomer name: %s") % _rep_(stereoisomerName));
        core::clasp_finish_output_t();
      }
      res = topology->buildResidueForMonomerName(stereoisomerName);
      LOG(BF("made res consistent with stereoisomer named res: %s") % _rep_(res));
      mol->addMatter(res);
      LOG(BF("Added matter"));
      monomersToResidues->hash_table_setf_gethash((*mi),res);
    }
//    core::write_bf_stream(BF("%s:%d Connecting residues\n") % __FILE__ % __LINE__ );
    gctools::Vec0<Coupling_sp>::iterator	ci;
    Monomer_sp	mon1, mon2;
    Residue_sp	res1, res2;
    for ( ci=this->_Couplings.begin(); ci!=this->_Couplings.end(); ci++ ) 
    {
	mon1 = (*ci)->getMonomer1();
	ASSERT(mon1.notnilp());
	mon2 = (*ci)->getMonomer2();
	ASSERT(mon2.notnilp());
        res1 = gc::As_unsafe<Residue_sp>(monomersToResidues->gethash(mon1));
        res2 = gc::As_unsafe<Residue_sp>(monomersToResidues->gethash(mon2));
//        core::write_bf_stream(BF("%s:%d  Connecting residues %s - %s\n") % __FILE__ % __LINE__ % _rep_(res1) % _rep_(res2));
    	(*ci)->doCoupling(res1,res2);
    };
    return mol;
}



void Oligomer_O::_gatherMultiMonomers(gctools::Vec0<Monomer_sp>& multiMonomers)
{
    gctools::Vec0<Monomer_sp>::iterator	mi;
    multiMonomers.clear();
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	if ( (*mi)->numberOfStereoisomers()>1 )
	{
	    multiMonomers.push_back(*mi);
	}
    }
    if ( multiMonomers.size() == 0 )
    {
	LOG(BF("There are no Monomers with more than 1 possibles") );
    }
}


#if 0
/*!
 * Perturb one MultiMonomer in the Oligomer.
 * If there are none then do nothing.
 */
void	Oligomer_O::perturbMonomerSequence()
{
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
#endif

CL_DEFMETHOD void	Oligomer_O::gotoSequence(core::Integer_sp index)
{
  vector<int>	digits, bases;
  gctools::Vec0<Monomer_sp>::iterator	mi;
  vector<int>::iterator	di;
  Bignum idx = core::clasp_to_mpz(index);
  for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
  {
    bases.push_back((*mi)->numberOfStereoisomers());
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
CL_DEFMETHOD void	Oligomer_O::firstSequence()
{
    gctools::Vec0<Monomer_sp>	multiMonomers;
    this->gotoSequence(core::Integer_O::create(0));
}


/*! Convert the current sequence into an index
 */
CL_DEFMETHOD Bignum Oligomer_O::currentSequenceIndex()
{
vector<int>	digits, bases;
gctools::Vec0<Monomer_sp>::iterator	mi;
Bignum			index = 0;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        bases.push_back((*mi)->numberOfStereoisomers());
	digits.push_back((*mi)->getMonomerIndex());
    }
    index = core::mixedBaseDigitsToBignum(bases,digits);
    return index;
}


CL_LISPIFY_NAME("monomerWithSequenceNumber");
CL_DEFMETHOD Monomer_sp Oligomer_O::monomerWithSequenceNumber(size_t sequence) const
{
  if (sequence <= this->_Monomers.size()) {
    return this->_Monomers[sequence];
  }
  SIMPLE_ERROR(BF("Monomer sequence number %lu must be less than %lu")% sequence % this->_Monomers.size());
}

/*!
 * Increment the monomer sequence.
 * Return -1 if there are too many sequences.
 */
CL_LISPIFY_NAME("numberOfSequences");
CL_DEFMETHOD Bignum Oligomer_O::numberOfSequences()
{
    gctools::Vec0<Monomer_sp>::iterator	mi;
vector<int>			bases;
Bignum				numSeq;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        bases.push_back((*mi)->numberOfStereoisomers());
    }
    numSeq = core::numberOfIndicesForMixedBase(bases);
    return numSeq;
}

/*!
 * Increment the monomer sequence.
 * Return true if it incremented, return false if it overflowed
 */
CL_DEFMETHOD bool Oligomer_O::incrementSequence()
{
    gctools::Vec0<Monomer_sp>::iterator	mi;
long unsigned int numSeq;
    numSeq = 1;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
        if ( (*mi)->incrementMonomerIndex() ) return true;
        numSeq *= (*mi)->numberOfStereoisomers();
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



void Oligomer_O::_fillMonomerAsString(Monomer_sp mon, stringstream& seq) const
{
    seq <<"( ";
    seq << "monomer '" << _rep_(mon->getId());
#if 1
    seq << " ( :part '" << _rep_(mon->currentStereoisomerName()) <<  ")";
#else
    if ( mon->getGroupName() == mon->currentStereoisomerName() )
    {
      seq << " '" << mon->getGroupName()->symbolNameAsString();
    } else {
      seq << " (GroupPart :group '" << mon->getGroupName()->symbolNameAsString() <<" :part '" <<mon->currentStereoisomerName()->symbolNameAsString() << " :isomer " << core::_rep_(mon->getIsomer()) <<  ")";
    }
#endif
#if 0
    adapt::SymbolSet_sp aliases = mon->getMonomerAliases();
    if ( aliases->size() > 0 )
    {
	seq << " :aliases '(" << aliases->asString() << " ) ";
    }
#endif
    seq <<" )";
}


void Oligomer_O::_fillSequenceAsStringForChildren(Monomer_sp rootMonomer, stringstream& seq) const
{
    gctools::Vec0<Coupling_sp>	outCouplings;
    gctools::Vec0<Coupling_sp>::iterator	oci;
    outCouplings = rootMonomer->getOutCouplings();
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
        Monomer_sp mon1 = (*oci)->getMonomer1();
        Monomer_sp mon2 = (*oci)->getMonomer2();
//        printf("%s:%d:%s About to build sequence \n", __FILE__, __LINE__, __FUNCTION__ );
//        printf("%s:%d:%s About to build sequence mon1 -> %s\n", __FILE__, __LINE__, __FUNCTION__, core::_rep_(mon1).c_str() );
        seq << "( link '"<<mon1->getId()->symbolNameAsString();
	seq <<" '"<<(*oci)->getName()->symbolNameAsString() << " ";
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

string Oligomer_O::__repr__() const {
  stringstream ss;
  ss << "#<" << this->className() << " ";
  if (this->_Monomers.size()>0) {
    ss << this->sequenceAsFileName();
  }
  ss << ">";
  return ss.str();
}

CL_LISPIFY_NAME("sequenceAsString");
CL_DEFMETHOD string	Oligomer_O::sequenceAsString() const
{
  Monomer_sp mon2 = this->rootMonomer();
  stringstream seq;
  seq << "(Oligomer  :parts (list " << std::endl;
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
      seq << _rep_(mon1->getId());
      seq <<" '"<<rc->getPlug1() << " ";
      seq << _rep_(mon2->getId());
      seq <<" '"<<rc->getPlug2() << " ";
      seq <<" ]"<<std::endl;
    }
  }
  seq << "))"<<std::endl;
  return seq.str();
}



void Oligomer_O::_fillSequenceAsFileNameForChildren(Monomer_sp rootMonomer, stringstream& seq) const
{
    gctools::Vec0<Coupling_sp>	outCouplings;
    gctools::Vec0<Coupling_sp>::iterator	oci;
    outCouplings = rootMonomer->getOutCouplings();
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
//        Monomer_sp mon1 = (*oci)->getMonomer1();
        Monomer_sp mon2 = (*oci)->getMonomer2();
	seq <<"-" << mon2->currentStereoisomerName()->symbolNameAsString();
    }
    for ( oci=outCouplings.begin(); oci!=outCouplings.end(); oci++ )
    {
	if ( (*oci).isA<RingCoupling_O>() ) continue;
        Monomer_sp mon2 = (*oci)->getMonomer2();
        this->_fillSequenceAsFileNameForChildren(mon2,seq);
    }
}

#if 0
adapt::ObjectSet_sp Oligomer_O::getAllAliases()
{
    adapt::ObjectSet_sp aliases = adapt::ObjectSet_O::create();
    Oligomer_O::monomerIterator mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	adapt::ObjectSet_sp oneMonomerAliases = (*mi)->getAllAliases();
	aliases->addObjects(oneMonomerAliases);
    }
    return aliases;
}
#endif

#if 0
adapt::SymbolSet_sp Oligomer_O::allMonomerAliases()
{
    adapt::SymbolSet_sp aliases = adapt::SymbolSet_O::create();
    Oligomer_O::monomerIterator mi;
    for ( mi=this->_Monomers.begin(); mi!=this->_Monomers.end(); mi++ )
    {
	adapt::SymbolSet_sp oneMonomerAliases = (*mi)->getMonomerAliases();
	aliases->insertSymbolSet(oneMonomerAliases);
    }
    return aliases;
}
#endif


CL_LISPIFY_NAME("sequenceAsFileName");
CL_DEFMETHOD string	Oligomer_O::sequenceAsFileName() const
{

    Monomer_sp mon2 = this->rootMonomer();
    stringstream seq;
    seq << mon2->currentStereoisomerName()->symbolNameAsString();
    this->_fillSequenceAsFileNameForChildren(mon2,seq);
//    seq << "_" << this->getName()->symbolName();
    string name = core::stripCharacters(seq.str(),"()");
    return name.substr(0,255);
}

CL_LISPIFY_NAME("allSpecificMonomerContexts");
CL_DEFMETHOD SpecificContextSet_sp Oligomer_O::allSpecificMonomerContexts()
{
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
{
    string seq = olig->sequenceAsString();
    return core::Str_O::create(seq);
}


void	Oligomer_O::_assembleFromParts(core::List_sp parts, CandoDatabase_sp bdb)
{
  core::HashTableEq_sp monomerMap = core::HashTableEq_O::create_default();
  for ( auto p : parts ) {
    OligomerPart_Base_sp oligPart = p->car<OligomerPart_Base_O>();
    Monomer_sp mon = oligPart->createMonomer(bdb);
    this->addMonomer(mon);
    monomerMap->setf_gethash(mon->getId(), mon);
    if ( oligPart.isA<OligomerPart_Link_O>() )
    {
      OligomerPart_Link_sp link = gc::As<OligomerPart_Link_sp>(oligPart);
      core::Symbol_sp	mon1Id = link->_Monomer1Id;
      core::Symbol_sp	mon2Id = link->_Monomer2->_MonomerId;
      if ( !monomerMap->contains(mon1Id) ) SIMPLE_ERROR(BF("Unknown monomer id: %s")%mon1Id);
      if ( !monomerMap->contains(mon2Id) ) SIMPLE_ERROR(BF("Unknown monomer id: %s")%mon2Id);
      Monomer_sp mon1 = monomerMap->gethash(mon1Id).as<Monomer_O>();
      Monomer_sp mon2 = monomerMap->gethash(mon2Id).as<Monomer_O>();
      this->couple(mon1,link->_Coupling,mon2);
    }
  }
}

CL_DEFMETHOD core::T_sp Oligomer_O::copyOligomer() const
{
  GC_COPY(Oligomer_O,newOligomer,*this);
  return newOligomer;
}

#define ARGS_Oligomer_O_make "(parts)"
#define DECL_Oligomer_O_make ""
#define DOCS_Oligomer_O_make "make Oligomer"
Oligomer_sp Oligomer_O::make(core::List_sp parts)
  {
      GC_ALLOCATE(Oligomer_O, me );
    if ( parts.notnilp() )
      {
	CandoDatabase_sp bdb = getCandoDatabase();
	me->_assembleFromParts(parts,bdb);
      }
    return me;
  };


}; // namespace chem
