#define	DEBUG_LEVEL_FULL
//
// (C) 2004 Christian E. Schafmeister
//




#include <string.h>
#include <iostream>
#include <vector>
//#include	<strings>
#include <stdio.h>
#include <clasp/core/common.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/matter.h>
#include <cando/geom/omatrix.h>
//#include "core/serialize.h"
#include <clasp/core/str.h>
#include <cando/chem/loop.h>
#include <clasp/core/cons.h>
#include <cando/geom/ovector3.h>
#include <cando/chem/restraint.h>
#include <clasp/core/environment.h>
#include <cando/geom/boundingBox.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/bond.h>
#include <cando/chem/angle.h>
#include <cando/chem/properTorsion.h>
#include <cando/chem/improperTorsion.h>
#include <cando/chem/atomIdMap.h>
#include <cando/chem/alias.h>
#include <cando/chem/symbolTable.h>
#include <clasp/core/translators.h>
#include <clasp/core/wrappers.h>

namespace chem
{





//
// Constructor
//

    void	Matter_O::initialize()
    {
	this->Base::initialize();
	this->_Id = 1; // UNDEFINED_ID;
	this->_NextContentId = 1;
	this->setContainedByNothing();
	this->_Restraints = _Nil<RestraintList_O>();
        this->_Properties = _Nil<core::T_O>();
    }

    string Matter_O::__repr__() const
    {
	stringstream ss;
	ss << this->className();
	return ss.str();
    }


//
// copy constructor
//
    Matter_O::Matter_O( const Matter_O& c ) : core::CxxObject_O(c)
    {
	this->_NextContentId = c._NextContentId;
	this->_Id = c._Id;
	this->name = c.name;
//	this->containerContainedBy = c.containerContainedBy;
	this->_Properties = c._Properties;
	this->_Restraints = _Nil<RestraintList_O>();
    }





//
// eraseContents
//
    void	Matter_O::eraseContents()
    {
	this->_contents.clear();
    }


//
// hasContentWithName
//
    bool	Matter_O::hasContentWithName(MatterName    sName )
    {_G();
	contentIterator		aCur;

	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    LOG(BF("Looking at(%s) for(%s)") % (*aCur)->getName().c_str() % sName.c_str()  );
	    if ( (*aCur)->getName() == sName ) {
		return true;
	    }
	}
	return false;
    }


    void Matter_O::accumulateRestraints(RestraintList_sp allRestraints) const
    {_OF();
	allRestraints->merge(this->_Restraints);
	for ( const_contentIterator ci=this->begin_contents(); ci!=this->end_contents(); ci++ )
	{
	    (*ci)->accumulateRestraints(allRestraints);
	}
    }


    RestraintList_sp Matter_O::allRestraints() const
    {
	GC_ALLOCATE(RestraintList_O, allRestraints );
	this->accumulateRestraints(allRestraints);
	return allRestraints;
    }



    void Matter_O::addRestraint(Restraint_sp restraint)
    {_OF();
	if ( this->_Restraints.nilp() )
	{
	    this->_Restraints = RestraintList_O::create();
	}
	this->_Restraints->addRestraint(restraint);
    }




    void Matter_O::resizeContents(int sz)
    {_OF();
	this->_contents.resize(sz,_Nil<Matter_O>());
    }


    void Matter_O::putMatter(int idx, Matter_sp matter)
    {_OF();
	ASSERTF(idx>=0 && idx<(int)this->_contents.size(),BF("Illegal putMatter index[%d] must be less than %d") % idx % this->_contents.size());
	this->_contents[idx] = matter;
	matter->setContainedBy(this->sharedThis<Matter_O>());
	matter->setId(idx);
    }

    string Matter_O::propertiesAsString() const
    {_OF();
	stringstream ss;
	ss << "Properties for matter[" << this->__repr__() << "]"<< std::endl;
        ss << _rep_(this->_Properties);
	return ss.str();
    }
	


    void Matter_O::clearProperty(core::Symbol_sp prop)
    {_G();
	core::core_rem_f(this->_Properties,prop);
    }

    void Matter_O::setProperty(core::Symbol_sp prop, core::T_sp val)
    {_G();
      core::core_put_f(this->_Properties,val,prop);
    }

    void Matter_O::setPropertyTrue(core::Symbol_sp prop)
    {_G();
      core::core_put_f(this->_Properties,_lisp->_true(),prop);
    }

    core::T_sp Matter_O::getProperty(core::Symbol_sp prop)
    {
      core::T_sp res = core::cl_getf(this->_Properties,prop,_Unbound<core::T_O>());
      if (res.unboundp()) {
        stringstream props;
        props << _rep_(this->_Properties);
        SIMPLE_ERROR(BF("You asked for an unknown property[%s] for matter[%s@%p] - the available properties are[%s]") % prop->__repr__() % this->__repr__() % this % props.str()  );
      }
      return res;
    }

    core::T_sp Matter_O::getPropertyOrDefault(core::Symbol_sp prop,core::T_sp defval)
    {_G();
      return core::cl_getf(this->_Properties,prop,defval);
    }

    bool Matter_O::hasProperty(core::Symbol_sp prop)
    {_G();
      return !core::cl_getf(this->_Properties,prop,_Unbound<core::T_O>()).unboundp();
    }



    Atom_sp	Matter_O::firstAtomWithName(MatterName name)
    {_OF();
	Loop l;
	l.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	while ( l.advance() )
	{
	    Atom_sp a = l.getAtom();
	    if ( a->getName() == name )
	    {
		return a;
	    }
	}
	SIMPLE_ERROR(BF("Could not find atom with name %s") % _rep_(name));
    }


    void Matter_O::calculateVirtualAtomPositions()
    {
	Loop l;
	l.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	while ( l.advance() )
	{
	    Atom_sp a = l.getAtom();
	    if ( a.isA<VirtualAtom_O>() )
	    {
		VirtualAtom_sp va = a.as<VirtualAtom_O>();
		va->calculatePosition();
	    }
	}
    }


    void Matter_O::connectAllCloseAtoms()
    {_G();
        gctools::Vec0<Atom_sp>	atoms;
	Loop		la;
	la.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	while ( la.advanceLoopAndProcess() )
	{
	    atoms.push_back(la.getAtom());
	}
        gctools::Vec0<Atom_sp>::iterator ai, aj;
	for ( ai=atoms.begin(); ai<atoms.end()-1; ai++ )
	{
	    Vector3 aipos = (*ai)->getPosition();
	    for ( aj=ai+1; aj!=atoms.end(); aj++ )
	    {
		Vector3 ajpos = (*aj)->getPosition();
		if ( fabs(ajpos.getX()-aipos.getX()) > 2.0 ) continue;
		if ( fabs(ajpos.getY()-aipos.getY()) > 2.0 ) continue;
		if ( fabs(ajpos.getZ()-aipos.getZ()) > 2.0 ) continue;
		if ( (*ai)->isBondedTo(*aj) ) continue;
		double distSquared = calculateDistanceSquared(aipos,ajpos);
		if ( distSquared > 4.0 ) continue;
		(*ai)->bondToSingle(*aj);
	    }
	}
    }








//
// contentWithName
//
    Matter_sp   Matter_O::contentWithName(MatterName    sName )
    {_G();
	contentIterator	aCur;

	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    LOG(BF("Looking at(%s) for(%s)") % (*aCur)->getName().c_str() % sName.c_str()  );
	    if ( (*aCur)->getName() == sName ) {
		return( (*aCur) );
	    }
	}
	LOG(BF("Matter(%s) with %d contents does not contain content with name(%s)") % this->name.c_str() % this->_contents.size() % sName.c_str()  );
	stringstream ss;
	ss << this->className() << " (" << _rep_(this->name) << ") does not contain name(" << _rep_(sName) << ")";
	SIMPLE_ERROR(BF(ss.str()));
    }


//
// contentWithName
//
    Matter_sp   Matter_O::contentWithNameOrNil(MatterName sName )
    {_G();
	contentIterator	aCur;

	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    LOG(BF("Looking at(%s) for(%s)") % (*aCur)->getName().c_str() % sName.c_str()  );
	    if ( (*aCur)->getName() == sName ) {
		return( (*aCur) );
	    }
	}
	return _Nil<Matter_O>();
    }



    void	Matter_O::makeAllAtomNamesInEachResidueUnique()
    {_G();
	contentIterator	aCur;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) 
	{
	    (*aCur)->makeAllAtomNamesInEachResidueUnique();
	}
    }

    void	Matter_O::fillInImplicitHydrogensOnCarbon()
    {_G();
	contentIterator	aCur;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) 
	{
	    (*aCur)->fillInImplicitHydrogensOnCarbon();
	}
    }

    void	Matter_O::randomizeAtomPositions()
    {_G();
	contentIterator	aCur;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) 
	{
	    (*aCur)->randomizeAtomPositions();
	}
    }


    void	Matter_O::perturbAtomPositions(double dist)
    {_G();
	contentIterator	aCur;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) 
	{
	    (*aCur)->perturbAtomPositions(dist);
	}
    }





//
// contentIndexWithName
//
    int	Matter_O::contentIndexWithName(MatterName sName )
    {_OF();
	contentIterator	aCur;
	int				i;
	i = 0;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++, i++ ) {
	    if ( (*aCur)->getName() == sName ) {
		return( i );
	    }
	}
	stringstream ss;
	ss << this->description() <<") does not contain name("<<_rep_(sName)<<")";
	SIMPLE_ERROR(BF("%s")%ss.str());
    }


//
// contentWithId
//
    Matter_sp   Matter_O::contentWithId( int lid )
    {_OF();
	contentIterator	aCur;
	Matter_sp			c;

	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    c = (*aCur);
	    LOG(BF("Looking at content with id(%d)") % c->getId() );
	    if ( c->getId() == lid ) {
		return( c );
	    }
	}
	// THROW an exception
	stringstream ss;
	ss << this->description() << ") does not contain id("<<lid<<")";
	SIMPLE_ERROR(BF("%s")%ss.str());
    }

//
// 	hasContentWithId
//
    bool	Matter_O::hasContentWithId( int lid )
    {_G();
	contentIterator	aCur;
	Matter_sp			c;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    c = (*aCur);
	    LOG(BF("Looking at content with id(%d)") % c->getId() );
	    if ( c->getId() == lid ) {
		return( true );
	    }
	}
	return( false);
    }


#if 0 //[
//
// contentWithStorageId
//
    Matter_sp   Matter_O::contentWithStorageId( int lid )
    {
	contentIterator	aCur;
	char				msg[255];
	char				str1[255];
	Matter_sp			c;

	LOG(BF("Looking in container(%s) type(%c) for content with storageId(%d)") % this->name.c_str() % this->containerType % lid );
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    c = (*aCur);
	    LOG(BF("Looking at content with id(%d)") % c->getId() );
	    if ( c->getStorageId() == lid ) {
		return( c );
	    }
	}

	// THROW an exception
	stringstream ss;
	ss << "Matter("<<this->name<<") type("<<this->containerType<<") does not contain id("<<lid<<")";
	SIMPLE_ERROR(BF("%s")%ss.str());
    }

//
// 	hasContentWithStorageId
//
    bool	Matter_O::hasContentWithStorageId( int lid )
    {
	contentIterator	aCur;
	Matter_sp			c;
	LOG(BF("Looking in container(%s) type(%c) for content with storageId(%d)") % this->name.c_str() % this->containerType % lid );
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++ ) {
	    c = (*aCur);
	    LOG(BF("Looking at content with id(%d)") % c->getId() );
	    if ( c->getStorageId() == lid ) {
		return( true );
	    }
	}
	return( false);
    }
#endif //]


//
// contentIndexWithId
//
    int	Matter_O::contentIndexWithId( int lid )
    {_OF();
	contentIterator	aCur;
	int				i;
	Matter_sp			c;
	i = 0;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++, i++ ) {
	    c = *aCur;
	    if ( c->getId() == lid ) {
		return( i );
	    }
	}
	SIMPLE_ERROR(BF("Matter %s does not contain content with id(%s)") % _rep_(this->name) % lid);
    }



//
// contentIndex
//
    int	Matter_O::contentIndex( Matter_sp cc)
    {_OF();
	contentIterator	aCur;
	int				i;
	i = 0;
	for ( aCur=this->_contents.begin();aCur!=this->_contents.end(); aCur++, i++ ) {
	    if ( (*aCur) == cc ) {
		return( i );
	    }
	}
        SIMPLE_ERROR(BF("Matter %s does not contain %s") % _rep_(this->name) % _rep_(cc->getName()));
    }


/*! Add the Matter as a child of this Matter
 */
    void	Matter_O::addMatter(Matter_sp cp )
    {_OF();
	Matter_sp	ctemp;
	ASSERTNOTNULL(cp);
	LOG(BF("Adding: %s of type: %c") % cp->getName().c_str() % cp->getMatterType()  );
	ctemp = this->sharedThis<Matter_O>();
	cp->setContainedBy(ctemp);
	LOG(BF("Adding %c:0x%08x to parent %c:0x%08x") % (cp->getMatterType()) % (cp.get()) % (this->getMatterType()) % (this ) );
    	// Always add the content to the end of the vector
	// A lot depends on Residues maintaining the order of Atoms
	// throughout the various passes of building databases
	this->_contents.push_back(cp);
	if ( this->_Id == UNDEFINED_ID )
	{
	    this->_Id = this->_NextContentId;
	    this->_NextContentId++;
	}
	LOG(BF("Finished adding") );
    }

//
// addMatterRetainId
//
// Add the container but retain the containers id
// This is used when we load a container from a file
//
    void	Matter_O::addMatterRetainId(Matter_sp cp )
    {_OF();
	Matter_sp	ctemp;
	LOG(BF("addMatterRetainId to %s adding container %s")
	    % this->description()
	    % cp->description() );
	ctemp = this->sharedThis<Matter_O>();
	cp->setContainedBy(ctemp);
	LOG(BF("Adding %c:0x%08x to parent %c:0x%08x") % (cp->getMatterType()) % (cp.get()) % (this->getMatterType()) % (this ) );
    	// Always add the content to the end of the vector
	// A lot depends on Residues maintaining the order of Atoms
	// throughout the various passes of building databases
	this->_contents.push_back(cp);
	if (this->_NextContentId <= cp->_Id ) 
	{
	    this->_NextContentId = cp->_Id+1;
	}
	LOG(BF("Finished adding") );
    }

//
//	setId
//
//	Set the id of the contained object and update
//	nextContentId if it is too small
//
    void	Matter_O::setId(int i)
    {
	if ( this->_NextContentId < i ) {
	    this->_NextContentId = i+1;
	}
	this->_Id = i;
    }


//
// applyTransformToAtoms
//


    void	Matter_O::applyTransformToAtoms( const Matrix& m )
    {_OF();
	contentIterator	a;
	LOG(BF("Transforming all atoms with:%s") % (m.asString().c_str() ) );
	for ( a=this->_contents.begin(); a!=this->_contents.end(); a++ ) {
	    (*a)->applyTransformToAtoms(m);
	}
    }




    void	Matter_O::invertStructureAndRestraints()
    {
	// First invert all of the contained objects and their restraints
	contentIterator	a;
	for ( a=this->_contents.begin(); a!=this->_contents.end(); a++ )
	{
	    (*a)->invertStructureAndRestraints();
	}
	this->invertStereochemistryOfRestraints();
    }

    void Matter_O::invertStereochemistryOfRestraints()
    {_OF();
	//
	//Then flip any restraints
	if ( this->_Restraints.notnilp() )
	{
	    RestraintList_O::iterator it;
	    for ( it=this->_Restraints->begin(); it!=this->_Restraints->end(); it++ )
	    {
		(*it)->invertStereochemistryOfRestraint();
	    }
	}
    }


    void	Matter_O::translateAllAtoms(Vector3 trans)
    {_G();
	Matrix transform;
	transform.translate(&trans);
	this->applyTransformToAtoms(transform);
    }


    bool	Matter_O::isContainedBy(Matter_sp container)
    {_G();
	Matter_sp outer = this->sharedThis<Matter_O>();
	while ( 1 )
	{
	    if ( outer == container ) return true;
	    if ( !outer->containedByValid() ) return false;
	    outer = outer->containedBy();
	}
    }


    bool	Matter_O::testConsistancy(const Matter_sp parentShouldBe )
    {_OF();
	contentIterator	a;
	Matter_sp			c;
	if ( this->containedByValid() ) {
	    if ( (parentShouldBe != Matter_sp(this->containedBy())) ) {
		LOG(BF("Matter_O::testConsistancy failed!") );
		LOG(BF("  My name/type = %s/%c") % this->getName().c_str() % this->getMatterType()  );
//	    LOG(BF("  My parent address = %x  was supposed to be: %x") % &(Matter_sp(this->containedBy))() % &parentShouldBe  );
	    }
	    return false;
	} else {
	    LOG(BF("  Parentage ok") );
	}
	c = this->sharedThis<Matter_O>();
	for ( a=this->_contents.begin(); a!=this->_contents.end(); a++ ) {
	    c  = *a;
	    if ( !(c->testConsistancy(this->sharedThis<Matter_O>())) ) {
		LOG(BF("Child failed consistancy test") );
		return false;
	    }
	}
	return true;
    }


    void	Matter_O::reparent(Matter_sp	newParent)
    {_OF();
	Matter_sp			oldParent;
	Matter_wp			wctemp;

//    oldParent = this->containedBy();
	this->setContainedBy(newParent);
	newParent->addMatter(this->sharedThis<Matter_O>());
    }


    int	Matter_O::totalNetResidueCharge()
    {_G();
	Loop	l;
	l.loopTopGoal(this->sharedThis<Matter_O>(),RESIDUES);
	int netCharge = 0;
	Residue_sp	r;
	while ( l.advanceLoopAndProcess() )
	{
	    r = l.getResidue();
	    netCharge += r->getNetCharge();
	}
	return netCharge;
    }



    void	Matter_O::setAllAtomMasks(int m)
    {_G();
	Loop	l;
	l.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	Atom_sp	a;
	while ( l.advanceLoopAndProcess() )
	{
	    a = l.getAtom();
	    a->setMask(m);
	}
    }

/*
  __BEGIN_DOC( candoScript.matter.setAtomAliasesForResiduesNamed, setAtomAliasesForResiduesNamed)
  \scriptCmd{extendAliases}{Cons::residuesAndInterestingAtomNames Cons::atomAliases}\par

  Lookup the residues in the Matter and set the atom aliases of the atoms.


  \begin{verbatim}
  setAtomAliasesForResiduesNamed (:
  (: 'glu(S) (aliasAtoms 'OE ) )
  (: 'glu(R) (aliasAtoms 'OE ) )
  )
  (atomAliases 'carbO )
  \end{verbatim}
  __END_DOC
*/
    void	Matter_O::setAtomAliasesForResiduesNamed(core::List_sp parts, core::List_sp atomAliases )
    {_G();
        FIX_ME();
#if 0
	if ( this->asSmartPtr().isA<Atom_O>() )
	{
	    SIMPLE_ERROR(BF("Needs Aggregate, Molecule or Residue"));
	}
	LOG(BF("setAtomAliasesForResiduesNamed with parts(%s) atomAliases(%s)") % parts->__repr__().c_str() % atomAliases->__repr__().c_str()  );
        gctools::SmallMap<MatterName,core::List_sp> residueNamesToAliasAtoms;
	for ( auto cur : parts ) {
          core::List_sp oneExtend = oCar(cur);
	    if ( oneExtend->length() != 2 )
	    {
		SIMPLE_ERROR(BF("Each entry must have two elements: "+oneExtend->__repr__() ));
	    }
	    core::Str_sp residueName = oneExtend->car<core::Str_O>();
	    core::List_sp aliasAtoms = oCadr(oneExtend);
	    LOG(BF("residueName(%s) aliasAtoms(%s)") % residueName->get().c_str() % aliasAtoms->__repr__().c_str() );
	    residueNamesToAliasAtoms[residueName->get()] = aliasAtoms;
	}
    	//
	// Now set the alias atoms of the residues
	//
	Loop lResidues;
	lResidues.loopTopGoal(this->sharedThis<Matter_O>(),RESIDUES);
	while ( lResidues.advanceLoopAndProcess() )
	{
	    Residue_sp res = lResidues.getResidue();
	    if ( residueNamesToAliasAtoms.count(res->getName())>0 )
	    {
		LOG(BF("%s:%d Setting aliases for residue(%s)") % (__FILE__) % (__LINE__) % (res->getName().c_str()) );
		core::List_sp aliasAtoms = residueNamesToAliasAtoms[res->getName()];
		LOG(BF("setting aliases for residueName(%s) aliasAtoms(%s) atomAliases(%s)") % res->getName().c_str() % aliasAtoms->__repr__().c_str() % atomAliases->__repr__().c_str() );
		res->setAliasesForAtoms(aliasAtoms,atomAliases);
	    } else
	    {
		LOG(BF("%s:%d Not setting aliases for residue(%s)") % (__FILE__) % (__LINE__) % (res->getName().c_str()) );
	    }
	}
#endif
    }




    Vector3	Matter_O::geometricCenter()
    {_OF();
	Vector3	sum;
	int	count;
	Atom_sp	a;
	Loop	l;
	Matter_sp	c;

	c = this->sharedThis<Matter_O>();
	sum.set(0,0,0);
	count = 0;
	l.loopTopGoal( c, ATOMS );
	while ((l.advanceLoopAndProcess())) {
	    a = l.getAtom();
	    sum = sum + a->getPosition();
	    count++;
	}
	LOG(BF("geometricCenterOfMatter number of atoms = %d") % count  );
	if ( count != 0 ) {
	    sum = sum*(1.0/count);
	}
	LOG(BF("Geometric center = %f, %f, %f") % sum.getX() % sum.getY() % sum.getZ()  );
	return sum;
    }


    geom::BoundingBox_sp Matter_O::boundingBox(double pad) 
    {_G();
	Vector3	sum;
	int	count;
	Atom_sp	a;
	Loop	l;
	Matter_sp	c;

	c = this->sharedThis<Matter_O>();
	sum.set(0,0,0);
	count = 0;
	geom::BoundingBox_sp bbox = geom::BoundingBox_O::create();
	l.loopTopGoal( c, ATOMS );
	while ((l.advanceLoopAndProcess())) {
	    a = l.getAtom();
	    bbox->expandToEncompassPoint(a->getPosition());
	    count++;
	}
	bbox->pad(pad);
	return bbox;
    }


    core::List_sp Matter_O::contentsAsCons()
    {_G();
	core::List_sp cur = _Nil<core::T_O>();
	contentIterator	it;
	for ( it=this->_contents.end()-1; it>=this->_contents.begin(); it-- )
	{
	    cur = core::Cons_O::create(*it,cur);
	}
	return cur;
    }



    core::Vector_sp	Matter_O::allAtoms()
    {_G();
        core::VectorObjectsWithFillPtr_sp	atoms = core::VectorObjectsWithFillPtr_O::create();
	Loop		la;
	la.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	while ( la.advanceLoopAndProcess() )
	{
	    atoms->vectorPushExtend(la.getAtom());
	}
	return atoms;
    }


    core::Vector_sp	Matter_O::allAtomsOfElement(Element element)
    {_G();
        core::VectorObjectsWithFillPtr_sp atoms = core::VectorObjectsWithFillPtr_O::create();
	Loop		la;
	Atom_sp		a;
	la.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	while ( la.advanceLoopAndProcess() )
	{
	    a = la.getAtom();
	    if ( a->getElement() == element )
	    {
		atoms->vectorPushExtend(a);
	    }
	}
	return atoms;
    }


    Residue_sp	Matter_O::aliasResidueOrNil(Alias_sp alias)
    {_G();
	Loop lRes;
	lRes.loopTopGoal(this->sharedThis<Matter_O>(),RESIDUES);
	while ( lRes.advance() )
	{
	    if ( lRes.getResidue()->recognizesMonomerAlias(alias->getMonomerAlias()) )
	    {
		return lRes.getResidue();
	    }
	}
	return _Nil<Residue_O>();
    }

    Residue_sp Matter_O::aliasResidue(Alias_sp alias)
    {_OF();
	Residue_sp a = this->aliasResidueOrNil(alias);
	if ( a.notnilp() ) return a;
	SIMPLE_ERROR(BF("could not find aliasResidue"));
    }

    Atom_sp	Matter_O::aliasAtomOrNil(Alias_sp alias)
    {
	Residue_sp res = this->aliasResidueOrNil(alias);
	if ( res.nilp() ) return _Nil<Atom_O>();
	return res->atomWithAliasOrNil(alias->getAtomAlias());
    }

    Atom_sp Matter_O::aliasAtom(Alias_sp alias)
    {_OF();
	Atom_sp a = this->aliasAtomOrNil(alias);
	if ( a.notnilp() ) return a;
	SIMPLE_ERROR(BF("could not find aliasAtom"));
    }





    core::List_sp Matter_O::allAtomsOfElementAsList(Element element)
    {_G();
      core::List_sp		list(_Nil<core::T_O>());
	Loop		la;
	Atom_sp		a;
	la.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
	while ( la.advanceLoopAndProcess() ) {
	    a = la.getAtom();
	    if ( a->getElement() == element )
	    {
		list = core::Cons_O::create(a,list);
	    }
	}
	return list;
    }

void	Matter_O::fields(core::Record_sp node )
{_G();
  node->field( INTERN_(kw,name), this->name);
  node->/*pod_*/field_if_not_default( INTERN_(kw,id), this->_Id, 0);
  node->field_if_not_nil( INTERN_(kw,restraints),this->_Restraints);
  node->field_if_not_nil( INTERN_(kw,properties),this->_Properties);
  node->field_if_not_empty( INTERN_(kw,contents), this->_contents);
  LOG(BF("Status") );
#if 0
  if ( node->loading() ) {
    _BLOCK_TRACEF(BF("serializing container contents - there are %d objects")% this->_contents.size() );
	    // Make sure all contents have us as a parent
    Matter_sp c = this->sharedThis<Matter_O>();
    gctools::Vec0<Matter_sp>::iterator ai;
    for ( ai=this->_contents.begin(); ai!=this->_contents.end(); ai++ ) {
      (*ai)->setContainedBy(c);
    }
  }
#endif
  LOG(BF("Status") );
}






/*! Subclass this method
 */
#ifdef RENDER
    geom::Render_sp Matter_O::rendered(core::List_sp kopts)
    {
	return _Nil<geom::Render_O>();
    };
#endif



    core::List_sp Matter_O::allAtomsAsCons(bool allowVirtualAtoms ) const
    {_OF();
	core::List_sp result = _Nil<core::T_O>();
	Loop l(this->const_sharedThis<Matter_O>(),ATOMS);
	while ( l.advance() )
	{
	    Atom_sp a1 = l.getAtom();
	    if ( !allowVirtualAtoms && a1->isVirtual() ) continue;
	    result = core::Cons_O::create(a1,result);
	}
	return result;
    }



    core::List_sp Matter_O::allBondsAsCons(bool allowVirtualAtoms ) const
    {_OF();
	core::List_sp result = _Nil<core::List_V>();
	Loop l(this->const_sharedThis<Matter_O>(),BONDS);
	while ( l.advance() )
	{
	    Bond_sp bond = l.getBond();
	    LOG(BF("Adding %s") % bond->description() );
	    Atom_sp a1 = bond->getAtom1();
	    Atom_sp a2 = bond->getAtom2();
	    if ( !allowVirtualAtoms && ( a1->isVirtual() || a2->isVirtual() ) ) continue;
	    result = core::Cons_O::create(bond,result);
	}
	return result;
    }



    core::List_sp Matter_O::allAnglesAsCons(bool allowVirtualAtoms ) const
    {_OF();
	core::List_sp result = _Nil<core::T_O>();
	Loop l(this->const_sharedThis<Matter_O>(),ANGLES);
	while ( l.advance() )
	{
	    Atom_sp a1 = l.getAtom1();
	    Atom_sp a2 = l.getAtom2();
	    Atom_sp a3 = l.getAtom3();
	    if ( !allowVirtualAtoms && ( a1->isVirtual()
					 || a2->isVirtual()
					 || a3->isVirtual()
		     ) ) continue;
	    Angle_sp entry = Angle_O::create(a1,a2,a3);
	    result = core::Cons_O::create(entry,result);
	}
	return result;
    }


    core::List_sp Matter_O::allImproperTorsionsAsCons(bool allowVirtualAtoms ) const
    {_OF();
	core::List_sp result = _Nil<core::List_V>();
	Loop l(this->const_sharedThis<Matter_O>(),IMPROPERS);
	{_BLOCK_TRACE("Iterating over loop");
	    int count=0;
	    while ( l.advance() )
	    {
		count++;
		Atom_sp a1 = l.getAtom1();
		Atom_sp a2 = l.getAtom2();
		Atom_sp a3 = l.getAtom3();
		Atom_sp a4 = l.getAtom4();
		if ( !allowVirtualAtoms && ( a1->isVirtual()
					     || a2->isVirtual()
					     || a3->isVirtual()
					     || a4->isVirtual() ) ) continue;
		ImproperTorsion_sp entry = ImproperTorsion_O::create(a1,a2,a3,a4);
		result = core::Cons_O::create(entry,result);
	    }
	    LOG(BF("There were %d improper torsions iterated over") % count );
	}
	return result;
    }




    core::List_sp Matter_O::allProperTorsionsAsCons(bool allowVirtualAtoms ) const
    {_OF();
	core::List_sp result = _Nil<core::T_O>();
	Loop l(this->const_sharedThis<Matter_O>(),PROPERS);
	while ( l.advance() )
	{
	    Atom_sp a1 = l.getAtom1();
	    Atom_sp a2 = l.getAtom2();
	    Atom_sp a3 = l.getAtom3();
	    Atom_sp a4 = l.getAtom4();
	    if ( !allowVirtualAtoms && ( a1->isVirtual()
					 || a2->isVirtual()
					 || a3->isVirtual()
					 || a4->isVirtual() ) ) continue;
	    ProperTorsion_sp entry = ProperTorsion_O::create(a1,a2,a3,a4);
	    result = core::Cons_O::create(entry,result);
	}
	return result;
    }




    Matter_sp Matter_O::copyDontRedirectAtoms()
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }

    void Matter_O::redirectAtoms()
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    void Matter_O::copyRestraintsDontRedirectAtoms(Matter_sp orig)
    {_OF();
	if ( orig->_Restraints.nilp() ) {
          this->_Restraints = _Nil<core::T_O>();
	} else {
          this->_Restraints = orig->_Restraints->copyDontRedirectAtoms();
	}
    }


    void Matter_O::redirectRestraintAtoms()
    {_OF();
	if ( this->_Restraints.notnilp() )
	{
	    this->_Restraints->redirectAtoms();
	}
    }


    AtomIdToAtomMap_sp Matter_O::buildAtomIdMap() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }

    Atom_sp Matter_O::atomWithAtomId(AtomId_sp atomId) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }






    void Matter_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Matter_O>()
//	    .def_raw("core:__init__",&Matter_O::__init__,"(self &key name)")
          .def("matter-copy",&Matter_O::copy)
	    .def( "setAtomAliasesForResiduesNamed", &Matter_O::setAtomAliasesForResiduesNamed)
	    .def("translateAllAtoms",&Matter_O::translateAllAtoms)
	    .def("connectAllCloseAtoms",&Matter_O::connectAllCloseAtoms)
	    .def("aliasAtom",&Matter_O::aliasAtom)
	    .def("aliasResidue",&Matter_O::aliasResidue)
	    .def("aliasAtomOrNil",&Matter_O::aliasAtomOrNil)
	    .def("aliasResidueOrNil",&Matter_O::aliasResidueOrNil)
//	    .def("getMatterType",&Matter_O::getMatterType)
	    .def("getId",&Matter_O::getId)
	    .def("numberOfAtoms",&Matter_O::numberOfAtoms)
	    .def("setName",&Matter_O::setName,"","","",false)
	    .def("getName",&Matter_O::getName_notConst,"","","",false)
	    .def("isAggregate",&Matter_O::isAggregate)
	    .def("isMolecule",&Matter_O::isMolecule)
	    .def("isResidue",&Matter_O::isResidue)
	    .def("isAtom",&Matter_O::isAtom)
	    .def("boundingBox",&Matter_O::boundingBox)
	    .def("contentsAsCons",&Matter_O::contentsAsCons)
	    .def("allAtoms",&Matter_O::allAtoms)
	    .def("invertStructureAndRestraints",&Matter_O::invertStructureAndRestraints)
	    .def("allAtomsOfElementAsList",&Matter_O::allAtomsOfElementAsList)
//	.def("get_StorageId",&Matter_O::get_StorageId)
	    .def("firstAtomWithName",&Matter_O::firstAtomWithName)
	    .def("addMatter",&Matter_O::addMatter)
          .def("containedBy",(Matter_sp (Matter_O::*)() const)&Matter_O::containedBy)
	    .def("setAllAtomMasks",&Matter_O::setAllAtomMasks)
	    .def("hasContentWithName",&Matter_O::hasContentWithName)
	    .def("contentWithName",&Matter_O::contentWithName)
	    .def("makeAllAtomNamesInEachResidueUnique",&Matter_O::makeAllAtomNamesInEachResidueUnique)
	    .def("contentIndexWithName",&Matter_O::contentIndexWithName)
	    .def("contentWithId",&Matter_O::contentWithId)
//	.def("contentWith_StorageId",&Matter_O::contentWith_StorageId)
	    .def("hasContentWithId",&Matter_O::hasContentWithId)
//	.def("hasContentWith_StorageId",&Matter_O::hasContentWith_StorageId)
	    .def("contentIndexWithId",&Matter_O::contentIndexWithId)
	    .def("contentIndex",&Matter_O::contentIndex)
	    .def("contentAt",&Matter_O::contentAt)
	    .def("contentSize",&Matter_O::contentSize)
	    .def("randomizeAtomPositions",&Matter_O::randomizeAtomPositions)
	    .def("perturbAtomPositions",&Matter_O::perturbAtomPositions)
	    .def("fillInImplicitHydrogensOnCarbon",&Matter_O::fillInImplicitHydrogensOnCarbon)
	    .def("applyTransformToAtoms",&Matter_O::applyTransformToAtoms)
	    .def("testConsistancy",&Matter_O::testConsistancy)
	    .def("geometricCenter", &Matter_O::geometricCenter )
//	    .def("getRestraints",&Matter_O::getRestraints)
	    .def("calculateVirtualAtomPositions",&Matter_O::calculateVirtualAtomPositions)
	    .def("clearProperty",&Matter_O::clearProperty)
	    .def("hasProperty",&Matter_O::hasProperty)
	    .def("setProperty",&Matter_O::setProperty)
	    .def("setPropertyTrue",&Matter_O::setPropertyTrue)
	    .def("Matter-getProperty",&Matter_O::getProperty)
	    .def("Matter-getPropertyOrDefault",&Matter_O::getPropertyOrDefault)
	    .def("propertiesAsString",&Matter_O::propertiesAsString)
	    .def("allAtomsAsCons",&Matter_O::allAtomsAsCons)
	    .def("allBondsAsCons",&Matter_O::allBondsAsCons)
	    .def("allAnglesAsCons",&Matter_O::allAnglesAsCons)
	    .def("allImproperTorsionsAsCons",&Matter_O::allImproperTorsionsAsCons)
	    .def("allProperTorsionsAsCons",&Matter_O::allProperTorsionsAsCons)
	    .def("buildAtomIdMap",&Matter_O::buildAtomIdMap)
	    .def("atomWithAtomId",&Matter_O::atomWithAtomId)
	    ;
    }

    void Matter_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Matter,"","",_lisp)
	    .def( "setAtomAliasesForResiduesNamed", &Matter_O::setAtomAliasesForResiduesNamed)
	    .def("translateAllAtoms",&Matter_O::translateAllAtoms)
	    .def("connectAllCloseAtoms",&Matter_O::connectAllCloseAtoms)
	    .def("aliasAtom",&Matter_O::aliasAtom)
	    .def("aliasResidue",&Matter_O::aliasResidue)
	    .def("aliasAtomOrNil",&Matter_O::aliasAtomOrNil)
	    .def("aliasResidueOrNil",&Matter_O::aliasResidueOrNil)
	    .def("getMatterType",&Matter_O::getMatterType)
	    .def("getId",&Matter_O::getId)
	    .def("numberOfAtoms",&Matter_O::numberOfAtoms)
	    .def("setName",&Matter_O::setName)
	    .def("getName",&Matter_O::getName_notConst)
	    .def("isAggregate",&Matter_O::isAggregate)
	    .def("isMolecule",&Matter_O::isMolecule)
	    .def("isResidue",&Matter_O::isResidue)
	    .def("isAtom",&Matter_O::isAtom)
	    .def("boundingBox",&Matter_O::boundingBox)
	    .def("contentsAsCons",&Matter_O::contentsAsCons)
	    .def("allAtoms",&Matter_O::allAtoms)
	    .def("invertStructureAndRestraints",&Matter_O::invertStructureAndRestraints)
	    .def("allAtomsOfElementAsList",&Matter_O::allAtomsOfElementAsList)
//	.def("get_StorageId",&Matter_O::get_StorageId)
	    .def("firstAtomWithName",&Matter_O::firstAtomWithName)
	    .def("addMatter",&Matter_O::addMatter)
	    .def("containedBy",&Matter_O::containedByLock)
	    .def("setAllAtomMasks",&Matter_O::setAllAtomMasks)
	    .def("hasContentWithName",&Matter_O::hasContentWithName)
	    .def("contentWithName",&Matter_O::contentWithName)
	    .def("makeAllAtomNamesInEachResidueUnique",&Matter_O::makeAllAtomNamesInEachResidueUnique)
	    .def("contentIndexWithName",&Matter_O::contentIndexWithName)
	    .def("contentWithId",&Matter_O::contentWithId)
//	.def("contentWith_StorageId",&Matter_O::contentWith_StorageId)
	    .def("hasContentWithId",&Matter_O::hasContentWithId)
//	.def("hasContentWith_StorageId",&Matter_O::hasContentWith_StorageId)
	    .def("contentIndexWithId",&Matter_O::contentIndexWithId)
	    .def("contentIndex",&Matter_O::contentIndex)
	    .def("contentAt",&Matter_O::contentAt)
	    .def("contentSize",&Matter_O::contentSize)
	    .def("randomizeAtomPositions",&Matter_O::randomizeAtomPositions)
	    .def("perturbAtomPositions",&Matter_O::perturbAtomPositions)
	    .def("fillInImplicitHydrogensOnCarbon",&Matter_O::fillInImplicitHydrogensOnCarbon)
	    .def("applyTransformToAtoms",&Matter_O::applyTransformToAtoms)
	    .def("testConsistancy",&Matter_O::testConsistancy)
	    .def("geometricCenter", &Matter_O::geometricCenter )
//	    .def("getRestraints",&Matter_O::getRestraints)
	    .def("calculateVirtualAtomPositions",&Matter_O::calculateVirtualAtomPositions)
	    .def("clearProperty",&Matter_O::clearProperty)
	    .def("hasProperty",&Matter_O::hasProperty)
	    .def("setProperty",&Matter_O::setProperty)
	    .def("setPropertyTrue",&Matter_O::setPropertyTrue)
	    .def("getProperty",&Matter_O::getProperty)
	    .def("getPropertyOrDefault",&Matter_O::getPropertyOrDefault)
	    .def("propertiesAsString",&Matter_O::propertiesAsString)
	    .def("allAtomsAsCons",&Matter_O::allAtomsAsCons)
	    .def("allBondsAsCons",&Matter_O::allBondsAsCons)
	    .def("allAnglesAsCons",&Matter_O::allAnglesAsCons)
	    .def("allImproperTorsionsAsCons",&Matter_O::allImproperTorsionsAsCons)
	    .def("allProperTorsionsAsCons",&Matter_O::allProperTorsionsAsCons)
	    .def("buildAtomIdMap",&Matter_O::buildAtomIdMap)
	    .def("atomWithAtomId",&Matter_O::atomWithAtomId)
	    ;
#endif
    }



    EXPOSE_CLASS(chem,Matter_O);

};  // namespace chem





