/*
    File: matter.cc
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




#include <string.h>
#include <iostream>
#include <iomanip>
#include <vector>
//#include	<strings>
#include <stdio.h>
#include <clasp/core/common.h>
#include <clasp/core/hashTableEq.h>
#include <cando/adapt/stringSet.h>
#include <cando/chem/matter.h>
#include <cando/geom/omatrix.h>
//#include "core/serialize.h"
#include <clasp/core/array.h>
#include <cando/chem/loop.h>
#include <clasp/core/cons.h>
#include <clasp/core/lispList.h>
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
#include <clasp/core/symbolTable.h>
#include <clasp/core/translators.h>
#include <clasp/core/wrappers.h>
#include <clasp/core/ql.h>

namespace chem
{

//
// Constructor
//


CL_DEFMETHOD Matter_sp Matter_O::contentAt(size_t i) const
{
  BOUNDS_ASSERT(i<this->_Contents.size());
  return this->_Contents[i];
}

string Matter_O::__repr__() const
{
  stringstream ss;
  ss << "#<" << this->className() << " " << _rep_(this->getName());
  ss << " :id " << this->_Id;
  ss << " 0x" << std::setbase(16) << core::lisp_general_badge(this->asSmartPtr());
  ss << ">";
  return ss.str();
}
//
// copy constructor
//
Matter_O::Matter_O( const Matter_O& c ) : core::CxxObject_O(c)
{
//  this->_NextContentId = c._NextContentId;
//  this->_Id = c._Id;
  this->_Name = c._Name;
//	this->containerContainedBy = c.containerContainedBy;
  this->_Properties = c._Properties; // You can't call allocators in ctors core::cl__copy_list(c._Properties);
  this->_RestraintList = nil<T_O>();
}

//
// eraseContents
//
void	Matter_O::eraseContents()
{
  this->_Contents.clear();
}

CL_LISPIFY_NAME("hasContent");
CL_DEFMETHOD bool	Matter_O::hasContent(Matter_sp maybe_child)
{
  contentIterator		aCur;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    LOG("Looking at(%s) for(%s)" , (*aCur)->getName().c_str() , sName.c_str()  );
    if ( (*aCur) == maybe_child ) {
      return true;
    }
  }
  return false;
}
//
// hasContentWithName
//
CL_LISPIFY_NAME("hasContentWithName");
CL_DEFMETHOD bool	Matter_O::hasContentWithName(MatterName    sName )
{
  contentIterator		aCur;

  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    LOG("Looking at(%s) for(%s)" , (*aCur)->getName().c_str() , sName.c_str()  );
    if ( (*aCur)->getName() == sName ) {
      return true;
    }
  }
  return false;
}


core::List_sp Matter_O::accumulateRestraints(core::List_sp allRestraints) const
{_OF();
  for ( auto cur : this->_RestraintList ) {
    Restraint_sp one = gc::As<Restraint_sp>(CONS_CAR(cur));
    allRestraints = core::Cons_O::create(one,allRestraints);
  }
  for ( const_contentIterator ci=this->begin_contents(); ci!=this->end_contents(); ci++ ) {
    allRestraints = (*ci)->accumulateRestraints(allRestraints);
  }
  return allRestraints;
}


CL_LISPIFY_NAME("allRestraints");
CL_DEFMETHOD core::List_sp Matter_O::allRestraints() const
{
  return this->accumulateRestraints(nil<core::T_O>());
}

CL_LISPIFY_NAME("clearRestraints");
CL_DEFMETHOD void Matter_O::clearRestraints()
{
  this->_RestraintList = nil<core::T_O>();
}

CL_LISPIFY_NAME("addRestraint");
CL_DEFMETHOD void Matter_O::addRestraint(Restraint_sp restraint)
{_OF();
  this->_RestraintList = core::Cons_O::create(restraint,this->_RestraintList);
}



CL_DEFMETHOD void Matter_O::resizeContents(int sz)
{_OF();
  this->_Contents.resize(sz,unbound<Matter_O>());
}


CL_DEFMETHOD void Matter_O::putMatter(int idx, Matter_sp matter)
{_OF();
  ASSERTF(idx>=0 && idx<(int)this->_Contents.size(),("Illegal putMatter index[%d] must be less than %d") , idx , this->_Contents.size());
  this->_Contents[idx] = matter;
}

CL_LISPIFY_NAME("propertiesAsString");
CL_DEFMETHOD string Matter_O::propertiesAsString() const
{_OF();
  stringstream ss;
  ss << "Properties for matter[" << this->__repr__() << "]"<< std::endl;
  ss << _rep_(this->_Properties);
  return ss.str();
}


bool Matter_O::applyPropertyToSlot(core::Symbol_sp prop, core::T_sp value)
{
  if ( prop == kw::_sym_name ) {
    this->_Name = gc::As<core::Symbol_sp>(value);
    return true;
  }
  return false;
}

void Matter_O::applyProperty(core::Symbol_sp prop, core::T_sp value)
{
//  printf("%s:%d applyProperty %s - %s\n", __FILE__, __LINE__, _rep_(prop).c_str(), _rep_(value).c_str());
  if ( !this->applyPropertyToSlot(prop,value) ) {
//    printf("%s:%d    using setProperty\n",__FILE__, __LINE__);
    this->setProperty(prop,value);
  }
}


void Matter_O::applyPropertyList(core::List_sp list)
{
  core::List_sp cur = list;
  while (cur.notnilp()) {
    core::Symbol_sp property = gc::As<core::Symbol_sp>(oCar(cur));
    cur = oCdr(cur);
    if ( cur.nilp() ) {
      SIMPLE_ERROR(("Incomplete property list: %s") , _rep_(list));
    }
    core::T_sp value = oCar(cur);
    this->applyProperty(property,value);
    cur = oCdr(cur);
  }
}

CL_LISPIFY_NAME("clearProperty");
CL_DEFMETHOD void Matter_O::clearProperty(core::Symbol_sp symbol)
{
  this->_Properties = core::core__rem_f(this->_Properties,symbol);
}

CL_DOCSTRING(R"dx(Set the property **symbol** of **this** (a chem:matter) to **value**.)dx")
CL_DEFMETHOD void Matter_O::setProperty(core::Symbol_sp symbol, core::T_sp value)
{
  this->_Properties = core::core__put_f(this->_Properties,value,symbol);
}

CL_DOCSTRING(R"dx(Set the property **symbol** of **this** (a chem:matter) to T.)dx")
CL_LISPIFY_NAME("setPropertyTrue");
CL_DEFMETHOD void Matter_O::setPropertyTrue(core::Symbol_sp symbol)
{
  this->_Properties = core::core__put_f(this->_Properties,_lisp->_true(),symbol);
}


CL_DOCSTRING(R"dx(Return the property **symbol** of **this** (a chem:matter) - if it isn't defined return NIL.)dx")
CL_LISPIFY_NAME("Matter-getProperty");
CL_DEFMETHOD core::T_sp Matter_O::getProperty(core::Symbol_sp symbol)
{
  core::T_sp res = core::cl__getf(this->_Properties,symbol,unbound<core::T_O>());
  if (res.unboundp()) {
    stringstream props;
    props << _rep_(this->_Properties);
    SIMPLE_ERROR(("You asked for an unknown property[%s] for matter[%s@%p] - the available properties are[%s]") , _rep_(symbol) , this->__repr__() , (void*)this , props.str()  );
  }
  return res;
}

CL_DOCSTRING(R"dx(Return the property **symbol** of **this** (a chem:matter) - if it isn't defined return **defval**.)dx")
CL_LISPIFY_NAME("Matter-getPropertyOrDefault");
CL_DEFMETHOD core::T_sp Matter_O::getPropertyOrDefault(core::Symbol_sp prop,core::T_sp defval)
{
  return core::cl__getf(this->_Properties,prop,defval);
}

CL_DOCSTRING(R"dx(Return T if the property **symbol** of **this** (a chem:matter) is defined.)dx")
CL_LISPIFY_NAME("hasProperty");
CL_DEFMETHOD bool Matter_O::hasProperty(core::Symbol_sp symbol)
{
  return !core::cl__getf(this->_Properties,symbol,unbound<core::T_O>()).unboundp();
}


CL_LISPIFY_NAME("firstAtomWithName");
CL_DEFMETHOD Atom_sp	Matter_O::firstAtomWithName(MatterName name)
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
  SIMPLE_ERROR(("Could not find atom with name %s") , _rep_(name));
}


CL_LISPIFY_NAME("calculateVirtualAtomPositions");
CL_DEFMETHOD void Matter_O::calculateVirtualAtomPositions()
{
  Loop l;
  Loop lres;
  lres.loopTopGoal(this->sharedThis<Matter_O>(),RESIDUES);
  while (lres.advance()) {
    Residue_sp res = lres.getResidue();
    l.loopTopGoal(res,ATOMS);
    while ( l.advance() ) {
      Atom_sp a = l.getAtom();
      if (gc::IsA<VirtualAtom_sp>(a)) {
        VirtualAtom_sp va = gc::As_unsafe<VirtualAtom_sp>(a);
        va->calculatePosition(res);
      }
    }
  }
}


CL_LISPIFY_NAME("connectAllCloseAtoms");
CL_DEFMETHOD void Matter_O::connectAllCloseAtoms()
{
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
      double distSquared = geom::calculateDistanceSquared(aipos,ajpos);
      if ( distSquared > 4.0 ) continue;
      (*ai)->bondToSingle(*aj);
    }
  }
}








//
// contentWithName
//
CL_LISPIFY_NAME("contentWithName");
CL_DEFMETHOD Matter_sp   Matter_O::contentWithName(MatterName    sName )
{
  contentIterator	aCur;

  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    LOG("Looking at(%s) for(%s)" , (*aCur)->getName().c_str() , sName.c_str()  );
    if ( (*aCur)->getName() == sName ) {
      return( (*aCur) );
    }
  }
  LOG("Matter(%s) with %d contents does not contain content with name(%s)" , this->_Name.c_str() , this->_Contents.size() , sName.c_str()  );
  stringstream ss;
  ss << this->className() << " (" << _rep_(this->_Name) << ") does not contain name(" << _rep_(sName) << ")";
  SIMPLE_ERROR((ss.str()));
}


//
// contentWithName
//
CL_DEFMETHOD core::T_sp   Matter_O::contentWithNameOrNil(MatterName sName )
{
  contentIterator	aCur;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    if ( (*aCur)->getName() == sName ) return( (*aCur) );
  }
  return nil<core::T_O>();
}



CL_LISPIFY_NAME("makeAllAtomNamesInEachResidueUnique");
CL_DEFMETHOD void	Matter_O::makeAllAtomNamesInEachResidueUnique()
{
  contentIterator	aCur;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ )
  {
    (*aCur)->makeAllAtomNamesInEachResidueUnique();
  }
}

CL_LISPIFY_NAME("fillInImplicitHydrogens");
CL_DEFMETHOD size_t	Matter_O::fillInImplicitHydrogens()
{
  contentIterator	aCur;
  size_t count = 0;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ )
  {
    count += (*aCur)->fillInImplicitHydrogens();
  }
  return count;
}

CL_LISPIFY_NAME("randomizeAtomPositions");
CL_DEFMETHOD void	Matter_O::randomizeAtomPositions()
{
  contentIterator	aCur;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ )
  {
    (*aCur)->randomizeAtomPositions();
  }
}


CL_LISPIFY_NAME("perturbAtomPositions");
CL_DEFMETHOD void	Matter_O::perturbAtomPositions(double dist)
{
  contentIterator	aCur;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ )
  {
    (*aCur)->perturbAtomPositions(dist);
  }
}





//
// contentIndexWithName
//
CL_LISPIFY_NAME("contentIndexWithName");
CL_DEFMETHOD int	Matter_O::contentIndexWithName(MatterName sName )
{_OF();
  contentIterator	aCur;
  int				i;
  i = 0;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++, i++ ) {
    if ( (*aCur)->getName() == sName ) {
      return( i );
    }
  }
  stringstream ss;
  ss << this->description() <<") does not contain name("<<_rep_(sName)<<")";
  SIMPLE_ERROR(("%s") , ss.str());
}

//
// contentWithId
//
CL_LISPIFY_NAME("contentWithId");
CL_DEFMETHOD Matter_sp   Matter_O::contentWithId( int lid )
{_OF();
  contentIterator	aCur;
  Matter_sp			c;

  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    c = (*aCur);
    LOG("Looking at content with id(%d)" , c->getId() );
    if ( c->getId() == lid ) {
      return( c );
    }
  }
	// THROW an exception
  stringstream ss;
  ss << this->description() << ") does not contain id("<<lid<<")";
  SIMPLE_ERROR(("%s") , ss.str());
}

//
// 	hasContentWithId
//
CL_LISPIFY_NAME("hasContentWithId");
CL_DEFMETHOD bool	Matter_O::hasContentWithId( int lid )
{
  contentIterator	aCur;
  Matter_sp			c;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    c = (*aCur);
    LOG("Looking at content with id(%d)" , c->getId() );
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

  LOG("Looking in container(%s) type(%c) for content with storageId(%d)" , this->_Name.c_str() , this->containerType , lid );
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    c = (*aCur);
    LOG("Looking at content with id(%d)" , c->getId() );
    if ( c->getStorageId() == lid ) {
      return( c );
    }
  }

	// THROW an exception
  stringstream ss;
  ss << "Matter("<<this->_Name<<") type("<<this->containerType<<") does not contain id("<<lid<<")";
  SIMPLE_ERROR(("%s") , ss.str());
}

//
// 	hasContentWithStorageId
//
bool	Matter_O::hasContentWithStorageId( int lid )
{
  contentIterator	aCur;
  Matter_sp			c;
  LOG("Looking in container(%s) type(%c) for content with storageId(%d)" , this->_Name.c_str() , this->containerType , lid );
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++ ) {
    c = (*aCur);
    LOG("Looking at content with id(%d)" , c->getId() );
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
CL_LISPIFY_NAME("contentIndexWithId");
CL_DEFMETHOD int	Matter_O::contentIndexWithId( int lid )
{_OF();
  contentIterator	aCur;
  int				i;
  Matter_sp			c;
  i = 0;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++, i++ ) {
    c = *aCur;
    if ( c->getId() == lid ) {
      return( i );
    }
  }
  SIMPLE_ERROR(("Matter %s does not contain content with id(%s)") , _rep_(this->_Name) , lid);
}



//
// contentIndex
//
CL_LISPIFY_NAME("contentIndex");
CL_DEFMETHOD int	Matter_O::contentIndex( Matter_sp cc)
{_OF();
  contentIterator	aCur;
  int				i;
  i = 0;
  for ( aCur=this->begin_contents();aCur!=this->end_contents(); aCur++, i++ ) {
    if ( (*aCur) == cc ) {
      return( i );
    }
  }
  SIMPLE_ERROR(("Matter %s does not contain %s") , _rep_(this->_Name) , _rep_(cc->getName()));
}


/*! Return the next _Id we would assign a child
 */
size_t Matter_O::nextId() const {
  size_t maxId = 0;
  for ( size_t idx = 0; idx<this->_Contents.size(); idx++ ) {
    size_t id = this->_Contents[idx]->_Id;
    if ( id > maxId ) maxId = id;
  }
  return maxId+1; // _Id's start counting at 1
}



/*! Add the Matter as a child of this Matter
 */
CL_LISPIFY_NAME("addMatter");
CL_DEFMETHOD Matter_mv	Matter_O::addMatter(Matter_sp cp )
{_OF();
  Matter_sp	ctemp;
  LOG("Adding: %s of type: %c" , cp->getName().c_str() , cp->getMatterType()  );
  ctemp = this->sharedThis<Matter_O>();
//  cp->setContainedBy(ctemp);
  LOG("Adding %c:0x%08x to parent %c:0x%08x" , (cp->getMatterType()) , (cp.get()) , (this->getMatterType()) , (this ) );
    	// Always add the content to the end of the vector
	// A lot depends on Residues maintaining the order of Atoms
	// throughout the various passes of building databases
  cp->_Id = this->nextId(); // Advance the _Id
  size_t index = this->_Contents.size();
  this->_Contents.push_back(cp);
  return Values(cp,core::make_fixnum(index));
}

//
// addMatterRetainId
//
// Add the container but retain the containers id
// This is used when we load a container from a file
//
void	Matter_O::addMatterRetainId(Matter_sp cp )
{_OF();
  IMPLEMENT_MEF("Get away from ids");
  Matter_sp	ctemp;
  LOG("addMatterRetainId to %s adding container %s"
      , this->description()
      , cp->description() );
  ctemp = this->sharedThis<Matter_O>();
//  cp->setContainedBy(ctemp);
  LOG("Adding %c:0x%08x to parent %c:0x%08x" , (cp->getMatterType()) , (cp.get()) , (this->getMatterType()) , (this ) );
    	// Always add the content to the end of the vector
	// A lot depends on Residues maintaining the order of Atoms
	// throughout the various passes of building databases
  this->_Contents.push_back(cp);
  LOG("Finished adding" );
}


// getId
int Matter_O::getId() const {
  return this->_Id;
}


//
// applyTransformToAtoms
//


CL_LISPIFY_NAME("applyTransformToAtoms");
CL_DEFMETHOD void	Matter_O::applyTransformToAtoms( const Matrix& m )
{_OF();
  contentIterator	a;
  LOG("Transforming all atoms with:%s" , (m.asString().c_str() ) );
  for ( a=this->begin_contents(); a!=this->end_contents(); a++ ) {
    (*a)->applyTransformToAtoms(m);
  }
    this->applyTransformToRestraints(m);
}

void Matter_O::applyTransformToRestraints(const Matrix& m)
{
  for ( auto cur : this->_RestraintList ) {
    Restraint_sp one = gc::As<Restraint_sp>(CONS_CAR(cur));
    one->applyTransform(m);
  }
}


CL_LISPIFY_NAME("invertStructureAndRestraints");
CL_DEFMETHOD void	Matter_O::invertStructureAndRestraints()
{
	// First invert all of the contained objects and their restraints
  contentIterator	a;
  for ( a=this->begin_contents(); a!=this->end_contents(); a++ )
  {
    (*a)->invertStructureAndRestraints();
  }
  this->invertStereochemistryOfRestraints();
}

void Matter_O::invertStereochemistryOfRestraints()
{_OF();
	//
	//Then flip any restraints
  for ( auto cur : this->_RestraintList ) {
    Restraint_sp one = gc::As<Restraint_sp>(CONS_CAR(cur));
    one->invertStereochemistryOfRestraint();
  }
}


CL_LISPIFY_NAME("translateAllAtoms");
CL_DEFMETHOD void Matter_O::translateAllAtoms(const Vector3& trans)
{
  Matrix transform;
  transform.translate(trans);
  this->applyTransformToAtoms(transform);
}



int	Matter_O::totalNetResidueCharge()
{
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



CL_LISPIFY_NAME("setAllAtomMasks");
CL_DEFMETHOD void	Matter_O::setAllAtomMasks(int m)
{
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
CL_LISPIFY_NAME("setAtomAliasesForResiduesNamed");
CL_DEFMETHOD void	Matter_O::setAtomAliasesForResiduesNamed(core::List_sp parts, core::List_sp atomAliases )
{
  FIX_ME();
#if 0
  if ( this->asSmartPtr().isA<Atom_O>() )
  {
    SIMPLE_ERROR(("Needs Aggregate, Molecule or Residue"));
  }
  LOG("setAtomAliasesForResiduesNamed with parts(%s) atomAliases(%s)" , _rep_(parts).c_str() , _rep_(atomAliases).c_str()  );
  gctools::SmallMap<MatterName,core::List_sp> residueNamesToAliasAtoms;
  for ( auto cur : parts ) {
    core::List_sp oneExtend = oCar(cur);
    if ( oneExtend->length() != 2 )
    {
      SIMPLE_ERROR(("Each entry must have two elements: "+_rep_(oneExtend) ));
    }
    core::String_sp residueName = oneExtend->car<core::String_O>();
    core::List_sp aliasAtoms = oCadr(oneExtend);
    LOG("residueName(%s) aliasAtoms(%s)" , residueName->get().c_str() , _rep_(aliasAtoms).c_str() );
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
      LOG("%s:%d Setting aliases for residue(%s)" , (__FILE__) , (__LINE__) , (res->getName().c_str()) );
      core::List_sp aliasAtoms = residueNamesToAliasAtoms[res->getName()];
      LOG("setting aliases for residueName(%s) aliasAtoms(%s) atomAliases(%s)" , res->getName().c_str() , _rep_(aliasAtoms).c_str() , _rep_(atomAliases).c_str() );
      res->setAliasesForAtoms(aliasAtoms,atomAliases);
    } else
    {
      LOG("%s:%d Not setting aliases for residue(%s)" , (__FILE__) , (__LINE__) , (res->getName().c_str()) );
    }
  }
#endif
}




CL_LISPIFY_NAME("geometricCenter");
CL_DEFMETHOD Vector3	Matter_O::geometricCenter()
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
  LOG("geometricCenterOfMatter number of atoms = %d" , count  );
  if ( count != 0 ) {
    sum = sum*(1.0/count);
  }
  LOG("Geometric center = %f, %f, %f" , sum.getX() , sum.getY() , sum.getZ()  );
  return sum;
}


CL_DOCSTRING(R"dx(Return the moment of inertia tensor where the masses of atoms are all identical. Return (values moment-of-inertia-matrix center-of-geometry))dx")
CL_LISPIFY_NAME("momentOfInertiaTensor");
CL_DEFMETHOD core::T_mv	Matter_O::momentOfInertiaTensor()
{_OF();
  Vector3	sum;
  int	count;
  Atom_sp	a;
  Loop	l;
  Matter_sp	c;
  double ixx = 0.0;
  double iyy = 0.0;
  double izz = 0.0;
  double ixy = 0.0;
  double iyz = 0.0;
  double ixz = 0.0;
  Vector3 center = this->geometricCenter();
  c = this->asSmartPtr();
  sum.set(0,0,0);
  count = 0;
  l.loopTopGoal( c, ATOMS );
  while ((l.advanceLoopAndProcess())) {
    a = l.getAtom();
    Vector3 pos = a->getPosition();
    double xpos = pos.getX()-center.getX();
    double ypos = pos.getY()-center.getY();
    double zpos = pos.getZ()-center.getZ();
    double xpossq = xpos*xpos;
    double ypossq = ypos*ypos;
    double zpossq = zpos*zpos;
    ixx += ypossq+zpossq;
    iyy += xpossq+zpossq;
    izz += xpossq+ypossq;
    ixy += - xpos*ypos;
    iyz += - ypos*zpos;
    ixz += - xpos*zpos;
    count++;
  }
  Matrix mres;
  mres.at(0, 0) = ixx;
  mres.at(0, 1) = ixy;
  mres.at(0, 2) = ixz;
  mres.at(1, 0) = ixy;
  mres.at(1, 1) = iyy;
  mres.at(1, 2) = iyz;
  mres.at(2, 0) = ixz;
  mres.at(2, 1) = iyz;
  mres.at(2, 2) = izz;
  geom::OMatrix_sp mm = geom::OMatrix_O::create();
  mm->setAll(mres);
  return Values(mm,geom::OVector3_O::createFromVector3(center));
}

CL_DOCSTRING(R"dx(Return the principle moments of inertia tensor where the masses of atoms are all identical. Return (values eigen-values eigen-vector-matrix center-of-geometry))dx")
CL_LISPIFY_NAME("momentOfGeometry");
CL_DEFMETHOD core::T_mv	Matter_O::momentOfGeometry()
{_OF();
  Vector3	sum;
  int	count;
  Atom_sp	a;
  Loop	l;
  Matter_sp	c;
  double ixx = 0.0;
  double iyy = 0.0;
  double izz = 0.0;
  double ixy = 0.0;
  double iyz = 0.0;
  double ixz = 0.0;
  Vector3 center = this->geometricCenter();
  c = this->asSmartPtr();
  sum.set(0,0,0);
  count = 0;
  l.loopTopGoal( c, ATOMS );
  while ((l.advanceLoopAndProcess())) {
    a = l.getAtom();
    Vector3 pos = a->getPosition();
    double xpos = pos.getX()-center.getX();
    double ypos = pos.getY()-center.getY();
    double zpos = pos.getZ()-center.getZ();
    double xpossq = xpos*xpos;
    double ypossq = ypos*ypos;
    double zpossq = zpos*zpos;
    ixx += ypossq+zpossq;
    iyy += xpossq+zpossq;
    izz += xpossq+ypossq;
    ixy += - xpos*ypos;
    iyz += - ypos*zpos;
    ixz += - xpos*zpos;
    count++;
  }
  Matrix mres;
  mres.at(0, 0) = ixx;
  mres.at(0, 1) = ixy;
  mres.at(0, 2) = ixz;
  mres.at(1, 0) = ixy;
  mres.at(1, 1) = iyy;
  mres.at(1, 2) = iyz;
  mres.at(2, 0) = ixz;
  mres.at(2, 1) = iyz;
  mres.at(2, 2) = izz;
  Vector4 evals;
  Matrix evecs;
  mres.eigenSystem(evals,evecs);
  core::List_sp lvals = core::Cons_O::createList(core::clasp_make_double_float(evals.getW()),
                                                 core::clasp_make_double_float(evals.getX()),
                                                 core::clasp_make_double_float(evals.getY()),
                                                 core::clasp_make_double_float(evals.getZ()));
  geom::OMatrix_sp mm = geom::OMatrix_O::create();
  mm->setAll(evecs);
  return Values(lvals,mm,geom::OVector3_O::createFromVector3(center));
}


CL_LISPIFY_NAME("matter-calculate-bounding-cuboid");
CL_DEFMETHOD geom::BoundingCuboid_sp Matter_O::matterCalculateBoundingCuboid(double pad)
{
  Vector3	sum;
  int	count;
  Atom_sp	a;
  Loop	l;
  Matter_sp	c;

  c = this->sharedThis<Matter_O>();
  sum.set(0,0,0);
  count = 0;
  geom::BoundingCuboid_sp bbox = geom::BoundingCuboid_O::create();
  l.loopTopGoal( c, ATOMS );
  while ((l.advanceLoopAndProcess())) {
    a = l.getAtom();
    bbox->expandToEncompassPoint(a->getPosition());
    count++;
  }
  bbox->pad(pad);
  return bbox;
}


CL_LISPIFY_NAME("contentsAsList");
CL_DEFMETHOD core::List_sp Matter_O::contentsAsList()
{
  ql::list result;
  contentIterator	it;
  for ( it=this->begin_contents(); it!=this->end_contents(); ++it )
    result << *it;
  return result.result();
}



CL_LISPIFY_NAME("allAtoms");
CL_DEFMETHOD core::Vector_sp	Matter_O::allAtoms()
{
  core::ComplexVector_T_sp	atoms = core::ComplexVector_T_O::make(8,nil<T_O>(),core::clasp_make_fixnum(0));
  Loop		la;
  la.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
  while ( la.advanceLoopAndProcess() )
  {
    atoms->vectorPushExtend(la.getAtom());
  }
  return atoms;
}


core::Vector_sp	Matter_O::allAtomsOfElement(Element element)
{
  core::ComplexVector_T_sp atoms = core::ComplexVector_T_O::make(8,nil<T_O>(),core::clasp_make_fixnum(0));
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


CL_LISPIFY_NAME("aliasResidueOrNil");
CL_DEFMETHOD Residue_sp	Matter_O::aliasResidueOrNil(Alias_sp alias)
{
  Loop lRes;
  lRes.loopTopGoal(this->sharedThis<Matter_O>(),RESIDUES);
  while ( lRes.advance() )
  {
    if ( lRes.getResidue()->recognizesMonomerAlias(alias->getMonomerAlias()) )
    {
      return lRes.getResidue();
    }
  }
  return nil<Residue_O>();
}

CL_LISPIFY_NAME("aliasResidue");
CL_DEFMETHOD Residue_sp Matter_O::aliasResidue(Alias_sp alias)
{_OF();
  Residue_sp a = this->aliasResidueOrNil(alias);
  if ( a.notnilp() ) return a;
  SIMPLE_ERROR(("could not find aliasResidue"));
}

CL_LISPIFY_NAME("aliasAtomOrNil");
CL_DEFMETHOD Atom_sp	Matter_O::aliasAtomOrNil(Alias_sp alias)
{
  Residue_sp res = this->aliasResidueOrNil(alias);
  if ( res.nilp() ) return nil<Atom_O>();
  return res->atomWithAliasOrNil(alias->getAtomAlias());
}

CL_LISPIFY_NAME("aliasAtom");
CL_DEFMETHOD Atom_sp Matter_O::aliasAtom(Alias_sp alias)
{_OF();
  Atom_sp a = this->aliasAtomOrNil(alias);
  if ( a.notnilp() ) return a;
  SIMPLE_ERROR(("could not find aliasAtom"));
}





CL_LISPIFY_NAME("allAtomsOfElementAsList");
CL_DEFMETHOD core::List_sp Matter_O::allAtomsOfElementAsList(Element element)
{
  core::List_sp		list(nil<core::T_O>());
  Loop		la;
  Atom_sp		a;
//  printf("%s:%d   At start list = |%s| list.nilp()--> %d  list@%p\n", __FILE__, __LINE__, _rep_(list).c_str(), list.nilp(), list.raw_());
  la.loopTopGoal(this->sharedThis<Matter_O>(),ATOMS);
  while ( la.advanceLoopAndProcess() ) {
    a = la.getAtom();
    if ( a->getElement() == element )
    {
      list = core::Cons_O::create(a,list);
//      printf("%s:%d  Into cons@%p  cdr@%p Pushed atom: %s\n", __FILE__, __LINE__, (void*)list.raw_(), (void*)oCdr(list).raw_(), _rep_(a).c_str());
    }
  }
  return list;
}

void	Matter_O::fields(core::Record_sp node )
{
  node->field( INTERN_(kw,name), this->_Name);
  node->/*pod_*/field_if_not_default( INTERN_(kw,id), this->_Id, 1);
  node->field_if_not_nil( INTERN_(kw,restraints),this->_RestraintList);
  node->field_if_not_nil( INTERN_(kw,properties),this->_Properties);
  node->field_if_not_empty( INTERN_(kw,contents), this->_Contents);
//  node->field_if_not_nil( INTERN_(kw,contained_by), this->containerContainedBy);
  LOG("Status" );
}






/*! Subclass this method
 */
#ifdef RENDER
geom::Render_sp Matter_O::rendered(core::List_sp kopts)
{
  return _Nil<geom::Render_O>();
};
#endif



CL_LISPIFY_NAME("allAtomsAsList");
CL_DEFMETHOD core::List_sp Matter_O::allAtomsAsList(bool allowVirtualAtoms ) const
{_OF();
  core::List_sp result = nil<core::T_O>();
  Loop l(this->const_sharedThis<Matter_O>(),ATOMS);
  while ( l.advance() )
  {
    Atom_sp a1 = l.getAtom();
    if ( !allowVirtualAtoms && a1->isVirtual() ) continue;
    result = core::Cons_O::create(a1,result);
  }
  return result;
}



CL_LISPIFY_NAME("allBondsAsList");
CL_DEFMETHOD core::List_sp Matter_O::allBondsAsList(bool allowVirtualAtoms ) const
{_OF();
  core::List_sp result = nil<core::List_V>();
  Loop l(this->const_sharedThis<Matter_O>(),BONDS);
  while ( l.advance() )
  {
    Bond_sp bond = l.getBond();
    LOG("Adding %s" , bond->description() );
    Atom_sp a1 = bond->getAtom1();
    Atom_sp a2 = bond->getAtom2();
    if ( !allowVirtualAtoms && ( a1->isVirtual() || a2->isVirtual() ) ) continue;
    result = core::Cons_O::create(bond,result);
  }
  return result;
}


core::HashTable_sp Matter_O::atomToResidueMap() {
  core::HashTableEq_sp map = core::HashTableEq_O::create_default();
  Loop lres(this->asSmartPtr(),RESIDUES);
  while (lres.advance()) {
    Residue_sp res = lres.getResidue();
    for (auto aa = res->begin_atoms(); aa!=res->end_atoms(); aa++ ) {
      map->setf_gethash(*aa,res);
    }
  }
  return map;
}


CL_LISPIFY_NAME("allAnglesAsList");
CL_DEFMETHOD core::List_sp Matter_O::allAnglesAsList(bool allowVirtualAtoms ) const
{_OF();
  core::List_sp result = nil<core::T_O>();
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


CL_LISPIFY_NAME("allImproperTorsionsAsList");
CL_DEFMETHOD core::List_sp Matter_O::allImproperTorsionsAsList(bool allowVirtualAtoms ) const
{_OF();
  core::List_sp result = nil<core::List_V>();
  Loop l(this->const_sharedThis<Matter_O>(),IMPROPERS);
  {
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
    LOG("There were %d improper torsions iterated over" , count );
  }
  return result;
}




CL_LISPIFY_NAME("allProperTorsionsAsList");
CL_DEFMETHOD core::List_sp Matter_O::allProperTorsionsAsList(bool allowVirtualAtoms ) const
{_OF();
  core::List_sp result = nil<core::T_O>();
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




Matter_sp Matter_O::copyDontRedirectAtoms(core::T_sp new_to_old)
{_OF();
  SUBCLASS_MUST_IMPLEMENT();
}

void Matter_O::redirectAtoms()
{_OF();
  SUBCLASS_MUST_IMPLEMENT();
}


void Matter_O::copyRestraintsDontRedirectAtoms(Matter_sp orig)
{_OF();
  this->_RestraintList = nil<core::T_O>();
  for ( auto cur : orig->_RestraintList ) {
    Restraint_sp restraint = gc::As<Restraint_sp>(CONS_CAR(cur));
    Restraint_sp copy_restraint = restraint->copyDontRedirectAtoms();
    this->_RestraintList = core::Cons_O::create(copy_restraint,this->_RestraintList);
  }
}


void Matter_O::redirectRestraintAtoms()
{
  for ( auto cur : this->_RestraintList ) {
    Restraint_sp restraint = gc::As<Restraint_sp>(CONS_CAR(cur));
    restraint->redirectAtoms();
  }
}


CL_LISPIFY_NAME("buildAtomIdMap");
CL_DEFMETHOD AtomIdMap_sp Matter_O::buildAtomIdMap() const
{
  SUBCLASS_MUST_IMPLEMENT();
}

CL_LISPIFY_NAME("atomWithAtomId");
CL_DEFMETHOD Atom_sp Matter_O::atomWithAtomId(const AtomId& atomId) const
{_OF();
  SUBCLASS_MUST_IMPLEMENT();
}

CL_LISPIFY_NAME("CHEM:MATTER-COPY");
CL_LAMBDA(original &optional new-to-old-map)
DOCGROUP(cando)
CL_DEFUN Matter_sp chem__matter_copy(Matter_sp orig, core::T_sp new_to_old) {
  Matter_sp copy = orig->copy(new_to_old);
  return copy;
}

};  // namespace chem
