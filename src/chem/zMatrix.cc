/*
    File: zMatrix.cc
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

#include <cando/chem/zMatrix.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/spanningLoop.h>
#include <cando/chem/loop.h>
#include <clasp/core/wrappers.h>

namespace chem
{












//
// Constructor
//
    void ZMatrixInternal_O::initialize()
    {
	this->Base::initialize();
	this->_Value = 0.0;
	this->_Constrain = false;
    }


#ifdef XML_ARCHIVE
    void	ZMatrixInternal_O::archiveBase(core::ArchiveP node)
    {
	node->archiveWeakPointer("WeakZMatrix",this->_WeakZMatrix);
	node->attribute("InternalName",this->_InternalName);
	node->attributeIfNotDefault("Constrain",this->_Constrain,false);
	node->attribute("AtomNew",this->_AtomNew);
	node->attributeIfNotDefault("Value",this->_Value,0.0);
    }
#endif




    Atom_sp	ZMatrixInternal_O::getAtomNew()
    {
	return this->getZMatrix()->_getAtomAtIndex(this->_AtomNew);
    };


    void	ZMatrixInternal_O::setAtomNew(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices)
    {
	this->_AtomNew = atomIndices[atom];
    }


    string ZMatrixInternal_O::__repr__() const
    {
	stringstream	ss;
	ss << this->_InternalName << " " << this->_AtomNew << " ";
	return ss.str();
    }










//
// Constructor
//

//
// Destructor
//


    string ZMatrixBondInternal_O::__repr__() const
    {
	stringstream	ss;
	ss << "Bond ";
	ss << this->ZMatrixInternal_O::__repr__();
	ss << " " << this->_AtomBond;
	return ss.str();
    }


    void	ZMatrixBondInternal_O::constrainIfNoAtomsWithMask(uint m)
    {
	Atom_sp	atomNew, atomBond;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	this->setConstrain(true);
	if ( atomNew->getMask()==m ||
	     atomBond->getMask()==m )
	{
	    this->setConstrain(false);
	}
    }

    void	ZMatrixBondInternal_O::extractInternal()
    {
	Atom_sp	atomNew, atomBond;
	double	value;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	value = geom::calculateDistance(atomNew->getPosition(),
				  atomBond->getPosition());
	this->setValue(value);
    }




CL_LISPIFY_NAME("getBondAtomZMatrixName");
CL_DEFMETHOD     core::T_sp	ZMatrixBondInternal_O::getBondAtomZMatrixName()
    {
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomBond);
    }


#ifdef XML_ARCHIVE
    void	ZMatrixBondInternal_O::archiveBase(core::ArchiveP node)
    {
	this->ZMatrixInternal_O::archiveBase(node);
	node->attribute("AtomBond",this->_AtomBond);
    }
#endif




    void	ZMatrixBondInternal_O::setup(Atom_sp atomNew, Atom_sp atomBond,
					     ZMatrixInternal_O::atomMap	atomIndices )
    {
	stringstream	name;
	this->setAtomNew(atomNew,atomIndices);
	name << "b" << atomIndices.size();
	this->setInternalName(name.str());
	this->_AtomBond = atomIndices[atomBond];
	this->setValue(0.0);
    }








//
// Constructor
//

//
// Destructor
//

    string ZMatrixAngleInternal_O::__repr__() const
    {
	stringstream	ss;
	ss << "Angle ";
	ss << this->ZMatrixInternal_O::__repr__();
	ss << " " << this->_AtomBond;
	ss << " " << this->_AtomAngle;
	return ss.str();
    }



    void	ZMatrixAngleInternal_O::constrainIfNoAtomsWithMask(uint m)
    {
	Atom_sp	atomNew, atomBond, atomAngle;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	atomAngle = this->getZMatrix()->_getAtomAtIndex(this->_AtomAngle);
	this->setConstrain(true);
	if ( atomNew->getMask()==m ||
	     atomBond->getMask()==m ||
	     atomAngle->getMask()==m )
	{
	    this->setConstrain(false);
	}
    }


    void	ZMatrixAngleInternal_O::extractInternal()
    {
	Atom_sp	atomNew, atomBond, atomAngle, atomDihedral;
	double	value;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	atomAngle = this->getZMatrix()->_getAtomAtIndex(this->_AtomAngle);
	value = geom::calculateAngle(atomNew->getPosition(),
			       atomBond->getPosition(),
			       atomAngle->getPosition() );
	this->setValue(value);
    }




#ifdef XML_ARCHIVE
    void	ZMatrixAngleInternal_O::archiveBase(core::ArchiveP node)
    {
	this->ZMatrixInternal_O::archiveBase(node);
	node->attribute("AtomBond",this->_AtomBond);
	node->attribute("AtomAngle",this->_AtomAngle);
    }
#endif

CL_LISPIFY_NAME("getBondAtomZMatrixName");
CL_DEFMETHOD     core::T_sp	ZMatrixAngleInternal_O::getBondAtomZMatrixName()
    {
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomBond);
    }

CL_LISPIFY_NAME("getAngleAtomZMatrixName");
CL_DEFMETHOD     core::T_sp	ZMatrixAngleInternal_O::getAngleAtomZMatrixName()
    {
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomAngle);
    }

    void	ZMatrixAngleInternal_O::setup(Atom_sp atomNew, Atom_sp atomBond, Atom_sp atomAngle,
					      ZMatrixInternal_O::atomMap	atomIndices )
    {
	stringstream	name;
	this->setAtomNew(atomNew,atomIndices);
	name << "a" << atomIndices.size();
	this->setInternalName(name.str());
	this->_AtomBond = atomIndices[atomBond];
	this->_AtomAngle = atomIndices[atomAngle];
	this->setValue(0.0);
    }










//
// Constructor
//

//
// Destructor
//

    string ZMatrixDihedralInternal_O::__repr__() const
    {
	stringstream	ss;
	ss << "Dihedral ";
	ss << this->ZMatrixInternal_O::__repr__();
	ss << " " << this->_AtomBond;
	ss << " " << this->_AtomAngle;
	ss << " " << this->_AtomDihedral;
	return ss.str();
    }

    void	ZMatrixDihedralInternal_O::constrainIfNoAtomsWithMask(uint m)
    {
	Atom_sp	atomNew, atomBond, atomAngle, atomDihedral;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	atomAngle = this->getZMatrix()->_getAtomAtIndex(this->_AtomAngle);
	atomDihedral = this->getZMatrix()->_getAtomAtIndex(this->_AtomDihedral);
	this->setConstrain(true);
	if ( atomNew->getMask()==m ||
	     atomBond->getMask()==m ||
	     atomAngle->getMask()==m ||
	     atomDihedral->getMask()==m )
	{
	    this->setConstrain(false);
	}
    }

    void	ZMatrixDihedralInternal_O::extractInternal()
    {
	Atom_sp	atomNew, atomBond, atomAngle, atomDihedral;
	double	value;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	atomAngle = this->getZMatrix()->_getAtomAtIndex(this->_AtomAngle);
	atomDihedral = this->getZMatrix()->_getAtomAtIndex(this->_AtomDihedral);
	value = geom::calculateDihedral(atomNew->getPosition(),
				  atomBond->getPosition(),
				  atomAngle->getPosition(),
				  atomDihedral->getPosition());
	this->setValue(value);
    }


#ifdef XML_ARCHIVE
    void	ZMatrixDihedralInternal_O::archiveBase(core::ArchiveP node)
    {
	this->ZMatrixInternal_O::archiveBase(node);
	node->attribute("AtomBond",this->_AtomBond);
	node->attribute("AtomAngle",this->_AtomAngle);
	node->attribute("AtomDihedral",this->_AtomDihedral);
    }
#endif

CL_LISPIFY_NAME("getBondAtomZMatrixName");
CL_DEFMETHOD     core::T_sp	ZMatrixDihedralInternal_O::getBondAtomZMatrixName()
    {
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomBond);
    }

CL_LISPIFY_NAME("getAngleAtomZMatrixName");
CL_DEFMETHOD     core::T_sp	ZMatrixDihedralInternal_O::getAngleAtomZMatrixName()
    {
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomAngle);
    }

CL_LISPIFY_NAME("getDihedralAtomZMatrixName");
CL_DEFMETHOD     core::T_sp	ZMatrixDihedralInternal_O::getDihedralAtomZMatrixName()
    {
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomDihedral);
    }




    void	ZMatrixDihedralInternal_O::setup(Atom_sp atomNew, Atom_sp atomBond,
						 Atom_sp atomAngle, Atom_sp atomDihedral,
						 ZMatrixInternal_O::atomMap	atomIndices )
    {
	stringstream	name;
	this->setAtomNew(atomNew,atomIndices);
	name << "d" << atomIndices.size();
	this->setInternalName(name.str());
	this->_AtomBond = atomIndices[atomBond];
	this->_AtomAngle = atomIndices[atomAngle];
	this->_AtomDihedral = atomIndices[atomDihedral];
	this->setValue(0.0);
    }









//
// Constructor
//


    ZMatrixEntry_sp ZMatrixEntry_O::create(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices)
    {
      auto 	entry  = gctools::GC<ZMatrixEntry_O>::allocate_with_default_constructor();
	stringstream name;
	name << atom->getElementAsString() << atomIndices.size();
	entry->_Atom = atom;
	entry->_ZMatrixAtomName = core::SimpleBaseString_O::make(name.str());
	return entry;
    }


    void	ZMatrixEntry_O::initialize()
    {
	this->Base::initialize();
	this->_ZMatrixAtomName = unbound<core::T_O>();
	this->_Atom = nil<Atom_O>();
	this->_Bond = nil<ZMatrixBondInternal_O>();
	this->_Angle = nil<ZMatrixAngleInternal_O>();
	this->_Dihedral = nil<ZMatrixDihedralInternal_O>();
    }

//
// Destructor
//

#ifdef XML_ARCHIVE
    void	ZMatrixEntry_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("ZMatrixAtomName",this->_ZMatrixAtomName);
	ASSERTNOTNULL(this->_Atom);
	node->attribute("Atom",this->_Atom);
	node->attributeIfNotNil("Bond",this->_Bond);
	node->attributeIfNotNil("Angle",this->_Angle);
	node->attributeIfNotNil("Dihedral",this->_Dihedral);
    }
#endif


CL_LISPIFY_NAME("get-zmatrix-atom-name-at-index");
CL_DEFMETHOD     core::T_sp ZMatrix_O::getZMatrixAtomNameAtIndex(uint i) const
    {
	return this->_getAtomZMatrixNameAtIndex(i);
    }








//
// Constructor
//
#if INIT_TO_FACTORIES

#define ARGS_ZMatrix_O_make "(matter root_atom)"
#define DECL_ZMatrix_O_make ""
#define DOCS_ZMatrix_O_make "make ZMatrix"
  ZMatrix_sp ZMatrix_O::make(Matter_sp matter, Atom_sp rootAtom)
  {
    auto  me  = gctools::GC<ZMatrix_O>::allocate_with_default_constructor();
      if ( rootAtom.nilp() )
      {
	  me->defineForMatter(matter);
      } else
      {
	  me->defineForMatterWithStartingAtom(matter,rootAtom);
      }
      return me;
  };

#else

    core::T_sp ZMatrix_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp)
    {
	Matter_sp matter = translate::from_object<Matter_sp>::convert(env->lookup(Pkg(),"matter"));
	Atom_sp rootAtom = translate::from_object<Atom_sp>::convert(env->lookup(Pkg(),"rootAtom"));
	if ( rootAtom.nilp() )
	{
	    this->defineForMatter(matter);
	} else
	{
	    this->defineForMatterWithStartingAtom(matter,rootAtom);
	}
	return _Nil<core::T_O>();
    }

#endif

//
// Destructor
//

#ifdef XML_ARCHIVE
    void	ZMatrix_O::archiveBase(core::ArchiveP node)
    {
	node->attribute("Matter",this->_Matter);
	node->archiveVector0("ZMatrix",this->_ZMatrix);
	node->archiveVector0("Internals",this->_Internals);
    }
#endif

    Atom_sp	ZMatrix_O::_getAtomAtIndex(uint i)
    {
	ASSERT_lessThan(i,this->_ZMatrix.size());
	return this->_ZMatrix[i]->getAtom();
    }

core::T_sp	ZMatrix_O::_getAtomZMatrixNameAtIndex(uint i) const
{
  LOG("Looking for atom at index: %d    _ZMatrix.size=%d" , i , this->_ZMatrix.size()  );
  ASSERT_lessThan(i,this->_ZMatrix.size());
  return this->_ZMatrix[i]->getTargetAtomZMatrixName();
}

    void	ZMatrix_O::defineForMatterWithStartingAtom(Matter_sp matter, Atom_sp atom)
    {
	Atom_sp				newAtom, bondToAtom, angleToAtom, dihedralToAtom;
	ZMatrixEntry_sp			entry;
	ZMatrixBondInternal_sp		bondInternal;
	ZMatrixAngleInternal_sp		angleInternal;
	ZMatrixDihedralInternal_sp	dihedralInternal;
	SpanningLoop_sp			span;
	ZMatrixInternal_O::atomMap 	atomIndices;
	uint				numberOfEntries;
	LOG("Starting on atom: %s" , atom->description().c_str()  );
	this->_Matter = matter;
	span = SpanningLoop_O::create();
	span->setTop(atom);
	// Build the spanning tree
	while ( span->advanceLoopAndProcess() )
	{
	    newAtom = span->getAtom();
	    
	    atomIndices[newAtom] = this->_ZMatrix.size();
	    entry = ZMatrixEntry_O::create(newAtom,atomIndices);
	    this->_ZMatrix.push_back(entry);
	    numberOfEntries = this->_ZMatrix.size();
	    if ( !span->isBackSpanValid(newAtom) )
	    {
		continue;
	    }
	    bondToAtom = gc::As_unsafe<Atom_sp>(span->getBackSpan(newAtom));
	    LOG("bondToAtom = %s" , bondToAtom->description().c_str() );
	    bondInternal = ZMatrixBondInternal_O::create(newAtom,bondToAtom,atomIndices,this->sharedThis<ZMatrix_O>());
	    this->_Internals.push_back(bondInternal);
	    entry->_Bond = bondInternal;
	    if ( !span->isBackSpanValid(bondToAtom) )
	    {
		LOG("Stopped at bond" );
		continue;
	    }
	    angleToAtom = gc::As_unsafe<Atom_sp>(span->getBackSpan(bondToAtom));
	    LOG("angleToAtom = %s" , angleToAtom->description().c_str() );
	    angleInternal = ZMatrixAngleInternal_O::create(newAtom,bondToAtom,
							   angleToAtom,atomIndices,this->sharedThis<ZMatrix_O>());
	    this->_Internals.push_back(angleInternal);
	    entry->_Angle = angleInternal;
            core::T_sp backspan = span->getBackSpan(angleToAtom);
	    if ( backspan.nilp() ) //!span->isBackSpanValid(angleToAtom) )
	    {
		LOG("Stopped at angle" );
		if ( numberOfEntries == 3 ) continue;
		if ( numberOfEntries == 4 || numberOfEntries == 5 )
		{
		    dihedralToAtom = this->_getAtomAtIndex(2);
		} else
		{
		    LOG("Stopped at angle but numberOfEntries=%d" , numberOfEntries  );
		    SIMPLE_ERROR(("There are too many entries to run out of a dihedral atom"));
		}
	    } else
	    {
		dihedralToAtom = gc::As_unsafe<Atom_sp>(span->getBackSpan(angleToAtom));
	    }
	    LOG("dihedralToAtom = %s" , dihedralToAtom->description().c_str() );
	    dihedralInternal = ZMatrixDihedralInternal_O::create(newAtom,bondToAtom,
								 angleToAtom,dihedralToAtom,atomIndices,this->sharedThis<ZMatrix_O>());
	    LOG("Appending dihedral" );
	    this->_Internals.push_back(dihedralInternal);
	    LOG("Setting entry->_Dihedral" );
	    entry->_Dihedral = dihedralInternal;
	}
	LOG("Finished building zMatrix" );
    }



    void	ZMatrix_O::defineForMatter(Matter_sp matter)
    {
	Atom_sp				atom;
	Loop				lAtoms;
	bool				gotPeripheralAtom;
	gotPeripheralAtom = false;
	lAtoms.loopTopGoal(matter,ATOMS);
	while ( lAtoms.advanceLoopAndProcess() )
	{
	    atom = lAtoms.getAtom();
	    if ( atom->numberOfBonds() == 1 )
	    {
		gotPeripheralAtom = true;
		break;
	    }
	}
	ASSERTP(gotPeripheralAtom,"There were no peripheral atoms in the matter");
	this->defineForMatterWithStartingAtom(matter,atom);
    }


CL_LISPIFY_NAME("extractInternals");
CL_DEFMETHOD     void	ZMatrix_O::extractInternals()
    {
	ZMatrix_O::internalIterator	it;
	for ( it=this->_Internals.begin(); it!=this->_Internals.end(); it++ )
	{
	    (*it)->extractInternal();
	}
    }


    void	ZMatrix_O::constrainInternalsWithNoAtomsWithMask(uint mask)
    {
	ZMatrix_O::internalIterator	it;
	for ( it=this->_Internals.begin(); it!=this->_Internals.end(); it++ )
	{
	    (*it)->constrainIfNoAtomsWithMask(mask);
	}

    }



CL_LISPIFY_NAME("entriesAsList");
CL_DEFMETHOD     core::List_sp ZMatrix_O::entriesAsList() const
    {
	return core::Cons_O::createFromVec0(this->_ZMatrix);
    };

CL_LISPIFY_NAME("zmatrix-entry-with-name");
CL_DEFMETHOD     ZMatrixEntry_sp ZMatrix_O::zMatrixEntryWithName(core::T_sp atomZMatrixName) const
    {
	for ( const_zMatrixEntryIterator it = this->begin_ZMatrixEntries();
	      it != this->end_ZMatrixEntries(); it++ )
	{
	    if ( (*it)->getTargetAtomZMatrixName()==atomZMatrixName) return *it;
	}
	return nil<ZMatrixEntry_O>();
    }


CL_LISPIFY_NAME("all-zmatrix-atom-names-as-string");
CL_DEFMETHOD     string ZMatrix_O::allZMatrixAtomNamesAsString() const
    {
	stringstream ss;
	for ( const_zMatrixEntryIterator it = this->begin_ZMatrixEntries();
	      it != this->end_ZMatrixEntries(); it++ )
	{
	    ss << (*it)->getTargetAtomZMatrixName() << " ";
	}
	return ss.str();
    }




};
