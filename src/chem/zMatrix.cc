         
#define	DEBUG_LEVEL_FULL

#include <cando/chem/zMatrix.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/spanningLoop.h>
#include <cando/chem/loop.h>
#include <clasp/core/wrappers.h>

namespace chem
{



    EXPOSE_CLASS(chem,ZMatrixInternal_O);
    EXPOSE_CLASS(chem,ZMatrixBondInternal_O);
    EXPOSE_CLASS(chem,ZMatrixAngleInternal_O);
    EXPOSE_CLASS(chem,ZMatrixDihedralInternal_O);
    EXPOSE_CLASS(chem,ZMatrixEntry_O);
    EXPOSE_CLASS(chem,ZMatrix_O);


    void ZMatrixInternal_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ZMatrixInternal_O>()
	    .def("getValue",&ZMatrixInternal_O::getValue)
	    ;
    }

    void ZMatrixInternal_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ZMatrixInternal,"","",_lisp)
	    .def("getValue",&ZMatrixInternal_O::getValue)
	    ;
#endif
    }






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
    {_G();
	node->archiveWeakPointer("WeakZMatrix",this->_WeakZMatrix);
	node->attribute("InternalName",this->_InternalName);
	node->attributeIfNotDefault("Constrain",this->_Constrain,false);
	node->attribute("AtomNew",this->_AtomNew);
	node->attributeIfNotDefault("Value",this->_Value,0.0);
    }
#endif




    Atom_sp	ZMatrixInternal_O::getAtomNew()
    {_G();
	return this->getZMatrix()->_getAtomAtIndex(this->_AtomNew);
    };


    void	ZMatrixInternal_O::setAtomNew(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices)
    {_G();
	this->_AtomNew = atomIndices[atom];
    }


    string ZMatrixInternal_O::__repr__() const
    {_G();
	stringstream	ss;
	ss << this->_InternalName << " " << this->_AtomNew << " ";
	return ss.str();
    }



    void ZMatrixBondInternal_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ZMatrixBondInternal_O>()
	    .def("getBondAtomZMatrixName",&ZMatrixBondInternal_O::getBondAtomZMatrixName)
	    .def("getBondAtomIndex",&ZMatrixBondInternal_O::getBondAtomIndex)
	    ;
    }

    void ZMatrixBondInternal_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ZMatrixBondInternal,"","",_lisp)
	    .def("getBondAtomZMatrixName",&ZMatrixBondInternal_O::getBondAtomZMatrixName)
	    .def("getBondAtomIndex",&ZMatrixBondInternal_O::getBondAtomIndex)
	    ;
#endif
    }






//
// Constructor
//

//
// Destructor
//


    string ZMatrixBondInternal_O::__repr__() const
    {_G();
	stringstream	ss;
	ss << "Bond ";
	ss << this->ZMatrixInternal_O::__repr__();
	ss << " " << this->_AtomBond;
	return ss.str();
    }


    void	ZMatrixBondInternal_O::constrainIfNoAtomsWithMask(uint m)
    {_G();
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
    {_G();
	Atom_sp	atomNew, atomBond;
	double	value;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	value = calculateDistance(atomNew->getPosition(),
				  atomBond->getPosition());
	this->setValue(value);
    }




    string	ZMatrixBondInternal_O::getBondAtomZMatrixName()
    {_G();
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomBond);
    }


#ifdef XML_ARCHIVE
    void	ZMatrixBondInternal_O::archiveBase(core::ArchiveP node)
    {_G();
	this->ZMatrixInternal_O::archiveBase(node);
	node->attribute("AtomBond",this->_AtomBond);
    }
#endif




    void	ZMatrixBondInternal_O::setup(Atom_sp atomNew, Atom_sp atomBond,
					     ZMatrixInternal_O::atomMap	atomIndices )
    {_G();
	stringstream	name;
	this->setAtomNew(atomNew,atomIndices);
	name << "b" << atomIndices.size();
	this->setInternalName(name.str());
	this->_AtomBond = atomIndices[atomBond];
	this->setValue(0.0);
    }




    void ZMatrixAngleInternal_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ZMatrixAngleInternal_O>()
	    .def("getBondAtomZMatrixName",&ZMatrixAngleInternal_O::getBondAtomZMatrixName)
	    .def("getAngleAtomZMatrixName",&ZMatrixAngleInternal_O::getAngleAtomZMatrixName)
	    .def("getBondAtomIndex",&ZMatrixAngleInternal_O::getBondAtomIndex)
	    .def("getAngleAtomIndex",&ZMatrixAngleInternal_O::getAngleAtomIndex)
	    ;
    }

    void ZMatrixAngleInternal_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ZMatrixAngleInternal,"","",_lisp)
	    .def("getBondAtomZMatrixName",&ZMatrixAngleInternal_O::getBondAtomZMatrixName)
	    .def("getAngleAtomZMatrixName",&ZMatrixAngleInternal_O::getAngleAtomZMatrixName)
	    .def("getBondAtomIndex",&ZMatrixAngleInternal_O::getBondAtomIndex)
	    .def("getAngleAtomIndex",&ZMatrixAngleInternal_O::getAngleAtomIndex)
	    ;
#endif
    }



//
// Constructor
//

//
// Destructor
//

    string ZMatrixAngleInternal_O::__repr__() const
    {_G();
	stringstream	ss;
	ss << "Angle ";
	ss << this->ZMatrixInternal_O::__repr__();
	ss << " " << this->_AtomBond;
	ss << " " << this->_AtomAngle;
	return ss.str();
    }



    void	ZMatrixAngleInternal_O::constrainIfNoAtomsWithMask(uint m)
    {_G();
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
    {_G();
	Atom_sp	atomNew, atomBond, atomAngle, atomDihedral;
	double	value;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	atomAngle = this->getZMatrix()->_getAtomAtIndex(this->_AtomAngle);
	value = calculateAngle(atomNew->getPosition(),
			       atomBond->getPosition(),
			       atomAngle->getPosition() );
	this->setValue(value);
    }




#ifdef XML_ARCHIVE
    void	ZMatrixAngleInternal_O::archiveBase(core::ArchiveP node)
    {_G();
	this->ZMatrixInternal_O::archiveBase(node);
	node->attribute("AtomBond",this->_AtomBond);
	node->attribute("AtomAngle",this->_AtomAngle);
    }
#endif

    string	ZMatrixAngleInternal_O::getBondAtomZMatrixName()
    {_G();
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomBond);
    }

    string	ZMatrixAngleInternal_O::getAngleAtomZMatrixName()
    {_G();
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomAngle);
    }

    void	ZMatrixAngleInternal_O::setup(Atom_sp atomNew, Atom_sp atomBond, Atom_sp atomAngle,
					      ZMatrixInternal_O::atomMap	atomIndices )
    {_G();
	stringstream	name;
	this->setAtomNew(atomNew,atomIndices);
	name << "a" << atomIndices.size();
	this->setInternalName(name.str());
	this->_AtomBond = atomIndices[atomBond];
	this->_AtomAngle = atomIndices[atomAngle];
	this->setValue(0.0);
    }



    void ZMatrixDihedralInternal_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ZMatrixDihedralInternal_O>()
	    .def("getBondAtomZMatrixName",&ZMatrixDihedralInternal_O::getBondAtomZMatrixName)
	    .def("getAngleAtomZMatrixName",&ZMatrixDihedralInternal_O::getAngleAtomZMatrixName)
	    .def("getDihedralAtomZMatrixName",&ZMatrixDihedralInternal_O::getDihedralAtomZMatrixName)
	    .def("getBondAtomIndex",&ZMatrixDihedralInternal_O::getBondAtomIndex)
	    .def("getAngleAtomIndex",&ZMatrixDihedralInternal_O::getAngleAtomIndex)
	    .def("getDihedralAtomIndex",&ZMatrixDihedralInternal_O::getDihedralAtomIndex)
	    ;
    }

    void ZMatrixDihedralInternal_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ZMatrixDihedralInternal,"","",_lisp)
	    .def("getBondAtomZMatrixName",&ZMatrixDihedralInternal_O::getBondAtomZMatrixName)
	    .def("getAngleAtomZMatrixName",&ZMatrixDihedralInternal_O::getAngleAtomZMatrixName)
	    .def("getDihedralAtomZMatrixName",&ZMatrixDihedralInternal_O::getDihedralAtomZMatrixName)
	    .def("getBondAtomIndex",&ZMatrixDihedralInternal_O::getBondAtomIndex)
	    .def("getAngleAtomIndex",&ZMatrixDihedralInternal_O::getAngleAtomIndex)
	    .def("getDihedralAtomIndex",&ZMatrixDihedralInternal_O::getDihedralAtomIndex)
	    ;
#endif
    }






//
// Constructor
//

//
// Destructor
//

    string ZMatrixDihedralInternal_O::__repr__() const
    {_G();
	stringstream	ss;
	ss << "Dihedral ";
	ss << this->ZMatrixInternal_O::__repr__();
	ss << " " << this->_AtomBond;
	ss << " " << this->_AtomAngle;
	ss << " " << this->_AtomDihedral;
	return ss.str();
    }

    void	ZMatrixDihedralInternal_O::constrainIfNoAtomsWithMask(uint m)
    {_G();
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
    {_G();
	Atom_sp	atomNew, atomBond, atomAngle, atomDihedral;
	double	value;
	atomNew = this->getAtomNew();
	atomBond = this->getZMatrix()->_getAtomAtIndex(this->_AtomBond);
	atomAngle = this->getZMatrix()->_getAtomAtIndex(this->_AtomAngle);
	atomDihedral = this->getZMatrix()->_getAtomAtIndex(this->_AtomDihedral);
	value = calculateDihedral(atomNew->getPosition(),
				  atomBond->getPosition(),
				  atomAngle->getPosition(),
				  atomDihedral->getPosition());
	this->setValue(value);
    }


#ifdef XML_ARCHIVE
    void	ZMatrixDihedralInternal_O::archiveBase(core::ArchiveP node)
    {_G();
	this->ZMatrixInternal_O::archiveBase(node);
	node->attribute("AtomBond",this->_AtomBond);
	node->attribute("AtomAngle",this->_AtomAngle);
	node->attribute("AtomDihedral",this->_AtomDihedral);
    }
#endif

    string	ZMatrixDihedralInternal_O::getBondAtomZMatrixName()
    {_G();
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomBond);
    }

    string	ZMatrixDihedralInternal_O::getAngleAtomZMatrixName()
    {_G();
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomAngle);
    }

    string	ZMatrixDihedralInternal_O::getDihedralAtomZMatrixName()
    {_G();
	return this->getZMatrix()->_getAtomZMatrixNameAtIndex(this->_AtomDihedral);
    }




    void	ZMatrixDihedralInternal_O::setup(Atom_sp atomNew, Atom_sp atomBond,
						 Atom_sp atomAngle, Atom_sp atomDihedral,
						 ZMatrixInternal_O::atomMap	atomIndices )
    {_G();
	stringstream	name;
	this->setAtomNew(atomNew,atomIndices);
	name << "d" << atomIndices.size();
	this->setInternalName(name.str());
	this->_AtomBond = atomIndices[atomBond];
	this->_AtomAngle = atomIndices[atomAngle];
	this->_AtomDihedral = atomIndices[atomDihedral];
	this->setValue(0.0);
    }




    void ZMatrixEntry_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ZMatrixEntry_O>()
	    .def("getTargetAtomZMatrixName",&ZMatrixEntry_O::getTargetAtomZMatrixName)
	    .def("getAtom",&ZMatrixEntry_O::getAtom)
	    .def("hasBond",&ZMatrixEntry_O::hasBond)
	    .def("getBond",&ZMatrixEntry_O::getBond)
	    .def("hasAngle",&ZMatrixEntry_O::hasAngle)
	    .def("getAngle",&ZMatrixEntry_O::getAngle)
	    .def("hasDihedral",&ZMatrixEntry_O::hasDihedral)
	    .def("getDihedral",&ZMatrixEntry_O::getDihedral)
	    ;
    }

    void ZMatrixEntry_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ZMatrixEntry,"","",_lisp)
	    .def("getTargetAtomZMatrixName",&ZMatrixEntry_O::getTargetAtomZMatrixName)
	    .def("getAtom",&ZMatrixEntry_O::getAtom)
	    .def("hasBond",&ZMatrixEntry_O::hasBond)
	    .def("getBond",&ZMatrixEntry_O::getBond)
	    .def("hasAngle",&ZMatrixEntry_O::hasAngle)
	    .def("getAngle",&ZMatrixEntry_O::getAngle)
	    .def("hasDihedral",&ZMatrixEntry_O::hasDihedral)
	    .def("getDihedral",&ZMatrixEntry_O::getDihedral)
	    ;
#endif
    }




//
// Constructor
//


    ZMatrixEntry_sp ZMatrixEntry_O::create(Atom_sp atom, ZMatrixInternal_O::atomMap atomIndices)
    {_G();
	GC_ALLOCATE(ZMatrixEntry_O,	entry );
	stringstream name;
	name << atom->getElementAsString() << atomIndices.size();
	entry->_Atom = atom;
	entry->_ZMatrixAtomName = name.str();
	return entry;
    }


    void	ZMatrixEntry_O::initialize()
    {
	this->Base::initialize();
	this->_ZMatrixAtomName = "";
	this->_Atom = _Nil<Atom_O>();
	this->_Bond = _Nil<ZMatrixBondInternal_O>();
	this->_Angle = _Nil<ZMatrixAngleInternal_O>();
	this->_Dihedral = _Nil<ZMatrixDihedralInternal_O>();
    }

//
// Destructor
//

#ifdef XML_ARCHIVE
    void	ZMatrixEntry_O::archiveBase(core::ArchiveP node)
    {_G();
	node->attribute("ZMatrixAtomName",this->_ZMatrixAtomName);
	ASSERTNOTNULL(this->_Atom);
	node->attribute("Atom",this->_Atom);
	node->attributeIfNotNil("Bond",this->_Bond);
	node->attributeIfNotNil("Angle",this->_Angle);
	node->attributeIfNotNil("Dihedral",this->_Dihedral);
    }
#endif


    string ZMatrix_O::getZMatrixAtomNameAtIndex(uint i) const
    {_OF();
	return this->_getAtomZMatrixNameAtIndex(i);
    }


    void ZMatrix_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<ZMatrix_O>()
//	    .def_raw("core:__init__",&ZMatrix_O::__init__,"(self matter &key rootAtom)")
	    .def("entriesAsCons",&ZMatrix_O::entriesAsCons)
	    .def("extractInternals",&ZMatrix_O::extractInternals)
	    .def("zmatrix-entry-with-name",&ZMatrix_O::zMatrixEntryWithName)
	    .def("all-zmatrix-atom-names-as-string",&ZMatrix_O::allZMatrixAtomNamesAsString)
	    .def("get-zmatrix-atom-name-at-index",&ZMatrix_O::getZMatrixAtomNameAtIndex)
	    ;
    }

    void ZMatrix_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,ZMatrix,"","",_lisp)
	    .def("entriesAsCons",&ZMatrix_O::entriesAsCons)
	    .def("extractInternals",&ZMatrix_O::extractInternals)
	    .def("zmatrix-entry-with-name",&ZMatrix_O::zMatrixEntryWithName)
	    .def("all-zmatrix-atom-names-as-string",&ZMatrix_O::allZMatrixAtomNamesAsString)
	    .def("get-zmatrix-atom-name-at-index",&ZMatrix_O::getZMatrixAtomNameAtIndex)
	    ;
#endif
    }





//
// Constructor
//
#if INIT_TO_FACTORIES

#define ARGS_ZMatrix_O_make "(matter root_atom)"
#define DECL_ZMatrix_O_make ""
#define DOCS_ZMatrix_O_make "make ZMatrix"
  ZMatrix_sp ZMatrix_O::make(Matter_sp matter, Atom_sp rootAtom)
  {_G();
      GC_ALLOCATE(ZMatrix_O, me );
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

    core::T_sp ZMatrix_O::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
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
    {_G();
	node->attribute("Matter",this->_Matter);
	node->archiveVector0("ZMatrix",this->_ZMatrix);
	node->archiveVector0("Internals",this->_Internals);
    }
#endif

    Atom_sp	ZMatrix_O::_getAtomAtIndex(uint i)
    {_G();
	ASSERT_lessThan(i,this->_ZMatrix.size());
	return this->_ZMatrix[i]->getAtom();
    }

    string	ZMatrix_O::_getAtomZMatrixNameAtIndex(uint i) const
    {_G();
	LOG(BF("Looking for atom at index: %d    _ZMatrix.size=%d") % i % this->_ZMatrix.size()  );
	ASSERT_lessThan(i,this->_ZMatrix.size());
	return this->_ZMatrix[i]->getTargetAtomZMatrixName();
    }

    void	ZMatrix_O::defineForMatterWithStartingAtom(Matter_sp matter, Atom_sp atom)
    {_G();
	Atom_sp				newAtom, bondToAtom, angleToAtom, dihedralToAtom;
	ZMatrixEntry_sp			entry;
	ZMatrixBondInternal_sp		bondInternal;
	ZMatrixAngleInternal_sp		angleInternal;
	ZMatrixDihedralInternal_sp	dihedralInternal;
	SpanningLoop_sp			span;
	ZMatrixInternal_O::atomMap 	atomIndices;
	uint				numberOfEntries;
	LOG(BF("Starting on atom: %s") % atom->description().c_str()  );
	this->_Matter = matter;
	span = SpanningLoop_O::create();
	span->setTop(atom);
	// Build the spanning tree
	while ( span->advanceLoopAndProcess() )
	{
	    newAtom = span->getAtom();
	    _BLOCK_TRACEF(BF("ZMatrix entry for: %s") % newAtom->description().c_str() );
	    atomIndices[newAtom] = this->_ZMatrix.size();
	    entry = ZMatrixEntry_O::create(newAtom,atomIndices);
	    this->_ZMatrix.push_back(entry);
	    numberOfEntries = this->_ZMatrix.size();
	    if ( ! newAtom->isBackSpanValid() )
	    {
		continue;
	    }
	    bondToAtom = newAtom->getBackSpan();
	    LOG(BF("bondToAtom = %s") % bondToAtom->description().c_str() );
	    bondInternal = ZMatrixBondInternal_O::create(newAtom,bondToAtom,atomIndices,this->sharedThis<ZMatrix_O>());
	    this->_Internals.push_back(bondInternal);
	    entry->_Bond = bondInternal;
	    if ( !bondToAtom->isBackSpanValid() )
	    {
		LOG(BF("Stopped at bond") );
		continue;
	    }
	    angleToAtom = bondToAtom->getBackSpan();
	    LOG(BF("angleToAtom = %s") % angleToAtom->description().c_str() );
	    angleInternal = ZMatrixAngleInternal_O::create(newAtom,bondToAtom,
							   angleToAtom,atomIndices,this->sharedThis<ZMatrix_O>());
	    this->_Internals.push_back(angleInternal);
	    entry->_Angle = angleInternal;
	    if ( !angleToAtom->isBackSpanValid() )
	    {
		LOG(BF("Stopped at angle") );
		if ( numberOfEntries == 3 ) continue;
		if ( numberOfEntries == 4 || numberOfEntries == 5 )
		{
		    dihedralToAtom = this->_getAtomAtIndex(2);
		} else
		{
		    LOG(BF("Stopped at angle but numberOfEntries=%d") % numberOfEntries  );
		    SIMPLE_ERROR(BF("There are too many entries to run out of a dihedral atom"));
		}
	    } else
	    {
		dihedralToAtom = angleToAtom->getBackSpan();
	    }
	    LOG(BF("dihedralToAtom = %s") % dihedralToAtom->description().c_str() );
	    dihedralInternal = ZMatrixDihedralInternal_O::create(newAtom,bondToAtom,
								 angleToAtom,dihedralToAtom,atomIndices,this->sharedThis<ZMatrix_O>());
	    LOG(BF("Appending dihedral") );
	    this->_Internals.push_back(dihedralInternal);
	    LOG(BF("Setting entry->_Dihedral") );
	    entry->_Dihedral = dihedralInternal;
	}
	LOG(BF("Finished building zMatrix") );
    }



    void	ZMatrix_O::defineForMatter(Matter_sp matter)
    {_G();
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


    void	ZMatrix_O::extractInternals()
    {_G();
	ZMatrix_O::internalIterator	it;
	for ( it=this->_Internals.begin(); it!=this->_Internals.end(); it++ )
	{
	    (*it)->extractInternal();
	}
    }


    void	ZMatrix_O::constrainInternalsWithNoAtomsWithMask(uint mask)
    {_G();
	ZMatrix_O::internalIterator	it;
	for ( it=this->_Internals.begin(); it!=this->_Internals.end(); it++ )
	{
	    (*it)->constrainIfNoAtomsWithMask(mask);
	}

    }



    core::List_sp ZMatrix_O::entriesAsCons() const
    {_OF();
	return core::Cons_O::createFromVec0(this->_ZMatrix);
    };

    ZMatrixEntry_sp ZMatrix_O::zMatrixEntryWithName(string const& atomZMatrixName) const
    {_OF();
	for ( const_zMatrixEntryIterator it = this->begin_ZMatrixEntries();
	      it != this->end_ZMatrixEntries(); it++ )
	{
	    if ( (*it)->getTargetAtomZMatrixName()==atomZMatrixName) return *it;
	}
	return _Nil<ZMatrixEntry_O>();
    }


    string ZMatrix_O::allZMatrixAtomNamesAsString() const
    {_OF();
	stringstream ss;
	for ( const_zMatrixEntryIterator it = this->begin_ZMatrixEntries();
	      it != this->end_ZMatrixEntries(); it++ )
	{
	    ss << (*it)->getTargetAtomZMatrixName() << " ";
	}
	return ss.str();
    }




};
