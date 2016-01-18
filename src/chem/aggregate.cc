#define	DEBUG_LEVEL_NONE


//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/common.h>
#include <cando/chem/aggregate.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <clasp/core/str.h>
#include <cando/chem/loop.h>
#include <cando/chem/moe.h>
#include <cando/chem/mol2.h>
//#include "core/serialize.h"
#include <clasp/core/numerics.h>
//#include "core/xmlSaveArchive.h"
//#include "core/xmlLoadArchive.h"
#include <cando/chem/chimera.h>
//#include "grPickableMatter.h"
#include <cando/chem/restraint.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/atomIdMap.h>
#include <clasp/core/wrappers.h>


namespace chem {

    string Aggregate_O::__repr__() const
    {
	return this->Base::__repr__();
    }

#ifdef XML_ARCHIVE
    Aggregate_sp	Aggregate_O::open(T_sp pathDesignator)
    {_G();
	Aggregate_sp	agg;
	int		dot;
	string		ext;
	dot = fileName.find_last_of('.');
	ext = fileName.substr(dot+1);
	if ( ext == "moe" ) {
	    agg = Aggregate_O::create();
	    moeReadAggregateFromFileName(agg,fileName.c_str(),lisp);
	} else if ( ext == "mol2" ) {
	    agg = Aggregate_O::create();
	    mol2ReadAggregateFromFileName(agg,fileName.c_str(),lisp);
	} else if ( ext == "xml" ) {
	    core::XmlLoadArchive_sp	archive;
	    archive = core::XmlLoadArchive_O::create();
	    archive->parse(fileName);
	    agg = downcast<Aggregate_O>(archive->get("aggregate"));
	} else {
	    SIMPLE_ERROR(BF("Cannot read aggregate from unknown fileType: %s") % fileName );
	}
	return agg;
    }
#endif



#ifdef OLD_SERIALIZE
    void	Aggregate_O::serialize(serialize::SNode node)
    {_G();
	this->Base::serialize(node);
    }
#endif


void Aggregate_O::fields(core::Record_sp node)
{
  this->Base::fields(node);
}


    AtomIdToAtomMap_sp Aggregate_O::buildAtomIdMap() const
    {_OF();
	AtomIdToAtomMap_sp atomIdMap = AtomIdToAtomMap_O::create();
	atomIdMap->resize(this->_contents.size());
	for ( int mid = 0; mid<(int)this->_contents.size(); mid++ )
	{
	    int numResidues = this->_contents[mid]->_contents.size();
	    atomIdMap->resize(mid,numResidues);
	    for ( int rid =0; rid<numResidues; rid++ )
	    {
		int numAtoms = this->_contents[mid]->_contents[rid]->_contents.size();
		atomIdMap->resize(mid,rid,numAtoms);
		for ( int aid=0; aid<numAtoms; aid++ )
		{
		    AtomId atomId(mid,rid,aid);
		    atomIdMap->set(atomId,this->_contents[mid]->_contents[rid]->_contents[aid].as<Atom_O>());
		}
	    }
	}
	return atomIdMap;
    }


    Atom_sp Aggregate_O::atomWithAtomId(AtomId_sp atomId) const
    {_OF();
	int molId = atomId->moleculeId();
	if ( molId >=0 && molId <=(int)this->_contents.size() )
	{
	    Molecule_sp molecule = this->_contents[molId].as<Molecule_O>();
	    return molecule->atomWithAtomId(atomId);
	}
	SIMPLE_ERROR(BF("Illegal moleculeId[%d] must be less than %d") % molId % this->_contents.size() );
    }


#if 0
    Atom_sp Aggregate_O::lookupAtom(const AtomId& atomId) const
    {_OF();
	ASSERTF((int)this->_contents.size()==this->_AtomIdMap.numberOfMolecules(),
		BF("The AtomIdMap is out of sync with the Aggregate contents - wrong number of molecules - use updateAtomIdMap to correct this"));
	this->_AtomIdMap.throwIfInvalidMoleculeId(atomId);
	int molId = atomId.moleculeId();
	ASSERTF((int)this->_contents[molId]->_contents.size()==this->_AtomIdMap.numberOfResidues(molId),
		BF("The AtomIdMap is out of sync with the Aggregate contents "
		   "- wrong number of residues in molecule[%d] - use updateAtomIdMap to correct this")
		% molId );
	this->_AtomIdMap.throwIfInvalidResidueId(atomId);
	int resId = atomId.residueId();
	ASSERTF((int)this->_contents[molId]->_contents[resId]->_contents.size()==this->_AtomIdMap.numberOfAtoms(molId,resId),
		BF("The AtomIdMap is out of sync with the Aggregate contents "
		   "- wrong number of atoms in residue[%d]s in molecule[%d] - use updateAtomIdMap to correct this")
		% resId % molId );
	this->_AtomIdMap.throwIfInvalidAtomId(atomId);
	return this->_AtomIdMap[atomId];
    }

#endif

    



    void	Aggregate_O::initialize()
    {
	this->Base::initialize();
    }





    void	Aggregate_O::duplicate(const Aggregate_O* a ) {
	*this = *a;
    }




Matter_sp Aggregate_O::copy()
{
  GC_COPY(Aggregate_O, newAgg, *this ); // = RP_Copy<Aggregate_O>(this);
  for ( contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ ) {
    Molecule_sp mol = gc::As<Molecule_sp>(*a);
    Molecule_sp mol_copy = gc::As<Molecule_sp>(mol->copy());
    newAgg->addMolecule(mol_copy);
  }
  newAgg->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
  newAgg->redirectRestraintAtoms();
  return newAgg;
}





    Matter_sp Aggregate_O::copyDontRedirectAtoms()
    {_OF();
	GC_COPY(Aggregate_O, newAgg , *this); // = RP_Copy<Aggregate_O>(this);
//    newAgg->duplicate(this);	// *newAgg = *this;
	for ( const_contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ )
	{
	    Molecule_sp mol = (*a).as<Molecule_O>();
	    newAgg->addMatter(mol->copyDontRedirectAtoms());
	}
	newAgg->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
	return newAgg;
    }



    void Aggregate_O::redirectAtoms()
    {_OF();
	LOG(BF("Aggregate_O::redirectAtoms START") );
	for ( contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ ) {
	    (*a)->redirectAtoms();
	}
	this->redirectRestraintAtoms();
	LOG(BF("Aggregate_O::redirectAtoms DONE") );
    }



    bool	Aggregate_O::equal(core::T_sp obj) const
    {_G();
	if ( this->eq(obj) ) return true;
	if ( !obj.isA<Aggregate_O>() ) return false;
	Aggregate_sp other = obj.as<Aggregate_O>();
	if ( other->getName() != this->getName() ) return false;
	if ( other->_contents.size() != this->_contents.size() ) return false;
	Matter_O::const_contentIterator tit,oit;
	for ( tit=this->_contents.begin(), oit=other->_contents.begin();
	      tit!=this->_contents.end(); tit++, oit++ )
	{
	    if ( ! (*tit)->equal(*oit) ) return false;
	}
	return true;
    }


    void Aggregate_O::transferCoordinates(Matter_sp obj)
    {_G();
	if ( !obj.isA<Aggregate_O>() ) 
	{
	    SIMPLE_ERROR(BF("You can only transfer coordinates to a Aggregate from another Aggregate"));
	}
	Aggregate_sp other = obj.as<Aggregate_O>();
	if ( other->_contents.size() != this->_contents.size() )
	{
	    SIMPLE_ERROR(BF("You can only transfer coordinates if the two Aggregates have the same number of contents"));
	}
	Matter_O::contentIterator tit,oit;
	for ( tit=this->_contents.begin(), oit=other->_contents.begin();
	      tit!=this->_contents.end(); tit++, oit++ )
	{
	    (*tit)->transferCoordinates(*oit);
	}
    }


CL_LISPIFY_NAME("atomWithId");
CL_DEFMETHOD     Atom_sp	Aggregate_O::atomWithId(int id)
    {
	Loop		la;
	Atom_sp		a;
	Aggregate_sp	me;
	me = this->sharedThis<Aggregate_O>();
	la.loopTopGoal(me,ATOMS);
	while ( la.advanceLoopAndProcess() ) {
	    a = la.getAtom();
	    if ( a->getId() == id ) {
		return a;
	    }
	}
	return Atom_sp();
    }



    bool	Aggregate_O::testForVdwOverlap(Aggregate_sp other)
    {
        gctools::Vec0<Atom_sp>	atoms;
	Atom_sp		a;
	Loop		la;
	Vector3		va, vb, vdiff;
	double		vdwa, vdwb, vdwSum, dist;

	//
	// First copy other atoms into a vararray for efficiency
	//
	la.loopTopGoal(other,ATOMS);
	while ( la.advanceLoopAndProcess() ) {
	    a = la.getAtom();
	    atoms.push_back(a);
	}
	//
	// Now check for a collision
	//
	la.loopTopGoal(this->sharedThis<Aggregate_O>(),ATOMS);
	while ( la.advanceLoopAndProcess() ) {
	    a = la.getAtom();
	    va = a->getPosition();
	    vdwa = a->getVdwRadius();
	    for ( gctools::Vec0<Atom_sp>::iterator oa=atoms.begin(); oa!=atoms.end();oa++){
		vb = (*oa)->getPosition();
		vdwb = (*oa)->getVdwRadius();
		vdiff = vb - va;
		vdwSum = vdwa+vdwb;
		if ( fabs(vdiff.getX()) > vdwSum ) continue;
		if ( fabs(vdiff.getY()) > vdwSum ) continue;
		if ( fabs(vdiff.getZ()) > vdwSum ) continue;
		dist = vdiff.length();
		if ( dist > vdwSum ) continue;
		return true;
	    }
	}
	return false;
    }

/*!
 * Return a list of atoms specified by the Chimera atom specifications
 * if atoms aren't found then substitute the atom specification as a string.
 */
CL_LISPIFY_NAME("atomsWithChimeraSpecifications");
CL_DEFMETHOD core::List_sp	Aggregate_O::atomsWithChimeraSpecifications(const string& specs)
{
  vector<string> specParts = core::split(specs," \n\t");
  core::List_sp result(_Nil<core::T_O>());
  contentIterator	mi,ri;
  Molecule_sp		mol;
  Residue_sp		res;
  Atom_sp			atom;
  uint			fileSequenceNumber;
  string			chain, atomName;
  for ( int idx(specParts.size()-1); idx>=0; --idx ) {
    parseChimeraAtomSpecification(specParts[idx],fileSequenceNumber,chain,atomName);
    MatterName chainSym = chemkw_intern(chain);
    MatterName atomSym = chemkw_intern(atomName);
//    printf("%s:%d Looking for chainSym[%s] atomSym[%s]\n", __FILE__, __LINE__, _rep_(chainSym).c_str(), _rep_(atomSym).c_str());
    bool foundAtom = false;
    for ( mi=this->getContents().begin(); mi!= this->getContents().end(); mi++ ) {
      mol = gc::As<Molecule_sp>(*mi);
//      printf("%s:%d mol->getName()= %s  matches chainSym[%s] = %d\n", __FILE__, __LINE__, _rep_(mol->getName()).c_str(),_rep_(chainSym).c_str(), mol->getName()==chainSym);
      if ( mol->getName() == chainSym ) {
        bool foundResidue = false;
        for ( ri = mol->getContents().begin(); ri!=mol->getContents().end(); ri++ ) {
          res = downcast<Residue_O>(*ri);
//          printf("%s:%d res->getFileSequenceNumber()= %d  matches fileSequenceNumber[%d] = %d\n", __FILE__, __LINE__, res->getFileSequenceNumber(), fileSequenceNumber,  res->getFileSequenceNumber()==fileSequenceNumber);
//          printf("%s:%d     res->getId() = %d\n", __FILE__, __LINE__, res->getId());
          if ( res->getFileSequenceNumber() == fileSequenceNumber ) {
            foundResidue = true;
//            printf("%s:%d res->hasAtomWithName(%s) = %d\n",  __FILE__, __LINE__, _rep_(atomSym).c_str(), res->hasAtomWithName(atomSym));
            if ( res->hasAtomWithName(atomSym) ) {
              atom = res->atomWithName(atomSym);
              result = core::Cons_O::create(atom,result);
//              printf("%s:%d  Adding atom to list: %s\n", __FILE__, __LINE__, _rep_(atom).c_str());
              foundAtom = true;
            }
          }
        }
          // Try the residue at the fileSequenceNumber index
        MatterVector& residues = mol->getContents();
//        printf("%s:%d Using molecule content array to lookup atom\n", __FILE__, __LINE__);
        if ( (fileSequenceNumber-1) < residues.size() ) {
          res = residues[fileSequenceNumber-1];
          if ( res->hasAtomWithName(atomSym) ) {
            atom = res->atomWithName(atomSym);
            result = core::Cons_O::create(atom,result);
//            printf("%s:%d  Adding atom from residue in content to list: %s\n", __FILE__, __LINE__, _rep_(atom).c_str());
            foundAtom = true;
          }
        }
      }
    }
    if (!foundAtom) result = core::Cons_O::create(core::Str_O::create(specParts[idx]),result);
  }
  return result;
}


//
//	removeMolecule
//
//	Remove the molecule
//
CL_LISPIFY_NAME("removeMolecule");
CL_DEFMETHOD     void Aggregate_O::removeMolecule( Molecule_sp a )
    {_OF();
	contentIterator	it;
	for ( it=this->getContents().begin(); it!= this->getContents().end(); it++ ) {
	    if ( downcast<Molecule_O>(*it) == a ) 
	    {
		this->eraseContent(it);
		return;
	    }
	}
	SIMPLE_ERROR(BF("Aggregate does not contain molecule: %s") % _rep_(a->getName()) );
    }




CL_LISPIFY_NAME("separateMolecules");
CL_DEFMETHOD     uint Aggregate_O::separateMolecules()
    {_G();
        gctools::SmallOrderedSet<Atom_sp>	unassignedAtoms;
        gctools::Vec0<gctools::SmallOrderedSet<Atom_sp> >	moleculeGroups;
	uint numMolecules = 0;
	{ _BLOCK_TRACEF(BF("Get all of the atoms into a set"));
	    Atom_sp a;
	    Loop lb;
	    Vector3 v,vd;
	    lb.loopTopAggregateGoal( this->sharedThis<Aggregate_O>(), ATOMS );
	    while ( lb.advanceLoopAndProcess() ) 
	    {
		unassignedAtoms.insert(lb.getAtom());
	    }
	    this->eraseContents();	// erase the old molecules
	}
	while ( unassignedAtoms.size() > 0 )
	{
	    Molecule_sp molecule = Molecule_O::create();
	    Residue_sp res = Residue_O::create();
	    molecule->addMatter(res);
	    Atom_sp top = *(unassignedAtoms.begin());
	    SpanningLoop_sp spanner = SpanningLoop_O::create(top);
	    while ( spanner->advance() )
	    {
		Atom_sp a = spanner->getAtom();
		unassignedAtoms.erase(a);
		res->addAtom(a);
	    }
	    this->addMolecule(molecule);
	    numMolecules++;
	}
	return numMolecules;
    }


CL_LISPIFY_NAME("firstMoleculeWithAtomNamed");
CL_DEFMETHOD     Molecule_sp Aggregate_O::firstMoleculeWithAtomNamed(MatterName name)
    {
	Atom_sp a = this->firstAtomWithName(name);
	Residue_sp res = a->containedBy().as<Residue_O>();
	Molecule_sp mol = res->containedBy().as<Molecule_O>();
	return mol;
    }


    void Aggregate_O::addMatter(Matter_sp matter)
    {
	if ( matter->isAssignableTo<Molecule_O>() ) this->Base::addMatter(matter);
	if ( matter->isAssignableTo<Aggregate_O>())
	{
	    for ( Matter_O::contentIterator it=matter->begin_contents(); it!=matter->end_contents(); it++ )
	    {
		this->Base::addMatter(*it);
	    }
	}
    }


    Atom_sp	Aggregate_O::firstAtomWithName( MatterName name )
    {_OF();
	Loop		lAtoms;
	Atom_sp		a;
	Aggregate_sp	agg;
	bool		gotAtom;

	agg = this->sharedThis<Aggregate_O>();
	gotAtom = false;
	lAtoms.loopTopGoal(agg,ATOMS);
	while ( lAtoms.advanceLoopAndProcess() ) 
	{
	    a = lAtoms.getAtom();
	    if ( a->getName() == name ) {
		gotAtom = true;
		break;
	    }
	}
	if ( !gotAtom ) {
          SIMPLE_ERROR(BF("Could not find atom with name: %s")% _rep_(name) );
	}
	return a;
    }



CL_LISPIFY_NAME("writeToFile");
CL_DEFMETHOD     void	Aggregate_O::writeToFile(const string& fileName)
    {_OF();
	IMPLEMENT_ME();
#ifdef XML_ARCHIVE
	Aggregate_sp	agg;
	int		dot;
	string		ext;
	dot = fileName.find_last_of('.');
	ext = fileName.substr(dot+1);
	agg = this->sharedThis<Aggregate_O>();
	if ( ext == "moe" ) {
	    moeWriteAggregateFileName(agg,fileName.c_str());
	} else if ( ext == "mol2" ) {
	    mol2WriteAggregateToFileName(agg,fileName);
	} else if ( ext == "oml" ) {
	    core::XmlSaveArchive_sp xml;
	    xml = core::XmlSaveArchive_O::create();
	    xml->put("aggregate",this->sharedThis<Aggregate_O>());
	    xml->saveAs(fileName);
	} else {
	    SIMPLE_ERROR(BF("Cannot write aggregate to unknown fileType: %s")% fileName );
	}
	return;
#endif
    }



CL_LISPIFY_NAME("firstMolecule");
CL_DEFMETHOD     Molecule_sp      Aggregate_O::firstMolecule()
    {
	Molecule_sp mol;
//    mol = downcast<Molecule_O>(this->contentAt(0));
	mol = downcast<Molecule_O>(this->contentAt(0));
	return mol;
    }

CL_LISPIFY_NAME("firstMoleculeName");
CL_DEFMETHOD     MatterName  Aggregate_O::firstMoleculeName()
    {_G();
	LOG(BF("allocating mol") );
	Molecule_sp mol;
	Matter_sp c;
	ASSERTP(this->_contents.size()>0,"Aggregate_O::firstMoleculeName contains no molecules");
	c = this->contentAt(0);
	LOG(BF("about to downcast mol") );
	mol = downcast<Molecule_O>(c);
	LOG(BF("about to return name") );
	return mol->getName();
    }


    uint	Aggregate_O::numberOfAtoms()
    {
	Atom_sp				a;
	Loop				lb;
	Vector3				v,vd;
        gctools::Vec0<Atom_sp>			atoms;
	int				numberOfAtoms;
	numberOfAtoms = 0;
	lb.loopTopAggregateGoal( this->sharedThis<Aggregate_O>(), MOLECULES );
	while ( lb.advanceLoopAndProcess() ) 
	{
	    numberOfAtoms += lb.getMolecule()->numberOfAtoms();
	}
	return numberOfAtoms;
    }


#define	MAX_ADJUST_CYCLES	10
CL_LISPIFY_NAME("perturbAtomPositions");
CL_DEFMETHOD     void	Aggregate_O::perturbAtomPositions(double dist)
    {_OF();
	Atom_sp				a;
	Loop				lb;
	double				xd,yd,zd;
	Vector3				v,vd;
        gctools::Vec0<Atom_sp>			atoms;
	int				numberOfAtoms;
	numberOfAtoms = this->numberOfAtoms();
	ASSERT_gt(numberOfAtoms,0);
	atoms.resize(numberOfAtoms);
	lb.loopTopAggregateGoal( this->sharedThis<Aggregate_O>(), ATOMS );
	int iAtom = 0;
	while ( lb.advanceLoopAndProcess() ) {
	    a = lb.getAtom();
	    atoms[iAtom] = a;
	    iAtom++;
	    xd = (core::randomNumber01()*2.0-1.0)*dist;
	    yd = (core::randomNumber01()*2.0-1.0)*dist;
	    zd = (core::randomNumber01()*2.0-1.0)*dist;
	    v = a->getPosition();
	    vd.set(xd,yd,zd);
	    v = v+vd;
	    a->setPosition(v);
	}
	bool	madeAdjustment;
	int		adjustCycles = 0;
	Vector3	vdiff, pos;
	double	diff;
        gctools::Vec0<Atom_sp>::iterator	af, an;
	while (1) {
	    madeAdjustment = false;
	    for ( af=atoms.begin();af!=atoms.end()-1;af++ ) {
		for ( an=af+1; an!=atoms.end(); an++ ) {
		    vdiff = (*af)->getPosition()-(*an)->getPosition();
		    diff = vdiff.length();
		    // If the atoms are right on top of each other then
		    // add a random vector to each
		    if ( diff<0.1 ) {
			_lisp->print(BF( "Aggregate_O::perturbAtomPositions>> Atoms are too close and being randomized!!!" ));
			xd = (core::randomNumber01()*2.0-1.0)*dist;
			yd = (core::randomNumber01()*2.0-1.0)*dist;
			zd = (core::randomNumber01()*2.0-1.0)*dist;
			v = (*af)->getPosition();
			vd.set(xd,yd,zd);
			v = v+vd;
			(*af)->setPosition(v);
			xd = (core::randomNumber01()*2.0-1.0)*dist;
			yd = (core::randomNumber01()*2.0-1.0)*dist;
			zd = (core::randomNumber01()*2.0-1.0)*dist;
			v = (*an)->getPosition();
			vd.set(xd,yd,zd);
			v = v+vd;
			(*an)->setPosition(v);
			madeAdjustment = true;
		    } else if ( diff < 0.7 ) {
			// If the atoms are less than 0.1 angstroms away from each other then
			// nudge them away from each other
			_lisp->print(BF( "Aggregate_O::perturbAtomPositions>> Atoms are close and being nudged out!!!" ));
			vdiff = vdiff.multiplyByScalar(0.5/diff);
			pos = (*af)->getPosition();
			pos = pos + vdiff;
			(*af)->setPosition(pos);
			pos = (*an)->getPosition();
			pos = pos - vdiff;
			(*an)->setPosition(pos);
			madeAdjustment = true;
		    }
		}
	    }
	    if ( !madeAdjustment) break;
	    adjustCycles++;
	    if ( adjustCycles>MAX_ADJUST_CYCLES ) {
		_lisp->print(BF("MADE MORE THAN %d ADJUSTMENT CYCLES WHEN PERTURBING ATOMS") % MAX_ADJUST_CYCLES);
		break;
	    }
	}
    }

#ifdef RENDER
    geom::Render_sp	Aggregate_O::rendered(core::List_sp options)
    {_G();
	GrPickableMatter_sp	rend;
	rend = GrPickableMatter_O::create();
	rend->setFromMatter(this->sharedThis<Aggregate_O>());
	return rend;
    }
#endif




#define ARGS_Aggregate_O_make "(&key (name :agg))"
#define DECL_Aggregate_O_make ""
#define DOCS_Aggregate_O_make "make Aggregate args: &key name"
CL_LAMBDA(&optional (name nil));
CL_LISPIFY_NAME(make-aggregate);
CL_DEFUN Aggregate_sp Aggregate_O::make(MatterName name)
    {_G();
        GC_ALLOCATE(Aggregate_O,me);
	me->setName(name);
	return me;
    };



    void Aggregate_O::exposeCando(core::Lisp_sp lisp)
    {
	core::class_<Aggregate_O>()
//	    .def_raw("core:__init__",&Aggregate_O::__init__,"(self &key name)")
//	    .def("deepCopy",&Aggregate_O::deepCopy)
//	    .def("getRestraints",&Aggregate_O::getRestraints)
	    .def("separateMolecules",&Aggregate_O::separateMolecules)
	    .def("firstMoleculeWithAtomNamed",&Aggregate_O::firstMoleculeWithAtomNamed)
	    .def("atomWithId",&Aggregate_O::atomWithId)
	    .def("addMolecule",&Aggregate_O::addMolecule)
	    .def("removeMolecule",&Aggregate_O::removeMolecule)
	    .def("firstMolecule",&Aggregate_O::firstMolecule)
	    .def("firstMoleculeName",&Aggregate_O::firstMoleculeName)
//	.def("open",&Aggregate_O::open)
	    .def("saveAs",&Aggregate_O::saveAs,"","","",false)
	    .def("writeToFile",&Aggregate_O::writeToFile)
	    .def("testAggregateConsistancy",&Aggregate_O::testAggregateConsistancy)
//	.def("parseFromXml",&Aggregate_O::parseFromXml)
//	.def("asXml",&Aggregate_O::asXml)
//	.def("asXmlWithCoordinates",&Aggregate_O::asXmlWithCoordinates)
	    .def("perturbAtomPositions",&Aggregate_O::perturbAtomPositions)
	    .def("atomsWithChimeraSpecifications", &Aggregate_O::atomsWithChimeraSpecifications )
	    .def("molecules",&Aggregate_O::molecules)
//	    .def("updateAtomIdMap",&Aggregate_O::updateAtomIdMap)
//	    .def("lookupAtom",&Aggregate_O::lookupAtom)
	    ;
//	Defun_maker(ChemPkg,Aggregate);
    }

    void Aggregate_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(ChemPkg,Aggregate,"","",_lisp)
////	    .def_raw("core:__init__",&Aggregate("(self &key name)")
//	    .def("deepCopy",&Aggregate_O::deepCopy)
//	    .def("getRestraints",&Aggregate_O::getRestraints)
	    .def("separateMolecules",&Aggregate_O::separateMolecules)
	    .def("firstMoleculeWithAtomNamed",&Aggregate_O::firstMoleculeWithAtomNamed)
	    .def("atomWithId",&Aggregate_O::atomWithId)
	    .def("addMolecule",&Aggregate_O::addMolecule)
	    .def("removeMolecule",&Aggregate_O::removeMolecule)
	    .def("firstMolecule",&Aggregate_O::firstMolecule)
	    .def("firstMoleculeName",&Aggregate_O::firstMoleculeName)
//	.def("open",&Aggregate_O::open)
	    .def("saveAs",&Aggregate_O::saveAs)
	    .def("writeToFile",&Aggregate_O::writeToFile)
	    .def("testAggregateConsistancy",&Aggregate_O::testAggregateConsistancy)
//	.def("parseFromXml",&Aggregate_O::parseFromXml)
//	.def("asXml",&Aggregate_O::asXml)
//	.def("asXmlWithCoordinates",&Aggregate_O::asXmlWithCoordinates)
	    .def("perturbAtomPositions",&Aggregate_O::perturbAtomPositions)
	    .def("atomWithChimeraSpecification", &Aggregate_O::atomWithChimeraSpecification )
	    .def("molecules",&Aggregate_O::molecules)
//	    .def("updateAtomIdMap",&Aggregate_O::updateAtomIdMap)
//	    .def("lookupAtom",&Aggregate_O::lookupAtom)
	    ;
#endif
    }


    EXPOSE_CLASS(chem,Aggregate_O);
}; //namespace chem


