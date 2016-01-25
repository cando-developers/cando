#define	DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//



//
//	Restraint & RestraintList
//
//	Maintains Restraints
//
//	Restraints are almost identical to MOE restraints,
//	they store 2,3 or 4 atom references, a _MinRad/_MaxRad value pair and a _Weight
//
//	Restraints can be of three types: distance, angle, dihedral
//	They can also be absolute or not (relative)
//
//	If they are absolute then the atom references are MOE ID's
//	If they are not absolute then the atom references are AtomHandles
//	within a residue
//
//
#include <stdlib.h>
//#include	<malloc.h>
#include <cando/chem/matter.h>
#include <cando/geom/ovector3.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/molecule.h>
#include <cando/chem/residue.h>
#include <cando/chem/atom.h>
#include <cando/chem/restraint.h>
//#include "core/archiveNode.h"
//#include "core/serialize.h"
//#include "core/archive.h"
#include <clasp/core/wrappers.h>



namespace chem {







#ifdef XML_ARCHIVE
    void	Restraint_O::archiveBase(core::ArchiveP node)
{
    // nothing
}
#endif

#ifdef OLD_SERIALIZE
    void	Restraint_O::serialize(serialize::SNode node)
{
    // nothing
}
#endif

void	RestraintAnchor_O::fields(core::Record_sp node)
{
//  this->Base::fields(node); // don't invoke parent fields
  node->field_if_not_nil( INTERN_(kw,atom), this->_Atom );
  node->/*pod_*/field_if_not_default( INTERN_(kw,position), this->_Pos, Vector3());
  node->/*pod_*/field_if_not_default( INTERN_(kw,weight), this->_Weight, 0.0);
};

#define ARGS_RestraintAnchor_O_make "(&key atom pos (weight 1.0))"
#define DECL_RestraintAnchor_O_make ""
#define DOCS_RestraintAnchor_O_make "make RestraintAnchor"
RestraintAnchor_sp RestraintAnchor_O::make(Atom_sp atom, const Vector3& pos, double weight )
  {
    GC_ALLOCATE_VARIADIC(RestraintAnchor_O, ra, atom, pos, weight );
    return ra;
  };


Restraint_sp	RestraintAnchor_O::copyDontRedirectAtoms()
{
    GC_COPY(RestraintAnchor_O, rest , *this); // = RP_Copy<RestraintAnchor_O>(this);
    rest->_Atom = this->getAtom();
    return rest;
}

void	RestraintAnchor_O::redirectAtoms()
{
    this->_Atom = this->getAtom()->getCopyAtom();
}


#ifdef XML_ARCHIVE
    void	RestraintAnchor_O::archiveBase(core::ArchiveP node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archivePlainObject<Vector3>("AnchorPos","AnchorPos",this->_AnchorPos);
    node->attribute("weight",this->_Weight);
}
#endif


#ifdef OLD_SERIALIZE
    void	RestraintAnchor_O::serialize(serialize::SNode node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->namedPOD("AnchorPos",this->_AnchorPos);
    node->attribute("weight",this->_Weight);
}
#endif



void	RestraintAnchor_O::invertStereochemistryOfRestraint()
{
    this->_Pos.set(this->_Pos.getX(),this->_Pos.getY(),this->_Pos.getZ()*-1.0);
}





RestraintChiral_O::RestraintChiral_O(const RestraintChiral_O& old) : Restraint_O(old)
{
    this->_A= old._A;
    this->_Chirality = old._Chirality;
}


Restraint_sp	RestraintChiral_O::copyDontRedirectAtoms()
{
    GC_COPY(RestraintChiral_O, rest , *this); // = RP_Copy<RestraintChiral_O>(this);
    rest->_A = this->getAtomA();
    return rest;
}

void	RestraintChiral_O::redirectAtoms()
{
    this->_A = this->getAtomA()->getCopyAtom();
}


#ifdef XML_ARCHIVE
    void	RestraintChiral_O::archiveBase(core::ArchiveP node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->attribute("Chirality",this->_Chirality);
}
#endif


#ifdef OLD_SERIALIZE
    void	RestraintChiral_O::serialize(serialize::SNode node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->attribute("Chirality",this->_Chirality);
}
#endif

void	RestraintChiral_O::invertStereochemistryOfRestraint()
{
    this->_Chirality *= -1.0;
}










RestraintDistance_O::RestraintDistance_O(const RestraintDistance_O& old) : Restraint_O(old)
{
    this->_A = old._A;
    this->_B = old._B;
    this->_Min = old._Min;
    this->_Max = old._Max;
    this->_Weight = old._Weight;
}

Restraint_sp	RestraintDistance_O::copyDontRedirectAtoms()
{
    GC_COPY(RestraintDistance_O, rest , *this); // = RP_Copy<RestraintDistance_O>(this);
    rest->_A = this->getAtomA();
    rest->_B = this->getAtomB();
    return rest;
}

void	RestraintDistance_O::redirectAtoms()
{
    this->_A = this->getAtomA()->getCopyAtom();
    this->_B = this->getAtomB()->getCopyAtom();
}






#ifdef XML_ARCHIVE
    void	RestraintDistance_O::archiveBase(core::ArchiveP node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archiveWeakPointer("AtomB",this->_B);
    node->attribute("Min",this->_Min);
    node->attribute("Max",this->_Max);
    node->attribute("Weight",this->_Weight);
}
#endif


#ifdef OLD_SERIALIZE
    void	RestraintDistance_O::serialize(serialize::SNode node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archiveWeakPointer("AtomB",this->_B);
    node->attribute("Min",this->_Min);
    node->attribute("Max",this->_Max);
    node->attribute("Weight",this->_Weight);
}
#endif

RestraintAngle_O::RestraintAngle_O(const RestraintAngle_O& old) : Restraint_O(old)
{
    this->_A = old._A;
    this->_B = old._B;
    this->_C = old._C;
    this->_MinRad = old._MinRad;
    this->_MaxRad = old._MaxRad;
    this->_Weight = old._Weight;
}

Restraint_sp	RestraintAngle_O::copyDontRedirectAtoms()
{
    GC_COPY(RestraintAngle_O, rest , *this); // = RP_Copy<RestraintAngle_O>(this);
    rest->_A = this->getAtomA();
    rest->_B = this->getAtomB();
    rest->_C = this->getAtomC();
    return rest;
}


void	RestraintAngle_O::redirectAtoms()
{
    this->_A = this->getAtomA()->getCopyAtom();
    this->_B = this->getAtomB()->getCopyAtom();
    this->_C = this->getAtomC()->getCopyAtom();
}

#ifdef XML_ARCHIVE
    void	RestraintAngle_O::archiveBase(core::ArchiveP node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archiveWeakPointer("AtomB",this->_B);
    node->archiveWeakPointer("AtomC",this->_C);
    node->attribute("MinRad",this->_MinRad);
    node->attribute("MaxRad",this->_MaxRad);
    node->attribute("Weight",this->_Weight);
}
#endif


#ifdef OLD_SERIALIZE
    void	RestraintAngle_O::serialize(serialize::SNode node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archiveWeakPointer("AtomB",this->_B);
    node->archiveWeakPointer("AtomC",this->_C);
    node->attribute("MinRad",this->_MinRad);
    node->attribute("MaxRad",this->_MaxRad);
    node->attribute("Weight",this->_Weight);
}
#endif


RestraintDihedral_O::RestraintDihedral_O(const RestraintDihedral_O& old) : Restraint_O(old)
{
    this->_A = old._A;
    this->_B = old._B;
    this->_C = old._C;
    this->_D = old._D;
    this->_MinDegrees = old._MinDegrees;
    this->_MaxDegrees = old._MaxDegrees;
    this->_Weight = old._Weight;
}

Restraint_sp	RestraintDihedral_O::copyDontRedirectAtoms()
{_OF();
    GC_COPY(RestraintDihedral_O, rest , *this); // = RP_Copy<RestraintDihedral_O>(this);
    rest->_A = this->getAtomA();
    rest->_B = this->getAtomB();
    rest->_C = this->getAtomC();
    rest->_D = this->getAtomD();
    return rest;
}

void	RestraintDihedral_O::redirectAtoms()
{
    this->_A = this->getAtomA()->getCopyAtom();
    this->_B = this->getAtomB()->getCopyAtom();
    this->_C = this->getAtomC()->getCopyAtom();
    this->_D = this->getAtomD()->getCopyAtom();
}

#ifdef XML_ARCHIVE
    void	RestraintDihedral_O::archiveBase(core::ArchiveP node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archiveWeakPointer("AtomB",this->_B);
    node->archiveWeakPointer("AtomC",this->_C);
    node->archiveWeakPointer("AtomD",this->_D);
    node->attribute("MinDegrees",this->_MinDegrees);
    node->attribute("MaxDegrees",this->_MaxDegrees);
    node->attribute("Weight",this->_Weight);
}
#endif


#ifdef OLD_SERIALIZE
    void	RestraintDihedral_O::serialize(serialize::SNode node)
{
    node->archiveWeakPointer("AtomA",this->_A);
    node->archiveWeakPointer("AtomB",this->_B);
    node->archiveWeakPointer("AtomC",this->_C);
    node->archiveWeakPointer("AtomD",this->_D);
    node->attribute("MinDegrees",this->_MinDegrees);
    node->attribute("MaxDegrees",this->_MaxDegrees);
    node->attribute("Weight",this->_Weight);
}
#endif

void	RestraintDihedral_O::invertStereochemistryOfRestraint()
{
    this->_MinDegrees = 360.0 - this->_MinDegrees;
    this->_MaxDegrees = 360.0 - this->_MaxDegrees;
}









RestraintFixedNonbond_sp RestraintFixedNonbond_O::create(Matter_sp matter)
{
    GC_ALLOCATE(RestraintFixedNonbond_O, nb );
    nb->setMatter(matter);
    return nb;
}


void RestraintFixedNonbond_O::initialize()
{
    this->_Matter = _Nil<Matter_O>();
}

#ifdef XML_ARCHIVE
    void RestraintFixedNonbond_O::archiveBase(core::ArchiveP node)
{
    node->attribute("matter",this->_Matter);
}
#endif


#ifdef OLD_SERIALIZE
    void RestraintFixedNonbond_O::serialize(serialize::SNode node)
{
    node->attribute("matter",this->_Matter);
}
#endif




//
//
//






#ifdef XML_ARCHIVE
    void	RestraintList_O::archiveBase(core::ArchiveP node)
{
    node->archiveVector0("Entries",this->_Restraints);
}
#endif


#ifdef OLD_SERIALIZE
    void	RestraintList_O::serialize(serialize::SNode node)
{
    node->namedVector0("Entries",this->_Restraints);
}
#endif


RestraintList_sp RestraintList_O::copyDontRedirectAtoms()
{
    RestraintList_sp newList = RestraintList_O::create();
    RestraintList_O::iterator it;
    for ( it=this->_Restraints.begin(); it!=this->_Restraints.end(); it++ )
    {
	newList->_Restraints.push_back((*it)->copyDontRedirectAtoms());
    }
    return newList;
}



void	RestraintList_O::redirectAtoms()
{
    RestraintList_O::iterator it;
    for ( it=this->_Restraints.begin(); it!=this->_Restraints.end(); it++ )
    {
	(*it)->redirectAtoms();
    }
}



void	RestraintList_O::merge(RestraintList_sp rl)
{
    gctools::Vec0<Restraint_sp>::iterator	ri;
    Restraint_sp			rr;
    if ( rl.nilp() ) return;
    for ( ri=rl->_Restraints.begin(); ri!=rl->_Restraints.end(); ri++ ) {
	rr = *ri;
	this->addRestraint(rr);
    }
}










































};

