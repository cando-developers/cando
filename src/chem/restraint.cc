/*
    File: restraint.cc
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



//
//	Restraint & RestraintVector
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

void	RestraintAnchor_O::fields(core::Record_sp node)
{
//  this->Base::fields(node); // don't invoke parent fields
  node->field( INTERN_(kw,atom), this->_Atom );
  node->/*pod_*/field_if_not_default( INTERN_(kw,position), this->_Pos, Vector3());
  node->/*pod_*/field_if_not_default( INTERN_(kw,weight), this->_Weight, 0.0);
};

#define ARGS_RestraintAnchor_O_make "(&key atom pos (weight 1.0))"
#define DECL_RestraintAnchor_O_make ""
#define DOCS_RestraintAnchor_O_make "make RestraintAnchor"
RestraintAnchor_sp RestraintAnchor_O::make(Atom_sp atom, const Vector3& pos, double weight )
{
  auto  ra = gctools::GC<RestraintAnchor_O>::allocate( atom, pos, weight );
  return ra;
};


Restraint_sp	RestraintAnchor_O::copyDontRedirectAtoms()
{
  auto  rest  = gctools::GC<RestraintAnchor_O>::copy( *this); // = RP_Copy<RestraintAnchor_O>(this);
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


void RestraintAnchor_O::applyTransform(const Matrix& m) {
  this->_Pos = (m)*this->_Pos;
}



RestraintChiral_O::RestraintChiral_O(const RestraintChiral_O& old) : Restraint_O(old)
{
  this->_A= old._A;
  this->_Chirality = old._Chirality;
}


Restraint_sp	RestraintChiral_O::copyDontRedirectAtoms()
{
  auto  rest  = gctools::GC<RestraintChiral_O>::copy( *this); // = RP_Copy<RestraintChiral_O>(this);
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








void	RestraintDistance_O::fields(core::Record_sp node)
{
//  this->Base::fields(node); // don't invoke parent fields
  node->field( INTERN_(kw,a), this->_A );
  node->field( INTERN_(kw,b), this->_B );
  node->/*pod_*/field_if_not_default( INTERN_(kw,r0), this->_R0, 0.0);
  node->/*pod_*/field_if_not_default( INTERN_(kw,k), this->_K, 0.0);
};


Restraint_sp	RestraintDistance_O::copyDontRedirectAtoms()
{
  auto  rest  = gctools::GC<RestraintDistance_O>::copy( *this); // = RP_Copy<RestraintDistance_O>(this);
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
  auto  rest  = gctools::GC<RestraintAngle_O>::copy( *this); // = RP_Copy<RestraintAngle_O>(this);
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
  auto  rest  = gctools::GC<RestraintDihedral_O>::copy( *this); // = RP_Copy<RestraintDihedral_O>(this);
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
  auto  nb  = gctools::GC<RestraintFixedNonbond_O>::allocate_with_default_constructor();
  nb->setMatter(matter);
  return nb;
}


void RestraintFixedNonbond_O::initialize()
{
  this->_Matter = nil<Matter_O>();
}


};
