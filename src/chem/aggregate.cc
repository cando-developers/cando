/*
    File: aggregate.cc
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

#include <clasp/core/common.h>
#include <cando/chem/aggregate.h>

#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <clasp/core/array.h>
#include <clasp/core/hashTable.h>
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

CL_DOCSTRING(R"doc(Create a chem:bounding-box using three lengths in angstroms and zero to three angles (in degrees).
When angles are missing they are assumed to be 90.0 degrees.)doc");
CL_LAMBDA(widths &key angles-degrees center);
CL_LISPIFY_NAME(make-bounding-box);
CL_DEF_CLASS_METHOD
BoundingBox_sp BoundingBox_O::make(core::List_sp widths, core::T_sp angles_degrees, core::T_sp center)
{
  Vector3 vwidths(0.0,0.0,0.0);
  if (core::cl__length(widths)!=3) {
    SIMPLE_ERROR(BF("Three widths are required to define bounding-box - you supplied %s") % _rep_(widths));
  }
  core::Number_sp num;
  num = gc::As<core::Number_sp>(core::oFirst(widths));
  double xwidth = core::clasp_to_double(num);
  num = gc::As<core::Number_sp>(core::oSecond(widths));
  double ywidth = core::clasp_to_double(num);
  num = gc::As<core::Number_sp>(core::oThird(widths));
  double zwidth = core::clasp_to_double(num);
  double xangle = 90.0;
  double yangle = 90.0;
  double zangle = 90.0;
  double xcenter = 0.0;
  double ycenter = 0.0;
  double zcenter = 0.0;
  if (angles_degrees.consp()) {
    if (core::cl__length(angles_degrees)>=1) xangle = core::clasp_to_double(gc::As<core::Number_sp>(core::oFirst(angles_degrees)));
    if (core::cl__length(angles_degrees)>=2) yangle = core::clasp_to_double(gc::As<core::Number_sp>(core::oSecond(angles_degrees)));
    if (core::cl__length(angles_degrees)>=3) zangle = core::clasp_to_double(gc::As<core::Number_sp>(core::oThird(angles_degrees)));
  }
  if (center.consp()) {
    if (core::cl__length(center)>=1) xcenter = core::clasp_to_double(gc::As<core::Number_sp>(core::oFirst(center)));
    if (core::cl__length(center)>=2) ycenter = core::clasp_to_double(gc::As<core::Number_sp>(core::oSecond(center)));
    if (core::cl__length(center)>=3) zcenter = core::clasp_to_double(gc::As<core::Number_sp>(core::oThird(center)));
  }
  Vector3 width(xwidth,ywidth,zwidth);
  Vector3 angle(xangle,yangle,zangle);
  Vector3 vcenter(xcenter,ycenter,zcenter);
  auto bb = gctools::GC<BoundingBox_O>::allocate(width,angle,vcenter);
  return bb;
}

BoundingBox_O::BoundingBox_O(const Vector3& widths, const Vector3& angles_degrees, const Vector3& center)
  : _Center(center),
    _Widths(widths),
    _AnglesDegrees(angles_degrees),
    _x_rwidth(1.0/widths.getX()),
    _y_rwidth(1.0/widths.getY()),
    _z_rwidth(1.0/widths.getZ()) {};

double BoundingBox_O::distance_squared_between_two_points(const Vector3& v1, const Vector3& v2)
{
  double dx = fabs(v1.getX()-v2.getX());
  dx -= static_cast<int>(dx*this->_x_rwidth+0.5)*this->_Widths.getX();
  double dy = fabs(v1.getY()-v2.getY());
  dy -= static_cast<int>(dy*this->_y_rwidth+0.5)*this->_Widths.getY();
  double dz = fabs(v1.getZ()-v2.getZ());
  dz -= static_cast<int>(dz*this->_z_rwidth+0.5)*this->_Widths.getZ();
  return dx*dx+dy*dy+dz*dz;
}

CL_DEFUN double chem__bounding_box_distance_squared_between_two_points(BoundingBox_sp bounding_box, const Vector3& v1, const Vector3& v2)
{
  return bounding_box->distance_squared_between_two_points(v1,v2);
}



double BoundingBox_O::distance_squared_between_two_atoms(Atom_sp a1, Atom_sp a2) {
  Vector3& pos1 = a1->getPositionRef();
  Vector3& pos2 = a2->getPositionRef();
  return this->distance_squared_between_two_points(pos1,pos2);
}

CL_DEFUN double chem__bounding_box_distance_squared_between_two_atoms(BoundingBox_sp bounding_box, Atom_sp a1, Atom_sp a2)
{
  return bounding_box->distance_squared_between_two_atoms(a1,a2);
}

CL_DEFMETHOD Vector3 BoundingBox_O::get_bounding_box_widths() const {
  return this->_Widths;
}

CL_DEFMETHOD Vector3 BoundingBox_O::get_bounding_box_rwidths() const {
  Vector3 vec(this->_x_rwidth,this->_y_rwidth,this->_z_rwidth);
  return vec;
}

CL_DEFMETHOD Vector3 BoundingBox_O::get_bounding_box_angles_degrees() const {
  return this->_AnglesDegrees;
}

CL_DEFMETHOD Vector3 BoundingBox_O::get_bounding_box_center() const
{
  return this->_Center;
}

CL_DEFMETHOD Vector3 BoundingBox_O::min_corner() const
{
  Vector3 halfWidth(this->_Widths.getX()*0.5,this->_Widths.getY()*0.5,this->_Widths.getZ()*0.5);
  Vector3 min_corner = this->_Center-halfWidth;
  return min_corner;
}

CL_DEFMETHOD Vector3 BoundingBox_O::max_corner() const
{
  Vector3 halfWidth(this->_Widths.getX()*0.5,this->_Widths.getY()*0.5,this->_Widths.getZ()*0.5);
  Vector3 max_corner = this->_Center+halfWidth;
  return max_corner;
}


CL_DEFMETHOD Vector3 BoundingBox_O::normalize_position(const Vector3& pos) const {
  Vector3 min_corner = this->min_corner();
  Vector3 temp_pos = pos - min_corner;
  double xi = temp_pos.getX() - floor(temp_pos.getX()*this->_x_rwidth)*this->_Widths.getX();
  double yi = temp_pos.getY() - floor(temp_pos.getY()*this->_y_rwidth)*this->_Widths.getY();
  double zi = temp_pos.getZ() - floor(temp_pos.getZ()*this->_z_rwidth)*this->_Widths.getZ();
  Vector3 pbc_pos(xi,yi,zi);
  return pbc_pos+min_corner;
}
  

CL_DEFMETHOD double BoundingBox_O::get_x_width() const
{
  return this->_Widths.getX();
}
CL_DEFMETHOD double BoundingBox_O::get_y_width() const
{
  return this->_Widths.getY();
}
CL_DEFMETHOD double BoundingBox_O::get_z_width() const
{
  return this->_Widths.getZ();
}

CL_DEFMETHOD double BoundingBox_O::get_x_angle_degrees() const
{
  return this->_AnglesDegrees.getX();
}
CL_DEFMETHOD double BoundingBox_O::get_y_angle_degrees() const
{
  return this->_AnglesDegrees.getY();
}
CL_DEFMETHOD double BoundingBox_O::get_z_angle_degrees() const
{
  return this->_AnglesDegrees.getZ();
}

CL_LISPIFY_NAME(BoundingBox_intersects);
CL_DEFMETHOD
bool BoundingBox_O::intersects(BoundingBox_sp other)
{
  if (this->cuboidp()||!other->cuboidp()) {
    double aHalfWidthZ = this->_Widths.getZ()/2.0;
    double bHalfWidthZ = other->_Widths.getZ()/2.0;
    double aminZ = this->_Center.getZ()-aHalfWidthZ;
    double bmaxZ = other->_Center.getZ()+bHalfWidthZ;
    if (!(aminZ<=bmaxZ)) return false;
    double amaxZ = this->_Center.getZ()+aHalfWidthZ;
    double bminZ = other->_Center.getZ()-bHalfWidthZ;
    if (!(amaxZ >= bminZ)) return false;
    double aHalfWidthY = this->_Widths.getY()/2.0;
    double bHalfWidthY = other->_Widths.getY()/2.0;
    double aminY = this->_Center.getY()-aHalfWidthY;
    double bmaxY = other->_Center.getY()+bHalfWidthY;
    if (!(aminY <= bmaxY)) return false;
    double amaxY = this->_Center.getY()+aHalfWidthY;
    double bminY = other->_Center.getY()-bHalfWidthY;
    if (!(amaxY >= bminY)) return false;
    double aHalfWidthX = this->_Widths.getX()/2.0;
    double bHalfWidthX = other->_Widths.getX()/2.0;
    double aminX = this->_Center.getX()-aHalfWidthX;
    double bmaxX = other->_Center.getX()+bHalfWidthX;
    if (!(aminX <= bmaxX)) return false;
    double amaxX = this->_Center.getX()+aHalfWidthX;
    double bminX = other->_Center.getX()-bHalfWidthX;
    return (amaxX >= bminX);
  }
  SIMPLE_ERROR(BF("Handle non-cuboid bounding-box"));
}

bool BoundingBox_O::cuboidp() const
{
  return (this->_AnglesDegrees.getX()>=89.9 && this->_AnglesDegrees.getX()<=90.1);
}

void BoundingBox_O::fields(core::Record_sp node)
{
  node->field(INTERN_(kw,Center),this->_Center);
  node->field(INTERN_(kw,Widths),this->_Widths);
  node->field(INTERN_(kw,AnglesDegrees),this->_AnglesDegrees);
}

CL_DEFMETHOD
core::T_mv BoundingBox_O::get_cuboid_rsize() const
{
  if (this->cuboidp()) {
    return Values(core::clasp_make_double_float(1.0/this->_Widths.getX()),
                  core::clasp_make_double_float(1.0/this->_Widths.getY()),
                  core::clasp_make_double_float(1.0/this->_Widths.getZ()));
  }
  SIMPLE_ERROR(BF("The bounding-box is not cuboid"));
}


void Aggregate_O::fields(core::Record_sp node)
{
  node->field_if_not_unbound(INTERN_(kw,force_field_name),this->_ForceFieldName);
  node->field_if_not_unbound(INTERN_(kw,boundingBox),this->_BoundingBox);
  this->Base::fields(node);
}


CL_DOCSTRING(R"doc(Return the bounding-box for the atom-table.)doc");
CL_LISPIFY_NAME(aggregate-bounding-box);
CL_DEFMETHOD BoundingBox_sp Aggregate_O::boundingBox() const
{
  if (this->_BoundingBox.boundp()) {
    return this->_BoundingBox;
  }
  SIMPLE_ERROR(BF("The bounding-box slot is unbound"));
}

CL_DOCSTRING(R"doc(Return T if the bounding-box is bound)doc");
CL_LISPIFY_NAME(aggregate-bounding-box-bound-p);
CL_DEFMETHOD bool Aggregate_O::boundingBoxBoundP() const
{
  return this->_BoundingBox.boundp();
}

CL_DOCSTRING(R"doc(Set the bounding-box)doc");
CL_LISPIFY_NAME(aggregate-set-bounding-box);
CL_DEFMETHOD void Aggregate_O::setBoundingBox(BoundingBox_sp boundingBox) {
  this->_BoundingBox = boundingBox;
}

CL_DOCSTRING(R"doc(Make the bounding-box unbound)doc");
CL_LISPIFY_NAME(aggregate-mak-unbound-bounding-box);
CL_DEFMETHOD void Aggregate_O::makUnboundBoundingBox() {
  this->_BoundingBox = _Unbound<BoundingBox_O>();
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


Atom_sp Aggregate_O::atomWithAtomId(const AtomId& atomId) const
    {_OF();
	int molId = atomId.moleculeId();
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




Matter_sp Aggregate_O::copy(core::T_sp new_to_old)
{
  auto  newAgg = gctools::GC<Aggregate_O>::copy( *this ); // = RP_Copy<Aggregate_O>(this);
  if (gc::IsA<core::HashTable_sp>(new_to_old)) {
    core::HashTable_sp new_to_old_ht = gc::As_unsafe<core::HashTable_sp>(new_to_old);
    new_to_old_ht->setf_gethash(newAgg,this->asSmartPtr());
  }
  for ( contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ ) {
    Molecule_sp mol = gc::As<Molecule_sp>(*a);
    Molecule_sp mol_copy = gc::As<Molecule_sp>(mol->copy(new_to_old));
    newAgg->addMolecule(mol_copy);
  }
  newAgg->copyRestraintsDontRedirectAtoms(this->asSmartPtr());
  newAgg->redirectRestraintAtoms();
  return newAgg;
}





Matter_sp Aggregate_O::copyDontRedirectAtoms(core::T_sp new_to_old)
{_OF();
  auto  newAgg  = gctools::GC<Aggregate_O>::copy( *this); // = RP_Copy<Aggregate_O>(this);
  if (gc::IsA<core::HashTable_sp>(new_to_old)) {
    core::HashTable_sp new_to_old_ht = gc::As_unsafe<core::HashTable_sp>(new_to_old);
    new_to_old_ht->setf_gethash(newAgg,this->asSmartPtr());
  }
//    newAgg->duplicate(this);	// *newAgg = *this;
  for ( const_contentIterator a=this->begin_contents(); a!=this->end_contents(); a++ )
  {
    Molecule_sp mol = (*a).as<Molecule_O>();
    newAgg->addMatter(mol->copyDontRedirectAtoms(new_to_old));
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



#if 0
bool	Aggregate_O::equal(core::T_sp obj) const
    {
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
#endif

    void Aggregate_O::transferCoordinates(Matter_sp obj)
    {
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

#if 0
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
#endif


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
    for ( mi=this->begin_molecules(); mi!= this->end_molecules(); mi++ ) {
      mol = gc::As<Molecule_sp>(*mi);
//      printf("%s:%d mol->getName()= %s  matches chainSym[%s] = %d\n", __FILE__, __LINE__, _rep_(mol->getName()).c_str(),_rep_(chainSym).c_str(), mol->getName()==chainSym);
      if ( mol->getName() == chainSym ) {
        bool foundResidue = false;
        for ( ri = mol->begin_residues(); ri!=mol->end_residues(); ri++ ) {
          res = (*ri).as<Residue_O>();
//          printf("%s:%d res->getFileSequenceNumber()= %d  matches fileSequenceNumber[%d] = %d\n", __FILE__, __LINE__, res->getFileSequenceNumber(), fileSequenceNumber,  res->getFileSequenceNumber()==fileSequenceNumber);
//          printf("%s:%d     res->getId() = %d\n", __FILE__, __LINE__, res->getId());
          if ( res->getFileSequenceNumber() == fileSequenceNumber ) {
            foundResidue = true;
//            printf("%s:%d res->hasAtomWithName(%s) = %d\n",  __FILE__, __LINE__, _rep_(atomSym).c_str(), res->hasAtomWithName(atomSym));
            if ( res->hasAtomWithName(atomSym) ) {
              core::T_mv atom_mv = res->atomWithName(atomSym);
              atom = gc::As<Atom_sp>(atom_mv);
              result = core::Cons_O::create(atom,result);
//              printf("%s:%d  Adding atom to list: %s\n", __FILE__, __LINE__, _rep_(atom).c_str());
              foundAtom = true;
            }
          }
        }
          // Try the residue at the fileSequenceNumber index
//        printf("%s:%d Using molecule content array to lookup atom\n", __FILE__, __LINE__);
        if ( (fileSequenceNumber-1) < mol->length() ) {
          res = gc::As<Residue_sp>(mol->contentAt(fileSequenceNumber-1));
          if ( res->hasAtomWithName(atomSym) ) {
            atom = gc::As_unsafe<Atom_sp>(res->atomWithName(atomSym));
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
	for ( it=this->begin_molecules(); it!= this->end_molecules(); it++ ) {
          if ((*it).as<Molecule_O>() == a ) 
	    {
		this->eraseContent(it);
		return;
	    }
	}
	SIMPLE_ERROR(BF("Aggregate does not contain molecule: %s") % _rep_(a->getName()) );
    }




CL_LISPIFY_NAME("separateMolecules");
CL_DEFMETHOD     uint Aggregate_O::separateMolecules()
    {
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
CL_DEFMETHOD     Molecule_sp Aggregate_O::firstMoleculeWithAtomNamed(core::Symbol_sp name)
{
  Loop lmol(this->asSmartPtr(),MOLECULES);
  while (lmol.advance()) {
    Molecule_sp mol = lmol.getMolecule();
    Loop latm(mol,ATOMS);
    while (latm.advance()) {
      if (latm.getAtom()->getName()==name) return mol;
    }
  }
  SIMPLE_ERROR(BF("Could not find molecule with atom named %s") % _rep_(name));
}


    void Aggregate_O::addMatter(Matter_sp matter)
    {
	if ( matter.isA<Molecule_O>() ) this->Base::addMatter(matter);
	if ( matter.isA<Aggregate_O>())
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
          mol2WriteAggregateToFileName(agg,fileName,true);
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
	mol = (this->contentAt(0)).as<Molecule_O>();
	return mol;
    }

CL_LISPIFY_NAME("firstMoleculeName");
CL_DEFMETHOD     MatterName  Aggregate_O::firstMoleculeName()
    {
	LOG(BF("allocating mol") );
	Molecule_sp mol;
	Matter_sp c;
	ASSERTP(this->_contents.size()>0,"Aggregate_O::firstMoleculeName contains no molecules");
	c = this->contentAt(0);
	LOG(BF("about to downcast mol") );
	mol = (c).as<Molecule_O>();
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
    {
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
CL_DEFUN Aggregate_sp Aggregate_O::make(core::Symbol_sp name)
    {
      auto me = gctools::GC<Aggregate_O>::allocate_with_default_constructor();
	me->setName(name);
	return me;
    };


CL_DEFMETHOD void Aggregate_O::setf_force_field_name(core::T_sp name) {
  this->_ForceFieldName = name;
}

CL_DEFMETHOD core::T_sp Aggregate_O::force_field_name() const {
  return this->_ForceFieldName;
}






}; //namespace chem
