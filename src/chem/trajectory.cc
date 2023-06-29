/*
    File: trajectory.cc
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

#include <clasp/core/foundation.h>
#include <clasp/core/hashTableEq.h>
#include <cando/chem/trajectory.h>
//#include "core/archiveNode.h"
//#include "core/archive.h"
#include <cando/chem/matter.h>
#include <cando/chem/atom.h>
#include <cando/chem/loop.h>
#include <cando/geom/coordinateArray.h>

// last include is wrappers.h
#include <clasp/core/wrappers.h>


namespace chem {

void	TrajectoryFrame_O::initialize()
{
    this->Base::initialize();
}

#ifdef XML_ARCHIVE
    void	TrajectoryFrame_O::archiveBase(core::ArchiveP node)
{
    this->Base::archiveBase(node);
    node->attribute("coords",this->_Coordinates);
}
#endif


void TrajectoryFrame_O::fillFromMatter(gctools::Vec0<Atom_sp>& atomList )
{
  this->_Coordinates = core::SimpleVector_float_O::make(atomList.size()*3);
  size_t ci;
  gctools::Vec0<Atom_sp>::iterator ai;
  for ( ai=atomList.begin(), ci = 0; ai!=atomList.end(); ai++, ci++ ) {
    const Vector3& pos = (*ai)->getPosition();
    (*this->_Coordinates)[ci*3+0] = pos.getX();
    (*this->_Coordinates)[ci*3+1] = pos.getY();
    (*this->_Coordinates)[ci*3+2] = pos.getZ();
  }
}



void	TrajectoryFrame_O::applyToMatter(gctools::Vec0<Atom_sp>& atomList )
{
  ASSERT_eq((*this->_Coordinates).length(),atomList.size());
  size_t ci;
  gctools::Vec0<Atom_sp>::iterator ai;
  for ( ai=atomList.begin(), ci = 0; ai!=atomList.end(); ai++, ci++ ) {
    Vector3 pos((*this->_Coordinates)[ci*3+0],
                (*this->_Coordinates)[ci*3+1],
                (*this->_Coordinates)[ci*3+2] );
    (*ai)->setPosition(pos);
  }
}

void	Trajectory_O::initialize()
{
    this->Base::initialize();
    this->_Matter = nil<Matter_O>();
    this->_Namespace = core::HashTableEq_O::create_default();
}

void	Trajectory_O::_setupAtomList(Matter_sp matter)
{
    this->_Matter = matter;
    this->_AtomList.clear();
    Loop lAtoms;
    lAtoms.loopTopGoal(this->_Matter,ATOMS);
    while ( lAtoms.advance() )
    {
	this->_AtomList.push_back(lAtoms.getAtom());
    }
    this->_Frames.clear();
}

Trajectory_sp Trajectory_O::makeTrajectory(Matter_sp matter) {
  auto traj = gctools::GC<Trajectory_O>::allocate();
  traj->_setupAtomList(matter);
  return traj;
}

void TrajectoryFrame_O::fields(core::Record_sp node)
{
  node->/*pod_*/field( INTERN_(kw,coords), this->_Coordinates);
}

void Trajectory_O::fields(core::Record_sp node)
{
  node->/*pod_*/field( INTERN_(kw,matter), this->_Matter);
  node->/*pod_*/field( INTERN_(kw,atomlist), this->_AtomList);
  node->/*pod_*/field( INTERN_(kw,frames), this->_Frames);
  node->/*pod_*/field( INTERN_(kw,namespace), this->_Namespace);
}

CL_LISPIFY_NAME("addFrame");
CL_DEFMETHOD TrajectoryFrame_sp Trajectory_O::addFrame(core::T_sp coord_object)
{
  TrajectoryFrame_sp frame = TrajectoryFrame_O::create();
  if (gc::IsA<Matter_sp>(coord_object)) {
    ASSERTP( coord_object == this->_Matter,"The matter argument must match the Matter used to define this trajectory");
    frame->fillFromMatter(this->_AtomList);
  } else if (gc::IsA<geom::SimpleVectorCoordinate_sp>(coord_object)) {
    auto vec = gc::As_unsafe<geom::SimpleVectorCoordinate_sp>(coord_object);
    if (this->_AtomList.size() != vec->length()) {
      SIMPLE_ERROR("The coord_object is the wrong size - it contains {} vectors - the topology represents {} atoms", vec->length(), this->_AtomList.size() );
    }
    core::SimpleVector_float_sp vec3 = core::SimpleVector_float_O::make(this->_AtomList.size()*3);
    for ( size_t ii = 0; ii<this->_AtomList.size(); ii++ ) {
      (*vec3)[ii*3+0] = (*vec)[ii].getX();
      (*vec3)[ii*3+1] = (*vec)[ii].getY();
      (*vec3)[ii*3+2] = (*vec)[ii].getZ();
    }
    frame->_Coordinates = vec3;
  } else if (gc::IsA<NVector_sp>(coord_object)) {
    auto vecd = gc::As_unsafe<NVector_sp>(coord_object);
    if (this->_AtomList.size()*3 != vecd->length()) {
      SIMPLE_ERROR("The coordinates object is the wrong size - it has {} doubles - the trajectory object has {} atoms and each atom needs x,y,z coordinate", vecd->length(), this->_AtomList.size() );
    }
    core::SimpleVector_float_sp vecf = core::SimpleVector_float_O::make(this->_AtomList.size()*3);
    for ( size_t ii = 0; ii<this->_AtomList.size()*3; ii++ ) {
      (*vecf)[ii] = (*vecd)[ii];
    }
    frame->_Coordinates = vecf;
  } else if (gc::IsA<core::SimpleVector_float_sp>(coord_object) ) {
    auto vec = gc::As_unsafe<core::SimpleVector_float_sp>(coord_object);
    if (this->_AtomList.size()*3 != vec->length()) {
      SIMPLE_ERROR("The coordinates object is the wrong size - it has {} doubles - the trajectory object has {} atoms and each atom needs x,y,z coordinate", vec->length(), this->_AtomList.size() );
    }
    frame->_Coordinates = vec;
  } else {
    SIMPLE_ERROR("Cannot add-frame with argument {} of type {}", coord_object, _rep_(core::lisp_instance_class(coord_object)));
  } 
  this->_Frames.push_back(frame);
  return frame;
}


CL_LISPIFY_NAME("numberOfTrajectoryFrames");
CL_DEFMETHOD uint Trajectory_O::numberOfTrajectoryFrames()
{
    return this->_Frames.size();
}


CL_LISPIFY_NAME("getTrajectoryFrame");
CL_DEFMETHOD TrajectoryFrame_sp Trajectory_O::getTrajectoryFrame(uint i)
{
    ASSERT_lt(i,this->_Frames.size());
    return this->_Frames[i];
}

#if 0
geom::Render_sp Trajectory_O::rendered(core::List_sp kargs)
{
    geom::FrameList_sp frames = geom::FrameList_O::create();
    for ( gctools::Vec0<TrajectoryFrame_sp>::iterator it=this->_Frames.begin(); it!=this->_Frames.end(); it++ )
    {
	this->applyTrajectoryFrameToMatter(*it);
	frames->append(this->_Matter->rendered(kargs));
    }
    return frames;
}
#endif

CL_LISPIFY_NAME("applyTrajectoryFrameToMatter");
CL_DEFMETHOD void Trajectory_O::applyTrajectoryFrameToMatter(TrajectoryFrame_sp f)
{
    f->applyToMatter(this->_AtomList);
}



};
