/*
    File: Joint.h
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
#ifndef	kinematics_Joint_H
#define kinematics_Joint_H

#include <clasp/core/foundation.h>
#include <cando/geom/vector3.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/dofType.h>
#include <cando/kinematics/stub.fwd.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/coordinateCalculators.h>
#include <cando/kinematics/joint.fwd.h>
#include <cando/chem/atomId.h>
//#include <cando/kinematics/pool.h>


#define PREPARE_ANGLE(ang) (MY_PI-ang)
#define FINAL_ANGLE(ang) (ang)

// Currently a maximum of 5 children hanging off of any Joint
// means we can handle octahedral coordination but no more

namespace kinematics
{


#if 0
#define KIN_LOG(...) if (chem::chem__verbose(3)) {core::clasp_write_string(fmt::format("{}:{}:{} - " , __FILE__ , __LINE__ , __FUNCTION__ )); core::clasp_write_string(fmt::format(__VA_ARGS__));};
#else
#define KIN_LOG(...)
#endif


#define	DEBUG_KIN_JOINT	1


class RootJointInfo;

/*! Joint stores a link to a parent Joint and a list of children Joints that
  are reference counted.  The children Joints can be of different types and the order of the 
  children is enforced as JumpJoints | BondedJoint Joints.
 */

FORWARD(Joint);
class Joint_O : public core::CxxObject_O
{
  LISP_CLASS(kinematics,KinPkg,Joint_O,"Joint",core::CxxObject_O);
public:
	//! Joint to the parent Joint (also used to contruct linked list of unused PoolMembers)
  Joint_sp 	        _Parent;
  core::T_sp            _Name;
  chem::AtomId	        _Id;
  int                   _PositionIndexX3;
  int                   _EndPositionIndexX3;
  core::List_sp         _Properties;
private:
	/*! Track my position in my owner's list of Joints with modified DOFs
	  -1 when my dofs have not changed since the last update_coords
	*/
  int		_DofRefoldIndex;
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node);
public:
	/*! Bonded Joints can have different numbers of children wrt JumpJoints */
  virtual int _maxNumberOfChildren() const = 0;
	/*! Return the current number of children */
  virtual int _numberOfChildren() const = 0;
	/*! Return a reference to the indexed child */
  virtual Joint_sp _child(int idx) = 0;
	/*! Return a reference to the indexed child */
  virtual Joint_sp _child(int idx) const = 0;
	/*! Set a value of a child */
  virtual void _setChild(int idx, Joint_sp Joint) = 0;
	/*! Delete the child at the given index */
  virtual void _releaseChild(int idx) = 0;
	/*! Insert the child at the given index - this does the
	  work of opening up a space and putting the new value in */
  virtual void _insertChild(int idx, Joint_sp Joint) = 0;
	/*! Insert the child at the given index - this does the work
	  of opening up a space and putting the new value in */
  virtual void _appendChild(Joint_sp Joint) = 0;

	/*! Destructors need to delete all Children */
  virtual void _releaseAllChildren() = 0;
  virtual core::List_sp jointChildren() const;
  Joint_sp onlyOtherChild(Joint_sp child) const;
  Joint_sp childWithName(core::T_sp name) const;

  CL_DEFMETHOD int positionIndexX3() const { return this->_PositionIndexX3; };



public:

  Joint_O() : _Parent(unbound<Joint_O>()), _Name(nil<core::T_O>()), _Id() {};
  Joint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable ) :
    _Parent(unbound<Joint_O>()),
    _Name(name),
    _Id(atomId),
    _Properties(nil<core::T_O>()) {
    this->_PositionIndexX3 = this->calculatePositionIndex(atomTable);
    this->_EndPositionIndexX3 = atomTable->getNVectorSize();
#ifdef DEBUG_ASSERT
    if ((this->_PositionIndexX3%3) || this->_PositionIndexX3 >= this->_EndPositionIndexX3) {
      printf("%s:%d:%s this->_PositionIndexX3 %d is not 3 aligned or >= this->_EndPositionIndexX3 %d\n", __FILE__, __LINE__, __FUNCTION__,
             this->_PositionIndexX3, this->_EndPositionIndexX3);
    }
#endif
  };

  int calculatePositionIndex(chem::AtomTable_sp atomTable);
  
  /*! Returns true if the joint represents an Joint */
  CL_DEFMETHOD virtual bool correspondsToJoint() const { return true; };

  CL_DEFMETHOD chem::AtomId id() const { return this->_Id;};
  CL_DEFMETHOD core::T_sp name() const;

  virtual core::Symbol_sp typeSymbol() const;
  virtual string asString() const;
  virtual string __repr__() const;

  /*! Set the parent and setup the ToExternal and ToInternal CoordinateCalculators */
  void setParent(Joint_sp parent);

	/*! Return a Joint_sp for the parent */
  CL_DEFMETHOD bool parentBoundP() const { return this->_Parent.boundp(); };
  CL_DEFMETHOD Joint_sp parent() const { if (this->_Parent.boundp()) return this->_Parent; SIMPLE_ERROR("parent of {} is not bound", _rep_(this->asSmartPtr())); };
  
	/*! Insert the child before the (before) index. */
  void insertChild( int before, Joint_sp child );

	/*! Add the child.
	  If the child is a JumpJoint then put it as the first child.
	  If it's a Bonded Joint then put it before all the existing BondedJoints.
	*/
  void addChild(Joint_sp child);

  virtual bool definedp() const;

  void setToInternal(core::Symbol_sp cc);
  core::Symbol_sp getToInternal();
  void setToExternal(core::Symbol_sp cc);
  core::Symbol_sp getToExternal();
	/*! Return the index of the first non-jump child */
  int firstNonJumpChildIndex() const;


	/*! Append the child. 
	  If the child is a JumpJoint then put it before all the other non-JumpJoints.
	  If it is a Bonded Joint then put it at the end.
	*/
  void appendChild(Joint_sp handle);

	/*! Erase the child from the Children using the index into the _Children array */
  void erase(int idx );

  void eraseChild(Joint_sp child);
  int indexOfChild(Joint_sp child);



  Joint_sp child(int idx) { return this->_child(idx);};
	/*! Return a reference to the indexed child */
  Joint_sp child(int idx) const { return this->_child(idx);};


	/*! Root nodes will return RootJointInfo structures */
  virtual RootJointInfo const* rootJointInfo() const;

	/*! For debugging dump the tree */
  void recursiveDumpChildrenIntoStringStream(const string& prefix,
                                             stringstream& out);

  void updateInternalCoord(chem::NVector_sp internals, chem::NVector_sp coords);
  virtual void _updateInternalCoord(chem::NVector_sp internals, chem::NVector_sp coords) { THROW_HARD_ERROR("Subclass must implement"); };
	/*! Update the internal coordinates */
  void updateInternalCoords(chem::NVector_sp internals, chem::NVector_sp coords);

	/*! Return true if this Joint is a JumpJoint (or subclass) */
  virtual bool isJump() const { return false;};

  CL_DEFMETHOD chem::AtomId atomId() const { return this->_Id; };

	/*! Return the i(th) non-jump Joint */
  Joint_sp getNonJumpJoint(int idx) const;

	/*! Return true if the stub is defined */
  virtual bool stubDefined() const;


	/*! Return the position */
  Vector3 position(chem::NVector_sp coords) const;

	/*! Set the position */
  void setPosition(chem::NVector_sp coords, const Vector3& pos);

  CL_LISPIFY_NAME("KIN:PROPERTIES");
  CL_DEFMETHOD core::List_sp getProperties() const { return this->_Properties; };
  void clearProperty(core::Symbol_sp propertySymbol);
  void setProperty(core::Symbol_sp propertySymbol, core::T_sp value);
  bool hasProperty (core::Symbol_sp propertySymbol);
  core::T_sp getProperty(core::Symbol_sp propertySymbol);
  core::T_sp getPropertyOrDefault(core::Symbol_sp propertySymbol,core::T_sp defaultValue);

  Vector3 getPosition2() const;
  void setPosition2(const Vector3& pos);

	/*! Return the input stub Joint */
  virtual Joint_sp inputStubJoint0() const
  {_OF();
    ASSERTF(this->parent().boundp(),("Parent isn't defined"));
    return this->parent();
  }

	/*! Return the input stub Joint */
  virtual Joint_sp inputStubJoint1() const = 0;

	/*! Return the input stub Joint */
  virtual Joint_sp inputStubJoint2() const = 0;

	/*! Return true if we want to keep the Dof fixed */
  bool keepDofFixed(DofType dof) const { return false;};

	/*! Walk the tree */
  void walkChildren(core::Function_sp callback_one_arg);

	/*! Walk the tree but stay in one residue */
  void walkResidueTree(int residueId, core::Function_sp callback_one_arg);


	/*! Update the external coordinates after calculating the input stub */
  virtual void updateXyzCoords(chem::NVector_sp internals, chem::NVector_sp coords);

  /*! Update the external coordinate of just this node */
  virtual void updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords);

	/*! Update the external coordinates using the input stub */
  virtual void _updateXyzCoords(chem::NVector_sp internals, chem::NVector_sp coords, Stub& stub) {THROW_HARD_ERROR("Subclass must implement");};

  virtual void _updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords,Stub& stub) {THROW_HARD_ERROR("Subclass must implement");};

  void _updateChildrenXyzCoords(chem::NVector_sp internals, chem::NVector_sp coords);

  void applyTransformToXyzCoordsRecursively(const Matrix& transform, core::Array_sp coords);

  template <typename DestVectorType, typename SourceVectorType>
  void extractXyzCoordinatesRecursivelyImpl(DestVectorType destination, Joint_sp sourceJoint, SourceVectorType source) {
    if (this->_Name != sourceJoint->_Name) {
      SIMPLE_ERROR("extractXyzCoordinatesRecursivelyImpl encountered a mismatch of joint names {} and {}", _rep_(this->_Name), _rep_(sourceJoint->_Name));
    }
    (*destination)[this->_PositionIndexX3] = (*source)[sourceJoint->_PositionIndexX3];
    (*destination)[this->_PositionIndexX3+1] = (*source)[sourceJoint->_PositionIndexX3+1];
    (*destination)[this->_PositionIndexX3+2] = (*source)[sourceJoint->_PositionIndexX3+2];
    if (this->_numberOfChildren() != sourceJoint->_numberOfChildren()) {
      SIMPLE_ERROR("There is a child size mismatch when extractXyzCoordinatesRecursiveImpl running this -> {}  sourceJoint-> {}", _rep_(this->asSmartPtr()), _rep_(sourceJoint));
    }
    for ( int ii=0; ii < this->_numberOfChildren(); ii++) {
      this->_child(ii)->extractXyzCoordinatesRecursivelyImpl(destination,sourceJoint->_child(ii),source);
    }
  }
  void extractXyzCoordinatesRecursively(core::Array_sp dest, Joint_sp sourceJoint, core::Array_sp source);

	/*! Ensure proper function of the output-sensitive refold subroutine
	  derived classes must invoke this function during their updateXyzCoords subroutines
	*/
  void noteXyzUpToDate()
  {
    this->_DofRefoldIndex = 0;
  }


	/*! Return the input stub */
  virtual CL_DEFMETHOD Stub getInputStub(chem::NVector_sp coords) const { THROW_HARD_ERROR("Subclass must implement"); };

	/*! Return the value of the DOF */
  virtual double dof(DofType const& dof) const { THROW_HARD_ERROR("SubClass must implement");};
  core::T_sp getParentOrNil() const;

};





};


#endif
