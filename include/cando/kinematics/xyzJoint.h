/*
    File: xyzJoint.h
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
//#ifndef	kinematics_xyzJoint_H
//#define kinematics_xyzJoint_H

#pragma once


#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/joint.h>



namespace kinematics
{

FORWARD(XyzJoint);
class XyzJoint_O : public Joint_O
{
  LISP_CLASS(kinematics,KinPkg,XyzJoint_O,"XyzJoint",Joint_O);
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node); 
  void initialize();
public:
  static const NodeType nodeType = xyzJoint;
public:
  Vector3       _AtomPosition;
  int           _AtomPositionIndex3 = -1; // >= 0 : source slot in a coords vector to refresh _AtomPositon from; -1: set directly
public:
  static XyzJoint_sp make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable, const Vector3& atomPosition );
public:
	// Child storage and all seven child methods now live on Joint_O over a growable
	// gctools::Vec0<Joint_sp>.  This class kept its own fixed array of 5 and a
	// _NumberOfChildren counter; both are gone, and so is _maxNumberOfChildren.
  XyzJoint_O() : Joint_O() {};

  XyzJoint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable, const Vector3& atomPosition )
      : Joint_O(atomId,name,atomTable)
      , _AtomPosition(atomPosition)
        {};

  XyzJoint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable, int atomPositionIndex3 )
      : Joint_O(atomId,name,atomTable)
      , _AtomPositionIndex3(atomPositionIndex3)
        {};

  virtual core::Symbol_sp typeSymbol() const;

	/*! Return the stubJoint1 */
  virtual Joint_sp inputStubJoint0() const;

	/*! Return the stubJoint2 */
  virtual Joint_sp inputStubJoint1() const;

	/*! Return the stubJoint3 */
  virtual Joint_sp inputStubJoint2() const;

  void updateAtomPositionFromCoords(chem::NVector_sp coords);
  
  virtual void _updateInternalCoord(chem::NVector_sp internals, chem::NVector_sp coords);

  bool keepDofFixed(DofType dof) const;


  string asString() const;

  void _updateChildrenXyzCoords(chem::NVector_sp internals, chem::NVector_sp coords);

	/*! Update the external coordinates using the input stub */
  virtual void _updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords,Stub& stub);

      /*! Get the stub and update the XYZ coordinate */
  void updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords);

  virtual Stub getInputStub(chem::NVector_sp coords) const;
  CL_DEFMETHOD virtual bool definedp(chem::NVector_sp internals) const;

	/*! Geta the value of the DOF */
  double dof(DofType const& dof) const;

  void setFixedPosition(const Vector3& pos);
  Vector3 getFixedPosition() const;

  CL_DEFMETHOD Vector3 transformedPos() const;
  CL_DEFMETHOD Vector3 untransformedPos() const;

};



};


namespace kinematics
{

FORWARD(StubJoint);
class StubJoint_O : public XyzJoint_O
{
  LISP_CLASS(kinematics,KinPkg,StubJoint_O,"StubJoint",XyzJoint_O);
public:
  bool fieldsp() const { return true; };
  void fields(core::Record_sp node); 
  void initialize();
public:
  static const NodeType nodeType = stubJoint;
public:
  Vector3       _ParentAtomPosition;
  Vector3       _GrandParentAtomPosition;
  Vector3       _GreatGrandParentAtomPosition;
  int           _ParentAtomPositionIndex3 = -1;
  int           _GrandParentAtomPositionIndex3 = -1;
  int           _GreatGrandParentAtomPositionIndex3 = -1;
public:
  static StubJoint_sp make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable
                           ,int atomPositionIndex3
                           ,int parentAtomPositionIndex3
                           ,int grandParentAtomPositionIndex3
                           ,int greatGrandParentAtomPositionIndex3
                           );
public:
  StubJoint_O() : XyzJoint_O() {};
  StubJoint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable
                            ,int atomPositionIndex3
                            ,int parentAtomPositionIndex3
                            ,int grandParentAtomPositionIndex3
                            ,int greatGrandParentAtomPositionIndex3
              )
      : XyzJoint_O(atomId,name,atomTable,atomPositionIndex3)
      , _ParentAtomPositionIndex3(parentAtomPositionIndex3), _GrandParentAtomPositionIndex3(grandParentAtomPositionIndex3), _GreatGrandParentAtomPositionIndex3(greatGrandParentAtomPositionIndex3) {}

  Vector3 transformedParentPos() const;
  Vector3 transformedGrandParentPos() const;
  Vector3 transformedGreatGrandParentPos() const;

  bool    anchoredFrameP()       const override { return true; }
  Vector3 anchorParentPos()      const override { return this->transformedParentPos(); }
  Vector3 anchorGrandParentPos() const override { return this->transformedGrandParentPos(); }
  
  Joint_sp parent() const { SIMPLE_ERROR("Never ask for the parent of {}", _rep_(this->asSmartPtr()) ); };
  void _updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords, Stub& stub);
  void updateStubPositionsFromCoords(chem::NVector_sp coords);
};



};
//#endif
