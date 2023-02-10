/*
    File: complexBondedJoint.h
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
#ifndef	kinematics_complexBondedJoint_H
#define kinematics_complexBondedJoint_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/joint.h>
#include <cando/kinematics/bondedJoint.h>



namespace kinematics
{

FORWARD(ComplexBondedJoint);
class ComplexBondedJoint_O : public BondedJoint_O
{
  LISP_CLASS(kinematics,KinPkg,ComplexBondedJoint_O,"ComplexBondedJoint",BondedJoint_O);
public:
  bool fieldsp() { return true; };
  void fields(core::Record_sp node); 
  void initialize();
public:
  static const int MaxInputStubJoints = 2;
  Joint_sp         _InputStubJoints[MaxInputStubJoints];
public:
  static ComplexBondedJoint_sp make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable );
public:
  ComplexBondedJoint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable ) : BondedJoint_O(atomId,name,atomTable) {};
  ComplexBondedJoint_O() {};
public:
  virtual Stub getInputStub(chem::NVector_sp coords) const;

  	/*! Return the stubJoint1 */
  virtual Joint_sp inputStubJoint0() const { return this->parent();}
  
	/*! Return the stubJoint2 */
  virtual Joint_sp inputStubJoint1() const { return this->_InputStubJoints[0]; };
  
	/*! Return the stubJoint3 */
  virtual Joint_sp inputStubJoint2() const { return this->_InputStubJoints[1]; };

  bool inputStubJoint0BoundP() const;
  void setInputStubJoint0(Joint_sp joint);
  void makeUnboundInputStubJoint0();
  bool inputStubJoint1BoundP() const;
  void setInputStubJoint1(Joint_sp joint);
  void makeUnboundInputStubJoint1();
  bool inputStubJoint2BoundP() const;
  void setInputStubJoint2(Joint_sp joint);
  void makeUnboundInputStubJoint2();

  virtual void _updateInternalCoord(chem::NVector_sp coords);
  void _updateChildrenXyzCoords(chem::NVector_sp coords);

};
};
#endif
