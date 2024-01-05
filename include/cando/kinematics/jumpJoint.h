/*
    File: jumpJoint.h
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
#ifndef	kinematics_jumpJoint_H
#define kinematics_jumpJoint_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/joint.h>
#include <cando/kinematics/jump.h>
#include <cando/chem/atomId.h>

namespace kinematics
{

  FORWARD(JumpJoint);
  class JumpJoint_O : public Joint_O {
    LISP_CLASS(kinematics,KinPkg,JumpJoint_O,"JumpJoint",Joint_O);
  public:
    static const NodeType nodeType = jumpJoint;
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
  public:
    core::T_sp _Orientation; // This is the stub for children
	/*! JumpJoints can have unlimited numbers of children */
    gc::Vec0< Joint_sp >	_Children;
    static JumpJoint_sp make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable, core::T_sp orientation );
  protected:
	/*! Bonded atoms can have different numbers of children wrt JumpJoints */
    virtual int _maxNumberOfChildren() const { return INT_MAX;};
	/*! Return the current number of children */
    virtual int _numberOfChildren() const {return this->_Children.size();};
	/*! Get the value of the indexed child */
    virtual Joint_sp _child(int idx) {return this->_Children[idx];};
    virtual Joint_sp _child(int idx) const {return this->_Children[idx];};
	/*! Set the value of a child */
    virtual void _setChild(int idx,Joint_sp v) {this->_Children[idx] = v;};
	/*! Delete the child at the given index */
    virtual void _releaseChild(int idx); //  {this->_Children.erase(this->_Children.begin()+idx);};
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
    virtual void _insertChild(int idx, Joint_sp c); //  {this->_Children.insert(this->_Children.begin()+idx,c);};
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
    virtual void _appendChild(Joint_sp c); //  {this->_Children.push_back(c);};

	/*! Delete all children for dtor */
    virtual void _releaseAllChildren();

  public:
    JumpJoint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable, core::T_sp orientation );

	/*! Empty ctor */
    JumpJoint_O() {};

    // Return the 'external' coordinate
    CL_DEFMETHOD core::T_sp getOrientation() const { return this->_Orientation; };
    CL_DEFMETHOD void setOrientation(core::T_sp orientation) { this->_Orientation = orientation; };

    virtual core::Symbol_sp typeSymbol() const;

    Stub getInputStub(chem::NVector_sp coords) const;
    virtual bool internalp() const;

    virtual void _updateInternalCoord(chem::NVector_sp coords);
    
	/*! Yes, this is a JumpJoint */
    bool isJump() const { return true;};

   	/*! Return the stubJoint0 */
    virtual Joint_sp inputStubJoint0() const { return unbound<Joint_O>(); };

	/*! Return the stubJoint1 */
    virtual Joint_sp inputStubJoint1() const { return unbound<Joint_O>(); };

	/*! Return the stubJoint2 */
    virtual Joint_sp inputStubJoint2() const { return unbound<Joint_O>(); };

    Matrix transform() const;

    bool keepDofFixed(DofType dof) const;

    virtual void updateXyzCoord(chem::NVector_sp coords);

    /*! Update just this joints position */
    virtual void _updateXyzCoord(chem::NVector_sp coords,Stub& stub);

    virtual void _updateChildrenXyzCoords(chem::NVector_sp coords);

	/*! Return the DOF */
    double dof(DofType const& dof) const;
	
  };



};
#endif
