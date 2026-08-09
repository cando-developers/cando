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
    static JumpJoint_sp make(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable );
  protected:
	// This class was already growable - Joint_O's Vec0 and its child methods are these,
	// promoted to the base, so all seven overrides are now redundant and gone.
  public:
    JumpJoint_O(const chem::AtomId& atomId, core::T_sp name, chem::AtomTable_sp atomTable );

	/*! Empty ctor */
    JumpJoint_O() {};

    virtual core::Symbol_sp typeSymbol() const;

    Stub getInputStub(chem::NVector_sp coords) const;


    virtual bool definedp(chem::NVector_sp internals ) const;

    virtual void _updateInternalCoord(chem::NVector_sp internals, chem::NVector_sp coords);
    
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

    virtual void updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords);

    /*! Update just this joints position */
    virtual void _updateXyzCoord(chem::NVector_sp internals, chem::NVector_sp coords,Stub& stub);

    virtual void _updateChildrenXyzCoords(chem::NVector_sp internals, chem::NVector_sp coords);

	/*! Return the DOF */
    double dof(DofType const& dof) const;
	
  };



};
#endif
