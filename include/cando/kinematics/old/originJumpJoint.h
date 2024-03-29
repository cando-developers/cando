/*
    File: originJumpJoint.h
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
#ifndef	kinematics_originJumpJoint_H
#define kinematics_originJumpJoint_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/joint.h>
#include <cando/kinematics/jump.h>
#include <cando/kinematics/jumpJoint.h>
#include <cando/chem/atomId.h>

namespace kinematics
{

  FORWARD(OriginJumpJoint);
  class OriginJumpJoint_O : public JumpJoint_O {
    LISP_CLASS(kinematics,KinPkg,OriginJumpJoint_O,"OriginJumpJoint",JumpJoint_O);
  public:
    static const NodeType nodeType = originJumpJoint;
  public:
	/*! Empty ctor */
  OriginJumpJoint_O(const chem::AtomId& atomId, core::T_sp name, const string& comment) : JumpJoint_O(atomId,name,comment) {};

    virtual bool correspondsToAtom() const { return false; };

    virtual core::Symbol_sp typeSymbol() const;
    Stub getStub() const;

	/*! Update the internal coordinates */
    virtual void updateInternalCoords(chem::NVector_sp coords,
                                      bool const recursive,
                                      JointTree_sp at);

    void _updateXyzCoord(chem::NVector_sp coords,Stub& stub);


  };



};
#endif
