/*
    File: originJumpJoint.cc
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
#define DEBUG_LEVEL_NONE



#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/symbolTable.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/originJumpJoint.h>

namespace kinematics
{

void OriginJumpJoint_O::updateInternalCoords(bool const recursive,
                                             JointTree_sp at)
{
  if ( recursive ) {
    for (int childIdx=0; childIdx<this->_numberOfChildren(); childIdx++ ) {
      this->_child(childIdx)->updateInternalCoords(true,at);
    }
  }
}



Stub OriginJumpJoint_O::getStub() const {
  Stub origin;
  origin._Transform.setToIdentity();
  return origin;
}




core::Symbol_sp OriginJumpJoint_O::typeSymbol() const
{
  return _sym_origin;
};


void OriginJumpJoint_O::_updateXyzCoord(chem::NVector_sp coords,Stub& stub)
{
  Vector3 origin(0.0,0.0,0.0);
  this->setPosition(coords,origin);
}



};
