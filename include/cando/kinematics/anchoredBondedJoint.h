
/*
    File: anchoredBondedJoint.h
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
#include <cando/kinematics/bondedJoint.h>


namespace kinematics
{

  FORWARD(AnchoredBondedJoint);
  class AnchoredBondedJoint_O : public BondedJoint_O
  {
    LISP_CLASS(kinematics,KinPkg,AnchoredBondedJoint_O,"AnchoredBondedJoint",BondedJoint_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
    void initialize();
  public:
    static const NodeType nodeType = anchoredBondedJoint;   // NEW enum value needed
  public:
    Vector3 _ParentAtomPosition, _GrandParentAtomPosition, _GreatGrandParentAtomPosition;
    int _ParentAtomPositionIndex3 = -1;
    int _GrandParentAtomPositionIndex3 = -1;
    int _GreatGrandParentAtomPositionIndex3 = -1;
  public:
    static AnchoredBondedJoint_sp make(const chem::AtomId& atomId, core::T_sp name,
                                       chem::AtomTable_sp atomTable,
                                       int parentAtomPositionIndex3=-1,
                                       int grandParentAtomPositionIndex3=-1,
                                       int greatGrandParentAtomPositionIndex3=-1);

  AnchoredBondedJoint_O() : BondedJoint_O() {};
  AnchoredBondedJoint_O(const chem::AtomId& atomId, core::T_sp name,
                        chem::AtomTable_sp atomTable,
                        int p3, int gp3, int ggp3)
      : BondedJoint_O(atomId,name,atomTable)          // <-- own index from ITS OWN atom-table
        , _ParentAtomPositionIndex3(p3)
        , _GrandParentAtomPositionIndex3(gp3)
        , _GreatGrandParentAtomPositionIndex3(ggp3) {}

    CL_DEFMETHOD Vector3 transformedParentPos() const;
    CL_DEFMETHOD Vector3 transformedGrandParentPos() const;
    CL_DEFMETHOD Vector3 transformedGreatGrandParentPos() const;
    CL_DEFMETHOD void    updateAnchorPositionsFromCoords(chem::NVector_sp coords);

    bool    anchoredFrameP()       const override { return true; }
    Vector3 anchorParentPos()      const override { return this->transformedParentPos(); }
    Vector3 anchorGrandParentPos() const override { return this->transformedGrandParentPos(); }

    virtual Stub getInputStub(chem::NVector_sp coords) const;
  };


};
