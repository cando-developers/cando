/*
    File: bondedJoint.h
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
#ifndef	kinematics_bondedJoint_H
#define kinematics_bondedJoint_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/joint.h>



namespace kinematics
{

FORWARD(BondedJoint);
    class BondedJoint_O : public Joint_O
    {
	LISP_CLASS(kinematics,KinPkg,BondedJoint_O,"BondedJoint",Joint_O);
    public:
        bool fieldsp() { return true; };
        void fields(core::Record_sp node); 
        void initialize();
    public:
	static const NodeType nodeType = bondedJoint;
	static const int MaxChildren = 5;
    public:
	int		_NumberOfChildren;
        // _Children have the value unbound unbound
	Joint_sp	_Children[MaxChildren];
	Real		_Phi;
	Real		_Theta;
	Real		_Distance;
	bool		_DofChangePropagatesToYoungerSiblings;
    public:
      static BondedJoint_sp make(const chem::AtomId& atomId, core::T_sp name);
    public:
	/*! Bonded atoms can have different numbers of children wrt JumpJoints */
	virtual int _maxNumberOfChildren() const { return MaxChildren;};
	/*! Return the current number of children */
	virtual int _numberOfChildren() const {return this->_NumberOfChildren;};
	/*! Return a reference to the indexed child */
	virtual Joint_sp _child(int idx) {return this->_Children[idx];};
	virtual Joint_sp _child(int idx) const {return this->_Children[idx];};
	/*! Set a child */
	virtual void _setChild(int idx, Joint_sp h) { this->_Children[idx] = h; };
	/*! Delete the child at the given index */
	virtual void _releaseChild(int idx);
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
	virtual void _insertChild(int idx, Joint_sp c );
	/*! Insert the child at the given index - this does the work of opening up a space and putting the new value in */
      virtual void _appendChild(Joint_sp c);
	/*! Delete all of the children for the destructor */
	virtual void _releaseAllChildren();

    public:
    BondedJoint_O() : Joint_O(), _NumberOfChildren(0), _Children{unbound<Joint_O>(),unbound<Joint_O>(),unbound<Joint_O>(),unbound<Joint_O>(),unbound<Joint_O>()} {};
    BondedJoint_O(const chem::AtomId& atomId, core::T_sp name) : Joint_O(atomId,name), _NumberOfChildren(0) {};

	virtual core::Symbol_sp typeSymbol() const;

	/*! Return the stubJoint1 */
      virtual Joint_sp inputStubJoint0() const;

	/*! Return the stubJoint2 */
      virtual Joint_sp inputStubJoint1() const;

	/*! Return the stubJoint3 */
      virtual Joint_sp inputStubJoint2() const;

      virtual void _updateInternalCoord();

	bool keepDofFixed(DofType dof) const;


	string asString() const;

      void _updateChildrenXyzCoords();

	/*! Update the external coordinates using the input stub */
      virtual void _updateXyzCoord(Stub& stub);

      /*! Get the stub and update the XYZ coordinate */
      void updateXyzCoord();
      
      virtual Stub getInputStub() const;

	/*! Geta the value of the DOF */
	double dof(DofType const& dof) const;

      CL_NAME(KIN:BONDED-JOINT/GET-DISTANCE);
      CL_DEFMETHOD double getDistance() const { return this->_Distance; }
      CL_NAME(KIN:BONDED-JOINT/GET-ANGLE);
      CL_DEFMETHOD double getTheta() const { return this->_Theta; };
      CL_NAME(KIN:BONDED-JOINT/GET-DIHEDRAL);
      CL_DEFMETHOD double getPhi() const { return this->_Phi; };

      CL_DEFMETHOD void setDistance(double distance) { this->_Distance = distance;}
      CL_DEFMETHOD void setTheta(double angle) { this->_Theta = angle;}
      CL_DEFMETHOD void setPhi(double dihedral) { this->_Phi = dihedral; };


    };



};
#endif
