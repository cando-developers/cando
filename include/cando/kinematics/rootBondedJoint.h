/*
    File: rootBondedJoint.h
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
#ifndef	kinematics_rootBondedJoint_H
#define kinematics_rootBondedJoint_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/chem/atomId.h>
#include <cando/kinematics/rootJointInfo.h>
#include <cando/kinematics/bondedJoint.h>



namespace kinematics
{

  FORWARD(RootBondedJoint);
    class RootBondedJoint_O : public BondedJoint_O
    {
	LISP_CLASS(kinematics,KinPkg,RootBondedJoint_O,"RootBondedJoint",BondedJoint_O);
    public:
	static const NodeType nodeType = rootBondedJoint;
    protected:
	/*! Store the Id of the Bond1 atom if we have one
	  otherwise it will be UndefinedUnsignedInt */
	RootJointInfo	_RootInfo;
    public:
    RootBondedJoint_O(const chem::AtomId& atomId,core::T_sp name, const string& comment) :
	BondedJoint_O(atomId,name,comment) {};

      CL_DEFMETHOD double getPhi() const { return this->_Phi; };
      
	virtual RootJointInfo const * rootJointInfo() const { return &this->_RootInfo;};


	/*! Set the Bond1Id for this atom */

	void setup(core::Symbol_sp constitutionName,
		   core::Symbol_sp topologyName,
		   chem::Plug_sp inPlug);

	virtual core::Symbol_sp typeSymbol() const;

    };



};
#endif
