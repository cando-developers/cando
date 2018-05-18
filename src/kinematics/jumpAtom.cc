/*
    File: jumpAtom.cc
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
#include <cando/kinematics/atomTree.h>
#include <clasp/core/symbolTable.h>
#include <cando/kinematics/stub.h>
#include <cando/kinematics/jumpAtom.h>

namespace kinematics
{

    void JumpJoint_O::_insertChild(int idx, Joint_sp c)
    {
	Joint_sp empty;
	this->_Children.insert(this->_Children.begin()+idx,empty);
	this->_Children[idx] = c;
    }


    void JumpJoint_O::_appendChild(Joint_sp c)
    {
	Joint_sp empty;
	int idx = this->_Children.size();
	this->_Children.push_back(empty);
	this->_Children[idx] = c;
    }


    void JumpJoint_O::_releaseChild(int idx)
    {
	this->_Children[idx].reset_();
	this->_Children.erase(this->_Children.begin()+idx);
    }


    void JumpJoint_O::_releaseAllChildren()
    {
	this->_Children.clear();
    }



    void JumpJoint_O::updateInternalCoords(Stub& stub,
					bool const recursive,
	AtomTree_sp at)
    {_OF();
	ASSERTF(stub.isOrthogonal(1e-3),BF("The incoming stub is not orthogonal"));
	Stub newStub(this->getStub(at));
	this->_Jump.fromStubs(stub,newStub);
	LOG(BF("Incoming stub:\n%s") % stub.asString() );
	LOG(BF("My new Stub:\n%s") % newStub.asString() );
	LOG(BF("Jump:\n%s") % this->_Jump.asString());
	if ( recursive )
	{
	    for (int childIdx=0; childIdx<this->_numberOfChildren(); childIdx++ )
	    {
		this->_child(childIdx).get()->updateInternalCoords(newStub,true,at);
	    }
	}
    }



    Joint_sp JumpJoint_O::stubAtom1() const
    {
      if (this->stubDefined()) {
        return this->asSmartPtr();
      }
      return this->parent();
    }

    Joint_sp JumpJoint_O::stubAtom2() const
    {
	if ( this->stubDefined() )
	{
	    return this->getNonJumpAtom(0);
	} else
	{
	    return this->parent()->stubAtom2();
	}
    }

    Joint_sp JumpJoint_O::stubAtom3(AtomTree_sp at) const
    {_OF();
	if ( this->stubDefined() )
	{
	    Joint_sp first(this->getNonJumpAtom(0));
	    ASSERT(first.boundp());
	    Joint_sp second(first->getNonJumpAtom(0));
	    if ( second.boundp() )
	    {
		return second;
	    } else
	    {
		return this->getNonJumpAtom(1);
	    }
	} else
	{
	    return this->parent()->stubAtom3(at);
	}
    }

    bool JumpJoint_O::keepDofFixed(DofType dof) const
    {
	return (this->parent().boundp());
    }



    /*! Update the external coordinates using the input stub */
    void JumpJoint_O::updateXyzCoords(Stub& stub,AtomTree_sp at)
    {_OF();
	ASSERTF(stub.isOrthogonal(1e-3),BF("Stub is not orthogonal - stub:\n%s") % stub.asString());
	Stub newStub;
	this->_Jump.makeJump(stub,newStub);
	this->position(newStub.translation());
	for ( int ii=0; ii<this->_numberOfChildren(); ii++ )
	{
	    this->_child(ii).get()->updateXyzCoords(newStub,at);
	}
	this->noteXyzUpToDate();
    }



    void JumpJoint_O::updateXyzCoords(AtomTree_sp at)
    {
	Stub stub(this->getInputStub(at));
	this->JumpJoint_O::updateXyzCoords(stub,at);
    }



    double JumpJoint_O::dof(DofType const& dof) const
    {_OF();
	if ( dof.isRigidBodyDof() )
	{
	    return _Jump.getRigidBodyDelta(dof,RigidBodyDirection::n2c);
	}
	SIMPLE_ERROR(BF("Illegal dof for JumpAtom - I can only return rigid body dofs"));
    }


    core::Symbol_sp JumpJoint_O::typeSymbol() const {_OF(); return _sym_jump;};



};
