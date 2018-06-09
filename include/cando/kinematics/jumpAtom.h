/*
    File: jumpAtom.h
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
#include <cando/kinematics/atom.h>
#include <cando/kinematics/jump.h>
#include <cando/chem/atomId.h>

namespace kinematics
{

  FORWARD(JumpJoint);
  class JumpJoint_O : public Joint_O {
    LISP_CLASS(kinematics,KinPkg,JumpJoint_O,"JumpAtom",Joint_O);
  public:
    static const NodeType nodeType = jumpAtom;

  protected:
    Jump			_Jump;
	/*! JumpAtoms can have unlimited numbers of children */
    gc::Vec0< Joint_sp >	_Children;
  protected:
	/*! Bonded atoms can have different numbers of children wrt JumpAtoms */
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
	/*! Empty ctor */
    JumpJoint_O() {};

  JumpJoint_O(const chem::AtomId& atomId, core::T_sp name, const string& comment) : Joint_O(atomId,name,comment) {};



    virtual core::Symbol_sp typeSymbol() const;


	/*! Update the internal coordinates */
    virtual void updateInternalCoords(Stub& stub,
                                      bool const recursive,
                                      AtomTree_sp at);


	/*! Yes, this is a JumpAtom */
    bool isJump() const { return true;};

	/*! Return the stubAtom1 */
    virtual Joint_sp stubAtom1() const;

	/*! Return the stubAtom2 */
    virtual Joint_sp stubAtom2() const;

	/*! Return the stubAtom3 */
    virtual Joint_sp stubAtom3(AtomTree_sp at) const;

    bool keepDofFixed(DofType dof) const;


	/*! Update the external coordinates */
    virtual void updateXyzCoords(AtomTree_sp at);


	/*! Update the external coordinates using the input stub */
    virtual void updateXyzCoords(Stub& stub,AtomTree_sp at);


	/*! Return the DOF */
    double dof(DofType const& dof) const;
	
  };



};
#endif
