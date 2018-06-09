/*
    File: bondedAtom.h
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
#include <cando/kinematics/atom.h>



namespace kinematics
{

    class BondedJoint_O : public Joint_O
    {
	LISP_CLASS(kinematics,KinPkg,BondedJoint_O,"BondedAtom",Joint_O);
    public:
        bool fieldsp() { return true; };
        void fields(core::Record_sp node); 
        void initialize();
    public:
	static const NodeType nodeType = bondedAtom;
	static const int MaxChildren = 5;
    protected:
	int		_NumberOfChildren;
        // _Children have the value 0x0 if unbound
	Joint_sp	_Children[MaxChildren];
	Real		_Phi;
	Real		_Theta;
	Real		_Distance;
	bool		_DofChangePropagatesToYoungerSiblings;
    public:
	/*! Bonded atoms can have different numbers of children wrt JumpAtoms */
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
	virtual void _appendChild(Joint_sp c) {this->_Children[this->_NumberOfChildren++] = c; };
	/*! Delete all of the children for the destructor */
	virtual void _releaseAllChildren();

    public:
    BondedJoint_O() : Joint_O(), _NumberOfChildren(0) {};
    BondedJoint_O(const chem::AtomId& atomId, core::T_sp name, const string& comment) : Joint_O(atomId,name,comment), _NumberOfChildren(0) {};

	virtual core::Symbol_sp typeSymbol() const;

	/*! Return the stubAtom1 */
	virtual Joint_sp stubAtom1() const { return this->asSmartPtr();}

	/*! Return the stubAtom2 */
	virtual Joint_sp stubAtom2() const { return this->parent();};

	/*! Return the stubAtom3 */
	virtual Joint_sp stubAtom3(AtomTree_sp tree) const;

	/*! Update the internal coordinates */
	virtual void updateInternalCoords(Stub& stub,
					  bool const recursive,
					  AtomTree_sp at);


	bool keepDofFixed(DofType dof,AtomTree_sp at) const;


	string asString() const;


	/*! Update the external coordinates */
	virtual void updateXyzCoords(AtomTree_sp at);


	/*! Update the external coordinates using the input stub */
	virtual void updateXyzCoords(Stub& stub,AtomTree_sp at);



	/*! Get the value of the DOF */
	double dof(DofType const& dof) const;


    };



};
#endif
