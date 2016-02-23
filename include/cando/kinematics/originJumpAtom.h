/*
    File: originJumpAtom.h
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
#ifndef	kinematics_originJumpAtom_H
#define kinematics_originJumpAtom_H

#include <clasp/core/foundation.h>
#include <cando/kinematics/atom.h>
#include <cando/kinematics/jump.h>
#include <cando/kinematics/jumpAtom.h>
#include <cando/chem/atomId.h>

namespace kinematics
{

    class OriginJumpAtom : public JumpAtom
    {
    public:
	static const NodeType nodeType = originJumpAtom;
    protected:
    public:
	/*! Empty ctor */
	OriginJumpAtom() {};

    OriginJumpAtom(const chem::AtomId& atomId,const string& comment) : JumpAtom(atomId,comment) {};



	virtual core::Symbol_sp typeSymbol() const;


	/*! The stubAtoms of an OriginJumpAtom are all myself */
	RefCountedAtomHandle stubAtom1() const { return this->_Me;};
	RefCountedAtomHandle stubAtom2() const { return this->_Me;};
	RefCountedAtomHandle stubAtom3(AtomTree_sp at) const { return this->_Me;};

	/*! Update the internal coordinates */
	virtual void updateInternalCoords(Stub& stub,
					  bool const recursive,
					  AtomTree_sp at);




    };



};
#endif
