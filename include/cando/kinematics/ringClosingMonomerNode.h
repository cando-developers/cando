/*
    File: ringClosingMonomerNode.h
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
#ifndef	_kinematicsRingClosingMonomerNode_H
#define _kinematicsRingClosingMonomerNode_H

#include <clasp/core/common.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

    FORWARD(RingClosingMonomerNode);
class RingClosingMonomerNode_O : public MonomerNode_O
{
    LISP_CLASS(kinematics,KinPkg,RingClosingMonomerNode_O,"RingClosingMonomerNode",MonomerNode_O);
#if INIT_TO_FACTORIES
 public:
    static RingClosingMonomerNode_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(RingClosingMonomerNode_O);
public:
	void initialize();
public:
	static RingClosingMonomerNode_sp create(int monomerId);
private:
	/*! Store weak pointers to the other ring closing monomers indexed
	  by the name of the ring closing plug */
	adapt::SymbolMap<RingClosingMonomerNode_O>	_RingClosingChildren;
public:
	// Functions here

};

}; /* kinematics */

TRANSLATE(kinematics::RingClosingMonomerNode_O);

#endif /* _kinematicsRingClosingMonomerNode_H */
