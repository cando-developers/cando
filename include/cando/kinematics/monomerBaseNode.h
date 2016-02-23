/*
    File: monomerBaseNode.h
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
#ifndef	_kinematicsMonomerBaseNode_H
#define _kinematicsMonomerBaseNode_H

#include <clasp/core/common.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/monomerBaseNode.fwd.h>

namespace kinematics
{


    class MonomerBaseNode_O : public core::T_O
{
    LISP_CLASS(kinematics,KinPkg,MonomerBaseNode_O,"MonomerBaseNode",core::T_O);
#if INIT_TO_FACTORIES
 public:
    static MonomerBaseNode_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(MonomerBaseNode_O);
public:
	void initialize();
public:
	// Functions here
private:
	// instance variables here
};

}; /* kinematics */

TRANSLATE(kinematics::MonomerBaseNode_O);

#endif /* _kinematicsMonomerBaseNode_H */
