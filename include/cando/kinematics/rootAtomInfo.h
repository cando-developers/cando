/*
    File: rootAtomInfo.h
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
#ifndef _kinematics_rootAtomInfo_H
#define _kinematics_rootAtomInfo_H

#include <clasp/core/common.h>
#include <cando/chem/plug.fwd.h>
#include <cando/chem/constitutionAtoms.fwd.h>

namespace kinematics
{

    /*! @class Stores information on the Residue that the containing Atom is the root of */

    class RootAtomInfo
    {
	friend class RootBondedAtom;
	friend class DelayedBondedAtom;
    private:
	core::Symbol_sp			_ConstitutionName;
	core::Symbol_sp			_TopologyName;
	chem::ConstitutionAtomIndex0N	_Bond1Id;
    public:
	void setup(core::Symbol_sp constitutionName,
		   core::Symbol_sp topologyName,
		   chem::Plug_sp inPlug);
    };


};

#endif // _kinematics_rootAtomInfo_H
