/*
    File: atomId.cc
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

#include <clasp/core/common.h>
#include <cando/chem/atomId.h>
#include <clasp/core/wrappers.h>



namespace chem
{

    string AtomId::asString() const
    {
	stringstream ss;
	ss << "Mol["<<this->_Molecule;
	ss << "]Res["<<this->_Residue;
	ss << "]Atom["<<this->_Atom<<"]";
	return ss.str();
    }






#if 0
CL_LISPIFY_NAME(make-atom-id);
DOCGROUP(cando);
CL_DEFUN AtomId_sp AtomId_O::make(int mol, int res, int atom)
    {
      auto  ai  = gctools::GC<AtomId_O>::allocate();
	ai->_AtomId._Molecule = mol;
	ai->_AtomId._Residue = res;
	ai->_AtomId._Atom = atom;
	return ai;
    }

#endif


};
