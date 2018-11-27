/*
    File: virtualAtom.h
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
#ifndef	VirtualAtom_H //[
#define VirtualAtom_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/atom.h>
#include <cando/chem/chemPackage.h>


namespace chem {


SMART(CalculatePosition);

SMART(VirtualAtom );
class VirtualAtom_O : public Atom_O
{
    LISP_CLASS(chem,ChemPkg,VirtualAtom_O,"VirtualAtom",Atom_O);
public:
	void	initialize();
        bool fieldsp() const {return true;};
        void fields(core::Record_sp node);
private:
		//! Calculate position using a procedure
	CalculatePosition_sp	_CalculatePosition;
public:
    static VirtualAtom_sp create(MatterName name, CalculatePosition_sp proc);
public:

	void calculatePosition();

	virtual bool	equal(core::T_sp obj) const;

	virtual bool isVirtual() const { return true;};


			// Constructor stuff
	virtual Matter_sp copy();
//	virtual Atom_sp	copyDropExternalResidueBonds();


	VirtualAtom_O( const VirtualAtom_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(VirtualAtom_O);
};


};

TRANSLATE(chem::VirtualAtom_O);

#endif //]
