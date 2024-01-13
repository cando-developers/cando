/*
    File: improperTorsion.h
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
#ifndef	_chem_ImproperTorsion_H
#define _chem_ImproperTorsion_H

#include <clasp/core/common.h>
#include <cando/chem/chemPackage.h>

namespace chem
{

    FORWARD(Atom);
    FORWARD(ImproperTorsion);

/*! ImproperTorsions store four atoms of an improper dihedral
    a1-a2-a3central-a4 where (central) is the central atom
    and a1,a2,a4 are all bonded to a3central and a1 and a2 are not
    bonded to each other.  This is the Amber way of defining impropers.
*/

    class ImproperTorsion_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,ImproperTorsion_O,"ImproperTorsion",core::CxxObject_O);
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(ImproperTorsion_O);
public:
	void initialize();

private: // instance variables here
	Atom_sp	_a1;
	Atom_sp	_a2;
	Atom_sp	_a3central;
	Atom_sp	_a4;
public:
	static ImproperTorsion_sp create(Atom_sp a1, Atom_sp a2, Atom_sp central, Atom_sp a4);
public: // Functions here

	void setup(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4);


	Atom_sp	getA1() const { return this->_a1;};
	Atom_sp	getA2() const { return this->_a2;};
	Atom_sp	getA3Central() const { return this->_a3central;};
	Atom_sp	getA4() const { return this->_a4;};

	string __repr__() const;

};

}; /* chem */


#endif /* _chem_ImproperTorsion_H */
