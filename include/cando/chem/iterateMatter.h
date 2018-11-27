/*
    File: iterateMatter.h
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
#ifndef	O_Loop_H //[
#define O_Loop_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/loop.h>
#include <clasp/core/iterator.h>


#include <cando/chem/chemPackage.h>


namespace chem 
{




SMART(IterateMatter);
class IterateMatter_O : public core::Iterator_O
{
    LISP_CLASS(chem,ChemPkg,IterateMatter_O,"IterateMatter",core::Iterator_O);
#if INIT_TO_FACTORIES
 public:
    static IterateMatter_sp make();
#else
    DECLARE_INIT();
#endif
    public:
static void lisp_initGlobals(core::Lisp_sp lisp);
public:
	void	initialize();
public:
	Matter_sp _Top;
	int	_Goal;
	Loop	_Loop;
	bool	_IsDone;
public:
    static IterateMatter_sp create(Matter_sp top, int goal );

public:

		//! Initialize the top of the loop and the goal
	void initTopAndGoal(Matter_sp top, int goal );
	void first();
	void next();
	core::T_sp currentObject();
	bool isDone();
	bool advance();
	Atom_sp getAtom();
	Atom_sp getBondA1();
	Atom_sp getBondA2();
	int	getBondOrder();

	IterateMatter_O( const IterateMatter_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateMatter_O);
};




SMART(IterateAtoms);
class IterateAtoms_O : public IterateMatter_O
{
    LISP_CLASS(chem,ChemPkg,IterateAtoms_O,"IterateAtoms",IterateMatter_O);
#if INIT_TO_FACTORIES
 public:
    static IterateAtoms_sp make(Matter_sp matter);
#else
    DECLARE_INIT();
#endif
public:
public:
    static IterateAtoms_sp create(Matter_sp top);
    static core::T_sp createForMatter(Matter_sp matter);
public:

    core::T_sp currentObject();

	IterateAtoms_O( const IterateAtoms_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateAtoms_O);
};




SMART(IterateResidues);
class IterateResidues_O : public IterateMatter_O
{
    LISP_CLASS(chem,ChemPkg,IterateResidues_O,"IterateResidues",IterateMatter_O);
public:
private:
public:
    static IterateResidues_sp create(Matter_sp top);
    static core::T_sp createForMatter(Matter_sp matter);

public:

    core::T_sp currentObject();


	DEFAULT_CTOR_DTOR(IterateResidues_O);
};



SMART(IterateBonds);
class IterateBonds_O : public IterateMatter_O
{
    LISP_CLASS(chem,ChemPkg,IterateBonds_O,"IterateBonds",IterateMatter_O);
public:
private:
public:
    static IterateBonds_sp make(Matter_sp top);
    static core::T_sp createForMatter(Matter_sp matter);

public:

    core::T_sp currentObject();
	Atom_sp	getAtom1();
	Atom_sp getAtom2();
	int getBondOrder();



	DEFAULT_CTOR_DTOR(IterateBonds_O);
};




};
TRANSLATE(chem::IterateMatter_O);
TRANSLATE(chem::IterateAtoms_O);
TRANSLATE(chem::IterateResidues_O);
TRANSLATE(chem::IterateBonds_O);
#endif //]
