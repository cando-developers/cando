/*
    File: candoScript.h
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
       
       
#ifndef	candoScript_H //[
#define candoScript_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>


#include <cando/chem/chemPackage.h>
#include <cando/geom/coordinateArray.h>

namespace chem {

extern	void	setupCandoPrimitives(core::LispPtr);
extern	void	setupPythonPrimitives(core::LispPtr);

 FORWARD(Monomer);
SMART(CandoDatabase);
SMART(OligomerPart_Base);
/*! A simple object to store and create a monomer within an oligomer
 */
SMART(OligomerPart_Base);
class OligomerPart_Base_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,OligomerPart_Base_O,"OligomerPart-Base",core::CxxObject_O);

public:
public:
	virtual Monomer_sp createMonomer(CandoDatabase_sp bdb) {_OF();SUBCLASS_MUST_IMPLEMENT();};
	DEFAULT_CTOR_DTOR(OligomerPart_Base_O);
};


SMART(OligomerPart_Monomer);
/*! A simple object to store and create a monomer within an oligomer
 */
SMART(OligomerPart_Monomer);
class OligomerPart_Monomer_O : public OligomerPart_Base_O
{
    LISP_CLASS(chem,ChemPkg,OligomerPart_Monomer_O,"OligomerPart-Monomer",OligomerPart_Base_O);

public:
public:	// Leave everything exposed to simplify getting values in and out
	core::Symbol_sp _MonomerId;
	core::Symbol_sp	_GroupName;
	core::List_sp	_MonomerAliases;
	string	_Comment;
public:
public:
    virtual Monomer_sp createMonomer(CandoDatabase_sp bdb);
    DEFAULT_CTOR_DTOR(OligomerPart_Monomer_O);
};


SMART(OligomerPart_Link);
/*! A simple object to store and create a monomer within an oligomer
 */
SMART(OligomerPart_Link);
class OligomerPart_Link_O : public OligomerPart_Base_O
{
    LISP_CLASS(chem,ChemPkg,OligomerPart_Link_O,"OligomerPart-Link",OligomerPart_Base_O);

public:
public:	// Leave everything exposed to simplify getting values in and out
	core::Symbol_sp		_Monomer1Id;
	core::Symbol_sp		_Coupling;
	OligomerPart_Monomer_sp	_Monomer2;
public:
public:
    virtual Monomer_sp createMonomer(CandoDatabase_sp bdb);
    DEFAULT_CTOR_DTOR(OligomerPart_Link_O);
};


extern void setCandoDatabase(CandoDatabase_sp bdb);


geom::SimpleVectorCoordinate_sp chem__make_simple_vector_coordinate_from_atom_list(core::List_sp atoms);






};
TRANSLATE(chem::OligomerPart_Base_O);
TRANSLATE(chem::OligomerPart_Monomer_O);
TRANSLATE(chem::OligomerPart_Link_O);
#endif //]
