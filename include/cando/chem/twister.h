/*
    File: twister.h
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
       
       
//
// (C) 2004 Christian E. Schafmeister
//




#ifndef	TWISTER_H
#define	TWISTER_H

#include <stdio.h>
#include <math.h>
#include <string>
#include <map>
#include <exception>
#include <clasp/core/common.h>
//#include "matter.h"
//#include "coordSys.h"



#include <cando/chem/chemPackage.h>


namespace chem {

SMART(Atom);
SMART(CoordinateSystem);


SMART(Twister);
/*!
 * Twister class can be set up to rotate a dihedral in an absolute sense
 * or around a bond in a relative sense.
 */
class Twister_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Twister_O,"Twister",core::CxxObject_O);
public:
	void initialize();
public:
//	void	archive(core::ArchiveP node);
private:
	Atom_sp			_FixedRef;
	Atom_sp			_Fixed;
	Atom_sp			_Movable;
	Atom_sp			_MovableRef;
    gctools::Vec0<Atom_sp>		_Atoms;	//!< Atoms to rotate

public:

private:
	void	_defineForDihedral(Atom_sp a1ref, Atom_sp a1, Atom_sp a2, Atom_sp a2ref, bool force2mobile);
public:


	void	setFixedAndMovable(Atom_sp fixed, Atom_sp movable);

	void	addAtom(Atom_sp atom);

		/*! Rotate the twister by angle(radians) relative
		 */
	void	rotate(double angle);

		/*! Rotate the twister to angle(radians) absolute
		 */
	void	rotateAbsolute(double angle);

		/*! Define a twister between a1-a2
		 */
	void	defineForBond(Atom_sp a1, Atom_sp a2);

		/*! Define a twister across a1-a2
		 * Does exactly what defineForBond but works
		 * when there isn't a bond between the atoms
		 */
	void	defineFixedAndMobile(Atom_sp fixed, Atom_sp mobile);


		/*! Define a twister between a1-a2
		 */
	void	defineForDihedral(Atom_sp a1ref, Atom_sp a1, Atom_sp a2, Atom_sp a2ref);


	string	description() const;


	DEFAULT_CTOR_DTOR(Twister_O);
};




SMART(TwisterDriver);
class TwisterDriver_O : public core::CxxObject_O
    {
	LISP_CLASS(chem,ChemPkg,TwisterDriver_O,"TwisterDriver",core::CxxObject_O);
public:
//	void	archive(core::ArchiveP node);

private:
	Aggregate_sp		_Aggregate;
        gctools::Vec0<Twister_sp>		_Twisters;

public:

public:
    void	addTwister(Twister_sp twister);
    Twister_sp	getTwister(uint i);

CL_NAME("numberOfTwisters");
CL_DEFMETHOD     uint	numberOfTwisters() { return this->_Twisters.size(); };

    void perturbConformation();
    string description() const;


    DEFAULT_CTOR_DTOR(TwisterDriver_O);
    };




};

TRANSLATE(chem::Twister_O);
TRANSLATE(chem::TwisterDriver_O);
#endif
