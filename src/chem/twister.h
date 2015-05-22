       
       
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
#include "core/common.h"
//#include "matter.h"
//#include "coordSys.h"



#include "chemPackage.h"


namespace chem {

SMART(Atom);
SMART(CoordinateSystem);


SMART(Twister);
/*!
 * Twister class can be set up to rotate a dihedral in an absolute sense
 * or around a bond in a relative sense.
 */
class Twister_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,Twister_O,"Twister");
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
class TwisterDriver_O : public core::T_O
    {
	LISP_BASE1(core::T_O);
	LISP_CLASS(chem,ChemPkg,TwisterDriver_O,"TwisterDriver");
public:
//	void	archive(core::ArchiveP node);

private:
	Aggregate_sp		_Aggregate;
        gctools::Vec0<Twister_sp>		_Twisters;

public:

public:
    void	addTwister(Twister_sp twister);
    Twister_sp	getTwister(uint i);

    uint	numberOfTwisters() { return this->_Twisters.size(); };

    void perturbConformation();
    string description() const;


    DEFAULT_CTOR_DTOR(TwisterDriver_O);
    };




};

TRANSLATE(chem::Twister_O);
TRANSLATE(chem::TwisterDriver_O);
#endif

