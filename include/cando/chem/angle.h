#ifndef	_chem_Angle_H
#define _chem_Angle_H

#include <clasp/core/common.h>
#include <cando/chem/chemPackage.h>

namespace chem
{

    FORWARD(Atom);
FORWARD(Angle);

/*! Angles store three atoms of an angle
    a1-a2-a3 where these atoms are sequentially bonded to each other.
*/

class Angle_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(chem,ChemPkg,Angle_O,"Angle");
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Angle_O);
public:
	void initialize();

private: // instance variables here
	Atom_sp	_a1;
	Atom_sp	_a2;
	Atom_sp	_a3;
public:
	static Angle_sp create(Atom_sp a1, Atom_sp a2, Atom_sp a3);
public: // Functions here
	void setup(Atom_sp a1, Atom_sp a2, Atom_sp a3);

	Atom_sp	getA1() const { return this->_a1;};
	Atom_sp	getA2() const { return this->_a2;};
	Atom_sp	getA3() const { return this->_a3;};

	string __repr__() const;
};

}; /* chem */

TRANSLATE(chem::Angle_O);

#endif /* _chem_Angle_H */


