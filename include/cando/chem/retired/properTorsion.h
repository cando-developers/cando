#ifndef	_chem_ProperTorsion_H
#define _chem_ProperTorsion_H

#include <clasp/core/common.h>
#include <cando/chem/chemPackage.h>

namespace chem
{

    FORWARD(Atom);
FORWARD(ProperTorsion);

/*! ProperTorsions store four atoms of a proper dihedral
    a1-a2-a3-a4 where these atoms are sequentially bonded to each other.
*/

class ProperTorsion_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ProperTorsion_O,"ProperTorsion");
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(ProperTorsion_O);
public:
	void initialize();

private: // instance variables here
	Atom_sp	_a1;
	Atom_sp	_a2;
	Atom_sp	_a3;
	Atom_sp	_a4;
public:
	static ProperTorsion_sp create(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4,core::Lisp_sp lisp);
public: // Functions here

	void setup(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4);

	Atom_sp	getA1() const { return this->_a1;};
	Atom_sp	getA2() const { return this->_a2;};
	Atom_sp	getA3() const { return this->_a3;};
	Atom_sp	getA4() const { return this->_a4;};

	string __repr__() const;

};

}; /* chem */

TRANSLATE(chem::ProperTorsion_O);

#endif /* _chem_ProperTorsion_H */


