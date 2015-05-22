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

    class ImproperTorsion_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ImproperTorsion_O,"ImproperTorsion");
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
	static ImproperTorsion_sp create(Atom_sp a1, Atom_sp a2, Atom_sp central, Atom_sp a4,core::Lisp_sp lisp);
public: // Functions here

	void setup(Atom_sp a1, Atom_sp a2, Atom_sp a3, Atom_sp a4);


	Atom_sp	getA1() const { return this->_a1;};
	Atom_sp	getA2() const { return this->_a2;};
	Atom_sp	getA3Central() const { return this->_a3central;};
	Atom_sp	getA4() const { return this->_a4;};

	string __repr__() const;

};

}; /* chem */

TRANSLATE(chem::ImproperTorsion_O);

#endif /* _chem_ImproperTorsion_H */


