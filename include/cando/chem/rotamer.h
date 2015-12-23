#ifndef	_chem_ROTAMER_H
#define _chem_ROTAMER_H

#include <clasp/core/common.h>
#include <cando/chem/chemPackage.h>

namespace chem
{



class RotamerAtom
{
    friend class Rotamer_O;
private:
//    void archive(core::ArchiveP node);
private:
    double	_Distance;
    double	_Angle;
    double	_Torsion;
};


class Rotamer_O : public core::CxxObject_O
{
    LISP_CLASS(chem,ChemPkg,Rotamer_O,"Rotamer",core::CxxObject_O);
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Rotamer_O);
public:
    void initialize();
	// Functions here
public:
    typedef std::vector<RotamerAtom>	RotamerAtoms;
private:
	// instance variables here
    /*! A list of distance/angle/dihedral values indexed by the ConstitutionAtomIndex0N */
    RotamerAtoms	_RotamerAtoms;
};

}; /* chem */

TRANSLATE(chem::Rotamer_O);

#endif /* _chem_ROTAMER_H */

