#ifndef	constitutionAtoms_FWD_H
#define	constitutionAtoms_FWD_H


#include "chemPackage.h"


namespace chem
{

    /*!Atoms within Constitutions have a unique index numbered 0-(N-1) for rapid lookup in tables */
    typedef int ConstitutionAtomIndex0N;

    FORWARD(Atom);
//    typedef map<Atom_sp, ConstitutionAtomIndex0N> MapAtomsToConstitutionAtomIndex0N;
    typedef gctools::SmallMap<Atom_sp,int> MapAtomsToConstitutionAtomIndex0N;

    FORWARD(ConstitutionBond);
    FORWARD(ConstitutionAtom);
    FORWARD(ConstitutionVirtualAtom);
    FORWARD(ConstitutionAtoms);

};
#endif
