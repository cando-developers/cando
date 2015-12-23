#ifndef	_BONDID_H
#define _BONDID_H

#include <clasp/core/common.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/bondId.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>



namespace kinematics
{



    class BondId_O : public core::T_O
{
    friend class BondedAtomTemplate_O;
    friend class RootBondedAtomTemplate_O;

    LISP_CLASS(kinematics,KinPkg,BondId_O,"BondId",core::T_O);
#if INIT_TO_FACTORIES
 public:
    static BondId_sp make();
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(BondId_O);
private: // instance variables here
    Atom_sp	_Parent;
    Atom_sp	_Child;

public: // Functions here
    void initialize();

public:
    static BondId_sp create(Atom_sp parent, Atom_sp child);
};

}; /* kinematics */

TRANSLATE(kinematics::BondId_O);

#endif /* _BONDID_H */

