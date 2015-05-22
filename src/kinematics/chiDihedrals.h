#ifndef	_kinematics_ChiDihedral_H
#define _kinematics_ChiDihedral_H

#include <clasp/core/common.h>
#include <cando/chem/constitution.fwd.h>
#include <clasp/core/multipleValues.fwd.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

FORWARD(ChiDihedral);
/*! Define a Chi-dihedral angle that is defined within a Rotamer
  Atom4 is built using a dihedral angle wrt Atom3-Atom2-Atom1.
  This internal coordinate needs to be defined/referenced within the AtomTemplate tree
  and the AtomTree.  Each atom index is the ConstitutionAtomIndex0N of the atom in the
  current residue.  A value of -1 refers to atoms in the previous residue.
 */
class ChiDihedral_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(kinematics,KinPkg,ChiDihedral_O,"ChiDihedral");
#if INIT_TO_FACTORIES
 public:
    static ChiDihedral_sp make(const int atom1Index, const int atom2Index, const int atom3Index, const int atom4Index);
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(ChiDihedral_O);
public:
	void initialize();

private: // instance variables here
	chem::ConstitutionIndex0N	_Atom1;
	chem::ConstitutionIndex0N	_Atom2;
	chem::ConstitutionIndex0N	_Atom3;
	chem::ConstitutionIndex0N	_Atom4;

public: // Functions here
	chem::ConstitutionIndex0N	atom1Index() const { return this->_Atom1;};
	chem::ConstitutionIndex0N	atom2Index() const { return this->_Atom2;};
	chem::ConstitutionIndex0N	atom3Index() const { return this->_Atom3;};
	chem::ConstitutionIndex0N	atom4Index() const { return this->_Atom4;};

	/*! Return all of the indices as integers
	  that represent the ConstitutionIndex0N for each atom
	  a negative value refers to an atom in a previous residue */
    core::T_mv atomIndices() const;
};




FORWARD(ChiList);
/*! Define a list of ChiDihedral angles that will be used to define
  a Rotamer.
 */
class ChiList_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(kinematics,KinPkg,ChiList_O,"ChiList");
#if INIT_TO_FACTORIES
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(ChiList_O);
public:
	void initialize();

private: // instance variables here
	core::Vector0<ChiDihedral_O>	_Chis;

public: // Functions here

	void append(ChiDihedral_sp chi);
};





}; /* kinematics */

TRANSLATE(kinematics::ChiDihedral_O);
TRANSLATE(kinematics::ChiList_O);

#endif /* _kin_ChiDihedral_H */


