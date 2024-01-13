/*
    File: chiDihedrals.h
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
  This internal coordinate needs to be defined/referenced within the JointTemplate tree
  and the JointTree.  Each atom index is the ConstitutionAtomIndex0N of the atom in the
  current residue.  A value of -1 refers to atoms in the previous residue.
 */
class ChiDihedral_O : public core::CxxObject_O
{
    LISP_CLASS(kinematics,KinPkg,ChiDihedral_O,"ChiDihedral",core::CxxObject_O);
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

public: // instance variables here
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
class ChiList_O : public core::CxxObject_O
{
    LISP_CLASS(kinematics,KinPkg,ChiList_O,"ChiList",core::CxxObject_O);
#if INIT_TO_FACTORIES
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(ChiList_O);
public:
	void initialize();

private: // instance variables here
	gctools::Vec0<ChiDihedral_sp>	_Chis;

public: // Functions here

	void append(ChiDihedral_sp chi);
};





}; /* kinematics */


#endif /* _kin_ChiDihedral_H */
