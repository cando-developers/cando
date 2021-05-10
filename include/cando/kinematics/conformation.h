/*
    File: conformation.h
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
#ifndef	_kinematicsConformation_H
#define _kinematicsConformation_H

#include <clasp/core/common.h>
#include <clasp/core/record.h>
#include <cando/chem/oligomer.fwd.h>
#include <clasp/core/symbolTable.h>
#include <cando/kinematics/kinematicsPackage.h>

namespace kinematics
{

  SYMBOL_SC_(KinPkg,resizeMoleculesEvent);
  SYMBOL_SC_(KinPkg,buildMoleculeUsingOligomerEvent);
  SYMBOL_EXPORT_SC_(KeywordPkg,newValue);
  SYMBOL_EXPORT_SC_(KeywordPkg,oldValue);
  SYMBOL_EXPORT_SC_(KeywordPkg,moleculeId);
  SYMBOL_EXPORT_SC_(KeywordPkg,residueId);

  FORWARD(FoldTree);
  FORWARD(JointTree);

  FORWARD(Conformation);
  class Conformation_O : public core::CxxObject_O
  {
    LISP_CLASS(kinematics,KinPkg,Conformation_O,"Conformation",core::CxxObject_O);
  public:
    static Conformation_sp make(core::List_sp oligomers);
    DEFAULT_CTOR_DTOR(Conformation_O);
  public:
    bool fieldsp() const { return true; };
    void fields(core::Record_sp node);
    void initialize();
  public:
	/*! Create the Conformation based on a single Oligomer */
    static Conformation_sp create(chem::Oligomer_sp oligomer);
  private:
#if 0
	/*! Add an Alchemist to the Conformation and return its MoleculeId */
    int addAlchemists(core::Cons_sp alchemists);
#endif
	/*! Build the JointTree that will build the Alchemists */
    void buildJointTree();
  public:	// Functions here
	/*! Return the number of molecules */
    int numberOfMolecules() const;
	/*! Resize the number of Molecules
	  @signal kin:resizeMoleculesEvent :data ( :newValue newNumberOfMolecules :oldValue oldNumberOfMolecules) */
    void resizeMolecules(int numMolecules);

	/*! Build a molecule as specified by the Oligomer.
	 It will use the first Monomer of each MultiMonomer as a starting Joint.
	@signal kin:rebuiltMolecule :data ( moleculeId ) */
    void buildMoleculeUsingOligomer(int moleculeId, chem::Oligomer_sp oligomer);

    CL_DEFMETHOD FoldTree_sp getFoldTree() const { return this->_FoldTree;};
    CL_DEFMETHOD JointTree_sp getJointTree() const { return this->_JointTree;};

  private:
	// instance variables here
	/*! Store the FoldTree */
    FoldTree_sp				_FoldTree;
	/*! Store the JointTree */
    JointTree_sp				_JointTree;
  };

}; /* kinematics */

TRANSLATE(kinematics::Conformation_O);

#endif /* _kinematicsConformation_H */
