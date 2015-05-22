#ifndef	_kinematicsConformation_H
#define _kinematicsConformation_H

#include <clasp/core/common.h>
#include <core/model.h>
#include <cando/chem/oligomer.fwd.h>
#include <cando/kinematics/symbolTable.h>
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
    FORWARD(AtomTree);

    FORWARD(Conformation);
    class Conformation_O : public core::Model_O
{
    LISP_BASE1(core::Model_O);
    LISP_CLASS(kinematics,KinPkg,Conformation_O,"Conformation");
#if INIT_TO_FACTORIES
 public:
    static Conformation_sp make(const core::Cons_sp& oligomers);
#else
    DECLARE_INIT();
#endif
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Conformation_O);
public:
	void initialize();
	static void lisp_initGlobals(core::Lisp_sp lisp);
public:
	/*! Create the Conformation based on a single Oligomer */
	static Conformation_sp create(chem::Oligomer_sp oligomer);

private:
#if 0
	/*! Add an Alchemist to the Conformation and return its MoleculeId */
	int addAlchemists(core::Cons_sp alchemists);
#endif
	/*! Build the AtomTree that will build the Alchemists */
	void buildAtomTree();
public:	// Functions here
	/*! Return the number of molecules */
	int numberOfMolecules() const;
	/*! Resize the number of Molecules
	  @signal kin:resizeMoleculesEvent :data ( :newValue newNumberOfMolecules :oldValue oldNumberOfMolecules) */
	void resizeMolecules(int numMolecules);

	/*! Build a molecule as specified by the Oligomer.
	 It will use the first Monomer of each MultiMonomer as a starting point.
	@signal kin:rebuiltMolecule :data ( moleculeId ) */
	void buildMoleculeUsingOligomer(int moleculeId, chem::Oligomer_sp oligomer);

	FoldTree_sp getFoldTree() const { return this->_FoldTree;};
	AtomTree_sp getAtomTree() const { return this->_AtomTree;};

private:
	// instance variables here
	/*! Store the FoldTree */
	FoldTree_sp				_FoldTree;
	/*! Store the AtomTree */
	AtomTree_sp				_AtomTree;
};

}; /* kinematics */

TRANSLATE(kinematics::Conformation_O);

#endif /* _kinematicsConformation_H */


