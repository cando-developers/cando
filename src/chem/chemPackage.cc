#include "core/common.h"
#include "core/package.h"
#include "chemPackage.h"
#include "candoScript.h"
//#include "candoDatabase.h"
#include "elements.h"
#include "mbbCoreTools.h"
#include "core/builtInClass.h"
//#define HEADER_INCLUDES
//#include "chem_initClasses_inc.h"
#include "ringFinder.h"
#include "largeSquareMatrix.h"
#include "alias.h"
#include "angle.h"
//#include "atomGrid.h"
#include "atomIdMap.h"
#include "atomId.h"
#include "atomIndexer.h"
#include "atomReference.h"
#include "energyAtomTable.h"
#include "bond.h"
#include "chemdraw.h"
#include "calculatePosition.h"
#include "candoDatabaseReference.h"
//#include "candoDatabase.h"
#include "chemInfo.h"
#include "cipPrioritizer.h"
#include "command.h"
#include "complexRestraints.h"
//#include "conformationCollection.h"
//#include "conformationExplorer.h"
#include "constitutionAtoms.h"
#include "coordSys.h"
#include "coupling.h"
#include "energyComponent.h"
#include "energyFunction.h"
#include "entity.h"
//#include "externalInterface.h"
#include "ffBaseDb.h"
#include "ffNonbondDb.h"
#include "forceField.h"
#include "frameRecognizer.h"
#include "improperTorsion.h"
#include "macroModel.h"
#include "matter.h"
#include "minimizerLog.h"
#include "minimizer.h"
#include "monomerContext.h"
#include "monomerCoordinates.h"
#include "nVector.h"
#include "numericalFunction.h"
#include "candoScript.h"
#include "oligomer.h"
#include "monomer.h"
#include "typeAssignmentRules.h"
#include "pdbMonomer.h"
#include "pdb.h"
#include "plug.h"
#include "properTorsion.h"
#include "readAmberParameters.h"
#include "restraint.h"
#include "rotamer.h"
#include "searchStatistics.h"
#include "spanningLoop.h"
#include "specificContext.h"
#include "statusTracker.h"
#include "stereochemistry.h"
#include "stereoisomerAtoms.h"
#include "structureComparer.h"
#include "structureList.h"
#include "superposeEngine.h"
#include "topology.h"
//#include "trainerArchive.h"
//#include "trainer.h"
//#include "jobHistory.h"
#include "trajectory.h"
#include "twister.h"
#include "virtualSphere.h"
#include "zMatrix.h"
#include "aggregate.h"
#include "atom.h"
#include "constitution.h"
#include "energyAnchorRestraint.h"
#include "energyAngle.h"
#include "energyChiralRestraint.h"
#include "energyDihedral.h"
#include "energyFixedNonbond.h"
#include "energyImproperRestraint.h"
#include "energyNonbond.h"
#include "energyStretch.h"
#include "entityNameSet.h"
#include "ffStretchDb.h"
#include "ffAngleDb.h"
#include "ffItorDb.h"
#include "ffPtorDb.h"
#include "ffTypesDb.h"
#include "ffVdwDb.h"
#include "iterateMatter.h"
#include "iterateRestraints.h"
#include "molecule.h"
//#include "grPickableMatter.h"
#include "residue.h"
#include "superposableConformationCollection.h"
#include "virtualAtom.h"
#include "monomerPack.h"
#include "representedEntityNameSet.h"

namespace chemkw {


#pragma GCC visibility push(default)
#define ChemKwPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include "symbols_scraped_inc.h"
#undef DO_SYMBOL
#undef ChemKwPkg_SYMBOLS
#pragma GCC visibility pop

};

namespace chem {


#pragma GCC visibility push(default)
#define ChemPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include "symbols_scraped_inc.h"
#undef DO_SYMBOL
#undef ChemPkg_SYMBOLS
#pragma GCC visibility pop

};



namespace chem
{
#define EXPOSE_TO_CANDO
#define Use_ChemPkg
#define EXTERN_REGISTER
#include "chem_initClasses_inc.h"
#undef EXTERN_REGISTER
#undef Use_ChemPkg
#undef EXPOSE_TO_CANDO
};




namespace chem
{


    const char* Chem_nicknames[] = {""};

//! The default CandoDatabase
    SYMBOL_SC_(ChemPkg,candoDatabase);
    SYMBOL_SC_(ChemPkg,AM1_BCC_ar5);
    SYMBOL_SC_(ChemPkg,AM1_BCC_ar6);
    SYMBOL_SC_(ChemPkg,AM1_BCC_ar7);

    extern void	Initialize_Mol2_TypeRules(core::Lisp_sp lisp);



    void ChemExposer::expose(core::Lisp_sp lisp, WhatToExpose what) const
    {_G();
	switch (what)
	{
	case candoClasses:
	{
#define ChemKwPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {chemkw::cname = _lisp->internWithPackageName(pkg,lispname); chemkw::cname->exportYourself(exportp);}
#include "symbols_scraped_inc.h"
#undef DO_SYMBOL
#undef ChemKwPkg_SYMBOLS


#define ChemPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {cname = _lisp->internWithPackageName(pkg,lispname); cname->exportYourself(exportp);}
#include "symbols_scraped_inc.h"
#undef DO_SYMBOL
#undef ChemPkg_SYMBOLS





#define ALL_STAGES
#define	Use_ChemPkg
#define INVOKE_REGISTER
#define LOOKUP_SYMBOL(pkg,name) _lisp->internWithPackageName(pkg,name)
#include "chem_initClasses_inc.h"
#undef INVOKE_REGISTER
#undef Use_ChemPkg
#undef ALL_STAGES
	}
	break;
	case candoFunctions:
	{
	    setupCandoPrimitives(_lisp);
#if 0
	    defNoWrapPackage(ChemPkg,"atomAliases", &core::prim_list ,_lisp->lisp());
	    defNoWrapPackage(ChemPkg,"aliasAtoms", &core::prim_list ,_lisp->lisp());
	    defNoWrapPackage(ChemPkg,"addPart", &core::prim_list ,_lisp->lisp());
	    defNoWrapPackage(ChemPkg,"parts", &core::prim_list ,_lisp->lisp());
	    defNoWrapPackage(ChemPkg,"monomerAlias", &core::prim_list ,_lisp->lisp());
#endif
	    exposeMbbCoreTools(_lisp);
	}
	break;
	case candoGlobals:
	{
	    Initialize_Mol2_TypeRules(_lisp);
	    initializeElementsAndHybridization();
//	    _lisp->defvar(_sym_candoDatabase,cdb);
	}
	break;
	case pythonClasses:
	{
#define _DBG(x)
#define EXPOSE_TO_PYTHON
#define	Use_ChemPkg
#include "chem_initClasses_inc.h"
#undef Use_ChemPkg
#undef EXPOSE_TO_PYTHON
#undef _DBG
	}
	break;
	case pythonFunctions:
	{
	    setupPythonPrimitives(_lisp);
	    // nothing
	}
	break;
	case pythonGlobals:
	{
	    // nothing currently
	}
	break;
	}
    }


#if 0
    CandoDatabase_sp getCandoDatabase(const core::Lisp_sp& lisp)
    {_G();
	return _sym_candoDatabase->symbolValue().as<CandoDatabase_O>();
    }

    /*! Intern a symbol in the ChemKwPkg and return it.
      The symbol is automatically exported */
    MatterName chemkw_intern(const string& name) {
        core::Package_sp chemkwPkg = _lisp->findPackage(ChemKwPkg);
        string upperName = name;
        for ( char* cP = upperName.c_str(); *cP; ++cP ) {
            *cP = toupper(*cp);
        }
        core::Symbol_sp sym = chemkwPkg->intern(upperName);
        chemkwPkg->export(sym);
        return sym;
    }
#endif
        

};




#if USE_INTRUSIVE_SMART_PTR==1
#define EXPAND_CLASS_MACROS
#define _CLASS_MACRO(_T_)				\
    STATIC_CLASS_INFO(_T_);			\
    INTRUSIVE_POINTER_REFERENCE_COUNT_ACCESSORS(_T_);
#include "chem_initClasses_inc.h"
#undef _CLASS_MACRO
#undef EXPAND_CLASS_MACROS
#endif


extern "C"
{
//    bool init_libchem()
    void ___user_chem()
    {
	static chem::ChemExposer* chemPkgP = new chem::ChemExposer(_lisp);
	_lisp->installPackage(chemPkgP);
    }

}



 
