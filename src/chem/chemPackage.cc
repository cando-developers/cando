#include <clasp/core/common.h>
#include <clasp/core/package.h>
#include <cando/chem/chemPackage.h>
#include <cando/chem/candoScript.h>
//#include "candoDatabase.h"
#include <cando/chem/elements.h>
#include <cando/chem/mbbCoreTools.h>
#include <clasp/core/builtInClass.h>
//#define HEADER_INCLUDES
//#include "chem_initClasses_inc.h"
#include <cando/chem/ringFinder.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/alias.h>
#include <cando/chem/angle.h>
//#include "atomGrid.h"
#include <cando/chem/atomIdMap.h>
#include <cando/chem/atomId.h>
#include <cando/chem/atomIndexer.h>
#include <cando/chem/atomReference.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/bond.h>
#include <cando/chem/chemdraw.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/candoDatabaseReference.h>
//#include "candoDatabase.h"
#include <cando/chem/chemInfo.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/command.h>
#include <cando/chem/complexRestraints.h>
//#include "conformationCollection.h"
//#include "conformationExplorer.h"
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/coordSys.h>
#include <cando/chem/coupling.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/entity.h>
//#include "externalInterface.h"
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/forceField.h>
#include <cando/chem/frameRecognizer.h>
#include <cando/chem/improperTorsion.h>
#include <cando/chem/macroModel.h>
#include <cando/chem/matter.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/minimizer.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/monomerCoordinates.h>
#include <cando/chem/nVector.h>
#include <cando/chem/numericalFunction.h>
#include <cando/chem/candoScript.h>
#include <cando/chem/oligomer.h>
#include <cando/chem/monomer.h>
#include <cando/chem/typeAssignmentRules.h>
#include <cando/chem/pdbMonomer.h>
#include <cando/chem/pdb.h>
#include <cando/chem/plug.h>
#include <cando/chem/properTorsion.h>
#include <cando/chem/readAmberParameters.h>
#include <cando/chem/restraint.h>
#include <cando/chem/rotamer.h>
#include <cando/chem/searchStatistics.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/specificContext.h>
#include <cando/chem/statusTracker.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <cando/chem/structureComparer.h>
#include <cando/chem/structureList.h>
#include <cando/chem/superposeEngine.h>
#include <cando/chem/topology.h>
//#include "trainerArchive.h"
//#include "trainer.h"
//#include "jobHistory.h"
#include <cando/chem/trajectory.h>
#include <cando/chem/twister.h>
#include <cando/chem/virtualSphere.h>
#include <cando/chem/zMatrix.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/atom.h>
#include <cando/chem/constitution.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyFixedNonbond.h>
#include <cando/chem/energyImproperRestraint.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyStretch.h>
#include <cando/chem/entityNameSet.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffVdwDb.h>
#include <cando/chem/iterateMatter.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/chem/molecule.h>
//#include "grPickableMatter.h"
#include <cando/chem/residue.h>
#include <cando/chem/superposableConformationCollection.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/monomerPack.h>
#include <cando/chem/representedEntityNameSet.h>

namespace chemkw {


#pragma GCC visibility push(default)
#define ChemKwPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef ChemKwPkg_SYMBOLS
#pragma GCC visibility pop

};

namespace chem {


#pragma GCC visibility push(default)
#define ChemPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef ChemPkg_SYMBOLS
#pragma GCC visibility pop

};



namespace chem
{
#define EXPOSE_TO_CANDO
#define Use_ChemPkg
#define EXTERN_REGISTER
#include <chem_initClasses_inc.h>
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
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef ChemKwPkg_SYMBOLS


#define ChemPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {cname = _lisp->internWithPackageName(pkg,lispname); cname->exportYourself(exportp);}
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef ChemPkg_SYMBOLS





#define ALL_STAGES
#define	Use_ChemPkg
#define INVOKE_REGISTER
#define LOOKUP_SYMBOL(pkg,name) _lisp->internWithPackageName(pkg,name)
#include <chem_initClasses_inc.h>
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
#include <chem_initClasses_inc.h>
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
#include <chem_initClasses_inc.h>
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



 
