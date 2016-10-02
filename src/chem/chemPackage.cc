/*
    File: chemPackage.cc
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
#include <clasp/core/common.h>
#include <clasp/core/package.h>
#include <clasp/core/pathname.h>
#include <clasp/core/compiler.h>
#include <cando/chem/chemPackage.h>
#include <cando/chem/candoScript.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/elements.h>
#include <cando/adapt/indexedObjectBag.h>
#include <cando/chem/mbbCoreTools.h>
#include <clasp/core/builtInClass.h>
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
//#include <cando/chem/candoDatabaseReference.h>
//#include "candoDatabase.h"
#include <cando/chem/chemInfo.h>
#include <cando/chem/loop.fwd.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/command.h>
#include <cando/chem/complexRestraints.h>
//#include "conformationCollection.h"
//#include "conformationExplorer.h"
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/coordSys.h>
#include <cando/chem/coupling.h>
#include <cando/chem/mol2.h>
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
#include <cando/chem/chimera.h>
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

SYMBOL_SHADOW_EXPORT_SC_(ChemPkg,atom);
SYMBOL_EXPORT_SC_(ChemPkg,R);
SYMBOL_EXPORT_SC_(ChemPkg,S);
SYMBOL_EXPORT_SC_(ChemPkg,UnknownConfiguration);
namespace chemkw {



};

namespace chem {
};




namespace chem
{
#if 0
#define EXPOSE_TO_CANDO
#define Use_ChemPkg
#define EXTERN_REGISTER
#include <clasp/core/initClasses.h>
#undef EXTERN_REGISTER
#undef Use_ChemPkg
#undef EXPOSE_TO_CANDO
#endif
};

#if 0
void chem_package_initializer()
{
//  printf("%s:%d Running chem_packageinitializer\n", __FILE__, __LINE__ );
                    // Create the CHEM-KW package
  std::list<std::string> nicknames = { "CK" };
  std::list<std::string> usePackages = { };
          // core::Package_sp chemkwpkg = _lisp->makePackage(ChemKwPkg,nicknames,usePackages);
          // chemkwpkg->setActsLikeKeywordPackage(true);

  core::T_sp pn;
  char* env = getenv("CANDO_LISP_SOURCE_DIR");
  if ( env != NULL ) {
    pn = core::cl__translate_logical_pathname(core::Str_O::create(std::string(env)+"/**/*.*"));
  } else {
    pn = core::cl__translate_logical_pathname(core::Str_O::create("APP-RESOURCES:CANDO;**;*.*"));
  }
  core::Cons_sp pts = core::Cons_O::createList(core::Str_O::create("cando:**;*.*"),pn);
  core::Cons_sp ptsList = core::Cons_O::createList(pts);
  core::core__pathname_translations(core::Str_O::create("CANDO"),_lisp->_true(),ptsList);
  chem::energyFunction_initializeSmarts();
  chem::initialize_mol2_type_rules();
  chem::setupCandoPrimitives(_lisp);
  chem::initialize_chimera();
  chem::initializeElementsAndHybridization();
}

core::Initializer global_chem_package_initializer(chem_package_initializer);
#endif




namespace chem
{


    const char* Chem_nicknames[] = {""};

//! The default CandoDatabase
    SYMBOL_EXPORT_SC_(ChemPkg,candoDatabase);
    SYMBOL_EXPORT_SC_(ChemPkg,AM1_BCC_ar5);
    SYMBOL_EXPORT_SC_(ChemPkg,AM1_BCC_ar6);
    SYMBOL_EXPORT_SC_(ChemPkg,AM1_BCC_ar7);

    extern void	Initialize_Mol2_TypeRules(core::Lisp_sp lisp);



    void ChemExposer_O::expose(core::Lisp_sp lisp, WhatToExpose what) const
    {
	switch (what)
	{
	case candoClasses:
	{
                    // Create the CHEM-KW package
          std::list<std::string> nicknames = { "CK" };
          std::list<std::string> usePackages = { };
          // core::Package_sp chemkwpkg = _lisp->makePackage(ChemKwPkg,nicknames,usePackages);
          // chemkwpkg->setActsLikeKeywordPackage(true);

          core::T_sp pn;
          char* env = getenv("CANDO_LISP_SOURCE_DIR");
          if ( env != NULL ) {
            pn = core::cl__translate_logical_pathname(core::Str_O::create(std::string(env)+"/**/*.*"));
          } else {
            pn = core::cl__translate_logical_pathname(core::Str_O::create("APP-RESOURCES:CANDO;**;*.*"));
          }
          core::Cons_sp pts = core::Cons_O::createList(core::Str_O::create("cando:**;*.*"),pn);
          core::Cons_sp ptsList = core::Cons_O::createList(pts);
          core::core__pathname_translations(core::Str_O::create("CANDO"),_lisp->_true(),ptsList);
          energyFunction_initializeSmarts();
          initialize_mol2_type_rules();
	}
	break;
	case candoFunctions:
	{
          setupCandoPrimitives(_lisp);
	}
	break;
	case candoGlobals:
	{
            initialize_chimera();
	    initializeElementsAndHybridization();
//	    _lisp->defvar(_sym_candoDatabase,cdb);
	}
	break;
	case pythonClasses:
	{
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

    CandoDatabase_sp getCandoDatabase()
    {
	return _sym_candoDatabase->symbolValue().as<CandoDatabase_O>();
    }


#if 0

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




core::Symbol_sp chemkw_intern(const string& symName)
{
  if ( symName == "" ) return _Nil<core::Symbol_O>();
  const string trimmed = core::trimWhiteSpace(symName);
  core::Package_sp chemkwPkg = _lisp->findPackage(ChemKwPkg);
  core::Str_sp strimmed = core::Str_O::create(trimmed);
  core::Symbol_sp sym = chemkwPkg->intern(strimmed);
  chemkwPkg->_export2(sym);
  return sym;
//  return chemkw_intern(trimmed);
}

core::Symbol_sp chemkw_intern(core::Str_sp symName)
{
  string s = symName->get();
  if ( s == "" ) return _Nil<core::Symbol_O>();
  return chemkw_intern(s);
}

void chem_initializer()
{
//    printf("%s:%d Running chem_initializer\n", __FILE__, __LINE__ );
}

core::Initializer global_ChemInitializer(chem_initializer);

};
