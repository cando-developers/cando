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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
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
#include <cando/chem/ringFinder.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/alias.h>
#include <cando/chem/angle.h>
#include <cando/chem/atomIdMap.h>
#include <cando/chem/atomId.h>
#include <cando/chem/atomIndexer.h>
#include <cando/chem/atomReference.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/bond.h>
#include <cando/chem/chemdraw.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/loop.fwd.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/command.h>
#include <cando/chem/complexRestraints.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/coordSys.h>
#include <cando/chem/coupling.h>
#include <cando/chem/mol2.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/entity.h>
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
#include <cando/chem/stereochemistry.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <cando/chem/structureComparer.h>
#include <cando/chem/structureList.h>
#include <cando/chem/superposeEngine.h>
#include <cando/chem/topology.h>
#include <cando/chem/trajectory.h>
#include <cando/chem/twister.h>
#include <cando/chem/virtualSphere.h>
#include <cando/chem/zMatrix.h>
#include <cando/chem/chimera.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/atom.h>
#include <cando/chem/constitution.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyPointToLineRestraint.h>
#include <cando/chem/energyOutOfZPlane.h>
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyFixedNonbond.h>
#include <cando/chem/energyDihedralRestraint.h>
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
#include <cando/chem/residue.h>
#include <cando/chem/superposableConformationCollection.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/monomerPack.h>
#include <cando/chem/representedEntityNameSet.h>

SYMBOL_SHADOW_EXPORT_SC_(ChemPkg, atom);
SYMBOL_EXPORT_SC_(ChemPkg, R);
SYMBOL_EXPORT_SC_(ChemPkg, S);
SYMBOL_EXPORT_SC_(ChemPkg, UnknownConfiguration);

SYMBOL_EXPORT_SC_(ChemPkg, STARcandoDatabaseSTAR);
SYMBOL_EXPORT_SC_(ChemPkg, constitutionForNameOrPdb);
SYMBOL_EXPORT_SC_(ChemPkg, recognizesNameOrPdb);
SYMBOL_EXPORT_SC_(ChemPkg, getMonomerNameForNameOrPdb);
SYMBOL_EXPORT_SC_(ChemPkg, monomerNameForNameOrPdb);
SYMBOL_EXPORT_SC_(ChemPkg, pdbNameForNameOrPdb);
SYMBOL_EXPORT_SC_(ChemPkg, STARverboseSTAR); // slow operations print extra output if this is true

namespace chem {

//! The default CandoDatabase
SYMBOL_EXPORT_SC_(ChemPkg, candoDatabase);
SYMBOL_EXPORT_SC_(ChemPkg, AM1_BCC_ar5);
SYMBOL_EXPORT_SC_(ChemPkg, AM1_BCC_ar6);
SYMBOL_EXPORT_SC_(ChemPkg, AM1_BCC_ar7);
SYMBOL_EXPORT_SC_(ChemPkg, STARamber_charge_conversion_18_DOT_2223STAR);

extern void Initialize_Mol2_TypeRules(core::LispPtr lisp);

void ChemExposer_O::expose(core::LispPtr lisp, WhatToExpose what) const {
  switch (what) {
  case candoClasses:
      _sym_STARamber_charge_conversion_18_DOT_2223STAR->defparameter(mk_double_float(18.2223));
      _sym_STARverboseSTAR->defparameter(nil<core::T_O>());
      _sym_STARdebug_octreeSTAR->defparameter(nil<core::T_O>());
      _sym_STARdebug_octreeSTAR->defparameter(nil<core::T_O>());
      _sym_STARcurrent_matchSTAR->defparameter(nil<core::T_O>());
      _sym_STARcurrent_aromaticity_informationSTAR->defparameter(unbound<core::T_O>());
      _sym_STARcurrent_ringsSTAR->defparameter(unbound<core::T_O>());
      break;
  case candoFunctions:
      //setupCandoPrimitives(_lisp);
      break;
  case candoGlobals:
      initialize_chimera();
      initializeElementsAndHybridization();
      break;
  default:
      break;
  }
}

CandoDatabase_sp getCandoDatabase() {
  if (_sym_STARcandoDatabaseSTAR->symbolValue().notnilp()) {
    return gc::As<CandoDatabase_sp>(_sym_STARcandoDatabaseSTAR->symbolValue());
  }
  SIMPLE_ERROR("*cando-database* is not defined");
}

core::Symbol_sp chemkw_intern(const string &symName) {
  if (symName == "")
    return nil<core::Symbol_O>();
  const string trimmed = core::trimWhiteSpace(symName);
  core::Package_sp chemkwPkg = _lisp->keywordPackage();
  core::SimpleBaseString_sp strimmed = core::SimpleBaseString_O::make(trimmed);
  if (trimmed == ",") {
    printf("%s:%d chemkw_intern of %s @%p\n", __FILE__, __LINE__, trimmed.c_str(), (void *)strimmed.raw_());
  }
  core::Symbol_sp sym = chemkwPkg->intern(strimmed);
  chemkwPkg->_export2(sym);
  return sym;
}

core::Symbol_sp chemkw_intern(core::String_sp symName) {
  ASSERT(cl__stringp(symName));
  string s = symName->get_std_string();
  return (s == "") ? nil<core::Symbol_O>() : chemkw_intern(s);
}

void chem_initializer() {
}

core::Initializer global_ChemInitializer(chem_initializer);

/*! Return true of the value of chem:*verbose* is not nil and >= level
 */
DOCGROUP(cando);
CL_DEFUN bool chem__verbose(size_t level) {
  core::T_sp verbose = _sym_STARverboseSTAR->symbolValue();
  return verbose.notnilp() && (!verbose.fixnump() || verbose.unsafe_fixnum() >= level);
}

}; // namespace chem
