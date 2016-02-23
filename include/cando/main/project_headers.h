/*
    File: project_headers.h
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

#include <llvm/ExecutionEngine/GenericValue.h>
#include <llvm/ExecutionEngine/SectionMemoryManager.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/CodeGen/LinkAllCodegenComponents.h>
#include <llvm/Bitcode/ReaderWriter.h>
#include <llvm/Support/raw_ostream.h>
#include <llvm/Support/MemoryBuffer.h>
#include <llvm/Support/FileSystem.h>
#include <llvm/ADT/Triple.h>
#include <llvm/IR/Module.h>
#include <llvm/IR/DerivedTypes.h>
#include <llvm/IR/Constants.h>
#include <llvm/IR/GlobalVariable.h>
#include <llvm/IR/Function.h>
#include <llvm/IR/CallingConv.h>
#include <llvm/IR/BasicBlock.h>
#include <llvm/IR/Instructions.h>
#include <llvm/Transforms/Instrumentation.h>
#include <llvm/Transforms/IPO.h>
#include <llvm/IR/InlineAsm.h>
#include <llvm/Support/FormattedStream.h>
#include <llvm/Support/TargetRegistry.h>
#include <llvm/Support/MathExtras.h>
#include <llvm/Pass.h>
#include <llvm/PassManager.h>
#include <llvm/ADT/SmallVector.h>
#include <llvm/IR/Verifier.h>
#include "llvm/IR/AssemblyAnnotationWriter.h" // will be llvm/IR



#include <cando/adapt/adaptPackage.fwd.h>
#include <cando/adapt/adaptPackage.h>
#include <cando/adapt/adapters.h>
#include <cando/adapt/indexedObjectBag.h>
#include <cando/adapt/iterateCons.fwd.h>
#include <cando/adapt/iterateCons.h>
#include <cando/adapt/mySaxInterface.h>
#include <cando/adapt/objectSet.fwd.h>
#include <cando/adapt/objectSet.h>
#include <cando/adapt/quickDom.fwd.h>
#include <cando/adapt/quickDom.h>
#include <cando/adapt/stringList.fwd.h>
#include <cando/adapt/stringList.h>
#include <cando/adapt/stringSet.fwd.h>
#include <cando/adapt/stringSet.h>
#include <cando/adapt/symbolList.fwd.h>
#include <cando/adapt/symbolList.h>
#include <cando/adapt/symbolMap.h>
#include <cando/adapt/symbolSet.fwd.h>
#include <cando/adapt/symbolSet.h>






#include <cando/chem/aggregate.fwd.h>
#include <cando/chem/aggregate.h>
#include <cando/chem/alchemist.fwd.h>
#include <cando/chem/alias.fwd.h>
#include <cando/chem/alias.h>
#include <cando/chem/angle.fwd.h>
#include <cando/chem/angle.h>
#include <cando/chem/atom.fwd.h>
#include <cando/chem/atom.h>
#include <cando/chem/atomId.fwd.h>
#include <cando/chem/atomId.h>
#include <cando/chem/atomIdMap.fwd.h>
#include <cando/chem/atomIdMap.h>
#include <cando/chem/atomIndexer.fwd.h>
#include <cando/chem/atomIndexer.h>
#include <cando/chem/atomReference.fwd.h>
#include <cando/chem/atomReference.h>
#include <cando/chem/bond.fwd.h>
#include <cando/chem/bond.h>
#include <cando/chem/calculatePosition.fwd.h>
#include <cando/chem/calculatePosition.h>
#include <cando/chem/candoDatabase.fwd.h>
#include <cando/chem/candoDatabase.h>
#include <cando/chem/candoDatabaseDependent.h>
#include <cando/chem/candoScript.fwd.h>
#include <cando/chem/candoScript.h>
#include <cando/chem/chemInfo.fwd.h>
#include <cando/chem/chemInfo.h>
#include <cando/chem/chemPackage.fwd.h>
#include <cando/chem/chemPackage.h>
#include <cando/chem/chemdraw.fwd.h>
#include <cando/chem/chemdraw.h>
#include <cando/chem/chimera.h>
#include <cando/chem/cipPrioritizer.fwd.h>
#include <cando/chem/cipPrioritizer.h>
#include <cando/chem/command.fwd.h>
#include <cando/chem/command.h>
#include <cando/chem/complexRestraints.fwd.h>
#include <cando/chem/complexRestraints.h>
#include <cando/chem/conformationCollection.fwd.h>
#include <cando/chem/conformationCollection.h>
#include <cando/chem/conformationExplorer.fwd.h>
#include <cando/chem/conformationExplorer.h>
#include <cando/chem/constitution.fwd.h>
#include <cando/chem/constitution.h>
#include <cando/chem/constitutionAtoms.fwd.h>
#include <cando/chem/constitutionAtoms.h>
#include <cando/chem/coordSys.fwd.h>
#include <cando/chem/coordSys.h>
#include <cando/chem/coupling.fwd.h>
#include <cando/chem/coupling.h>
#include <cando/chem/elements.fwd.h>
#include <cando/chem/elements.h>
#include <cando/chem/energyAnchorRestraint.fwd.h>
#include <cando/chem/energyAnchorRestraint.h>
#include <cando/chem/energyAngle.fwd.h>
#include <cando/chem/energyAngle.h>
#include <cando/chem/energyAtomTable.fwd.h>
#include <cando/chem/energyAtomTable.h>
#include <cando/chem/energyChiralRestraint.fwd.h>
#include <cando/chem/energyChiralRestraint.h>
#include <cando/chem/energyComponent.fwd.h>
#include <cando/chem/energyComponent.h>
#include <cando/chem/energyDihedral.fwd.h>
#include <cando/chem/energyDihedral.h>
#include <cando/chem/energyFixedNonbond.fwd.h>
#include <cando/chem/energyFixedNonbond.h>
#include <cando/chem/energyFunction.fwd.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/energyImproperRestraint.fwd.h>
#include <cando/chem/energyImproperRestraint.h>
#include <cando/chem/energyNonbond.fwd.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyStretch.fwd.h>
#include <cando/chem/energyStretch.h>
#include <cando/chem/entity.fwd.h>
#include <cando/chem/entity.h>
#include <cando/chem/entityNameSet.fwd.h>
#include <cando/chem/entityNameSet.h>
#include <cando/chem/ffAngleDb.fwd.h>
#include <cando/chem/ffAngleDb.h>
#include <cando/chem/ffBaseDb.fwd.h>
#include <cando/chem/ffBaseDb.h>
#include <cando/chem/ffItorDb.fwd.h>
#include <cando/chem/ffItorDb.h>
#include <cando/chem/ffNonbondDb.fwd.h>
#include <cando/chem/ffNonbondDb.h>
#include <cando/chem/ffPtorDb.fwd.h>
#include <cando/chem/ffPtorDb.h>
#include <cando/chem/ffStretchDb.fwd.h>
#include <cando/chem/ffStretchDb.h>
#include <cando/chem/ffTypesDb.fwd.h>
#include <cando/chem/ffTypesDb.h>
#include <cando/chem/ffVdwDb.fwd.h>
#include <cando/chem/ffVdwDb.h>
#include <cando/chem/forceField.fwd.h>
#include <cando/chem/forceField.h>
#include <cando/chem/frameRecognizer.fwd.h>
#include <cando/chem/frameRecognizer.h>
#include <cando/chem/gamessInterface.fwd.h>
#include <cando/chem/gaussianInterface.fwd.h>
#include <cando/chem/hold.h>
#include <cando/chem/improperTorsion.fwd.h>
#include <cando/chem/improperTorsion.h>
#include <cando/chem/indirectAtomCoordinateReference.h>
#include <cando/chem/iterateMatter.fwd.h>
#include <cando/chem/iterateMatter.h>
#include <cando/chem/iterateRestraints.fwd.h>
#include <cando/chem/iterateRestraints.h>
#include <cando/chem/largeSquareMatrix.fwd.h>
#include <cando/chem/largeSquareMatrix.h>
#include <cando/chem/linearAlgebra.h>
#include <cando/chem/loop.fwd.h>
#include <cando/chem/loop.h>
#include <cando/chem/macroModel.h>
#include <cando/chem/matter.fwd.h>
#include <cando/chem/matter.h>
#include <cando/chem/mbbCoreTools.h>
#include <cando/chem/minimizer.fwd.h>
#include <cando/chem/minimizer.h>
#include <cando/chem/minimizerLog.fwd.h>
#include <cando/chem/minimizerLog.h>
#include <cando/chem/moe.h>
#include <cando/chem/mol2.h>
#include <cando/chem/molecule.fwd.h>
#include <cando/chem/molecule.h>
#include <cando/chem/monomer.fwd.h>
#include <cando/chem/monomer.h>
#include <cando/chem/monomerContext.fwd.h>
#include <cando/chem/monomerContext.h>
#include <cando/chem/monomerCoordinates.fwd.h>
#include <cando/chem/monomerCoordinates.h>
#include <cando/chem/monomerGroup.fwd.h>
#include <cando/chem/monomerPack.fwd.h>
#include <cando/chem/monomerPack.h>
#include <cando/chem/nVector.fwd.h>
#include <cando/chem/nVector.h>
#include <cando/chem/numericalFunction.fwd.h>
#include <cando/chem/numericalFunction.h>
#include <cando/chem/oligomer.fwd.h>
#include <cando/chem/oligomer.h>
#include <cando/chem/pdb.fwd.h>
#include <cando/chem/pdb.h>
#include <cando/chem/pdbMonomer.fwd.h>
#include <cando/chem/pdbMonomer.h>
#include <cando/chem/plug.fwd.h>
#include <cando/chem/plug.h>
#include <cando/chem/properTorsion.fwd.h>
#include <cando/chem/properTorsion.h>
#include <cando/chem/randomGenerators.h>
#include <cando/chem/readAmberParameters.fwd.h>
#include <cando/chem/readAmberParameters.h>
#include <cando/chem/representedEntityNameSet.h>
#include <cando/chem/residue.fwd.h>
#include <cando/chem/residue.h>
#include <cando/chem/restraint.fwd.h>
#include <cando/chem/restraint.h>
#include <cando/chem/ringFinder.fwd.h>
#include <cando/chem/ringFinder.h>
//#include <cando/chem/rotamer.fwd.h>
#include <cando/chem/rotamer.h>
#include <cando/chem/searchStatistics.fwd.h>
#include <cando/chem/searchStatistics.h>
#include <cando/chem/spanningLoop.fwd.h>
#include <cando/chem/spanningLoop.h>
#include <cando/chem/sparseMatrix.h>
#include <cando/chem/specificContext.fwd.h>
#include <cando/chem/specificContext.h>
#include <cando/chem/spline.h>
#include <cando/chem/statusTracker.fwd.h>
#include <cando/chem/statusTracker.h>
#include <cando/chem/stereochemistry.fwd.h>
#include <cando/chem/stereochemistry.h>
#include <cando/chem/stereoisomerAtoms.fwd.h>
#include <cando/chem/stereoisomerAtoms.h>
#include <cando/chem/structureComparer.fwd.h>
#include <cando/chem/structureComparer.h>
#include <cando/chem/structureList.fwd.h>
#include <cando/chem/structureList.h>
#include <cando/chem/superposableConformationCollection.fwd.h>
#include <cando/chem/superposableConformationCollection.h>
#include <cando/chem/superposeEngine.fwd.h>
#include <cando/chem/superposeEngine.h>
#include <clasp/core/symbolTable.h>
#include <cando/chem/topology.fwd.h>
#include <cando/chem/topology.h>
//#include <cando/chem/torsionDriver.h>
#include <cando/chem/trajectory.fwd.h>
#include <cando/chem/trajectory.h>
#include <cando/chem/twister.fwd.h>
#include <cando/chem/twister.h>
#include <cando/chem/typeAssignmentRules.fwd.h>
#include <cando/chem/typeAssignmentRules.h>
#include <cando/chem/units.h>
#include <cando/chem/units_scrape_flag.h>
#include <cando/chem/virtualAtom.fwd.h>
#include <cando/chem/virtualAtom.h>
#include <cando/chem/virtualSphere.fwd.h>
#include <cando/chem/virtualSphere.h>
#include <cando/chem/zMatrix.h>

#include <cando/geom/boundingBox.fwd.h>
#include <cando/geom/boundingBox.h>
#include <cando/geom/color.fwd.h>
#include <cando/geom/color.h>
#include <cando/geom/coordinateArray.fwd.h>
#include <cando/geom/coordinateArray.h>
#include <cando/geom/geomPackage.fwd.h>
#include <cando/geom/geomPackage.h>
#include <cando/geom/matrix.h>
#include <cando/geom/omatrix.fwd.h>
#include <cando/geom/omatrix.h>
#include <cando/geom/otherPackageClasses.h>
#include <cando/geom/ovector2.fwd.h>
#include <cando/geom/ovector2.h>
#include <cando/geom/ovector3.fwd.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/symbolTable.h>
#include <cando/geom/vector2.h>
#include <cando/geom/vector3.h>
#include <cando/geom/vector4.h>

#include <cando/units/dimension.fwd.h>
#include <cando/units/dimension.h>
#include <cando/units/kinFoundation.h>
#include <cando/units/namedUnit.h>
#include <cando/units/otherPackageClasses.h>
#include <cando/units/quantity.fwd.h>
#include <cando/units/quantity.h>
#include <clasp/core/symbolTable.h>
#include <cando/units/unit.h>
#include <cando/units/unitsFoundation.h>
#include <cando/units/unitsPackage.h>
