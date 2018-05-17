#define USE_KINEMATICS
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

#ifdef GC_INTERFACE_INCLUDES
#include <llvm/ExecutionEngine/GenericValue.h>
#include <llvm/ExecutionEngine/SectionMemoryManager.h>
#include <llvm/IR/LLVMContext.h>
#include <llvm/CodeGen/LinkAllCodegenComponents.h>
#include <llvm/Bitcode/BitcodeReader.h>
#include <llvm/Bitcode/BitcodeWriter.h>
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
#include <llvm/IR/LegacyPassManager.h>
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
#include <cando/chem/scoringFunction.h>
#include <cando/chem/energyFunction.h>
#include <cando/chem/energyImproperRestraint.fwd.h>
#include <cando/chem/energyImproperRestraint.h>
#include <cando/chem/energyNonbond.fwd.h>
#include <cando/chem/energyNonbond.h>
#include <cando/chem/energyStretch.fwd.h>
#include <cando/chem/energyStretch.h>
#include <cando/chem/rigidBodyEnergyFunction.h>
#include <cando/chem/energyRigidBodyStaple.h>
#include <cando/chem/energyRigidBodyNonbond.h>
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
#include <cando/chem/octree.h>

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

#ifdef USE_KINEMATICS
#include <cando/kinematics/aggregateNode.h>
#include <cando/kinematics/atom.fwd.h>
#include <cando/kinematics/atom.h>
#include <cando/kinematics/atomHandle.h>
#include <cando/kinematics/atomTemplate.fwd.h>
#include <cando/kinematics/atomTemplate.h>
#include <cando/kinematics/atomTree.fwd.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/atomo.h>
#include <cando/kinematics/bondId.fwd.h>
#include <cando/kinematics/bondId.h>
#include <cando/kinematics/bondedAtom.h>
#include <cando/kinematics/chainNode.fwd.h>
#include <cando/kinematics/chainNode.h>
#include <cando/kinematics/chiDihedrals.fwd.h>
#include <cando/kinematics/chiDihedrals.h>
#include <cando/kinematics/conformation.fwd.h>
#include <cando/kinematics/conformation.h>
#include <cando/kinematics/delayedBondedAtom.h>
#include <cando/kinematics/dofId.h>
#include <cando/kinematics/dofType.h>
#include <cando/kinematics/foldTree.h>
#include <cando/kinematics/jump.h>
#include <cando/kinematics/jumpAtom.h>
#include <cando/kinematics/kinFoundation.h>
#include <cando/kinematics/kinematicsPackage.h>
#include <cando/kinematics/monomerBaseNode.fwd.h>
#include <cando/kinematics/monomerBaseNode.h>
#include <cando/kinematics/monomerId.h>
#include <cando/kinematics/monomerIdMap.h>
#include <cando/kinematics/monomerNode.fwd.h>
#include <cando/kinematics/monomerNode.h>
#include <cando/kinematics/originJumpAtom.h>
#include <cando/kinematics/otherPackageClasses.h>
#include <cando/kinematics/pool.h>
#include <cando/kinematics/ringClosingMonomerNode.h>
#include <cando/kinematics/rootAtomInfo.h>
#include <cando/kinematics/rootBondedAtom.h>
#include <cando/kinematics/rotamer.h>
#include <cando/kinematics/stub.fwd.h>
#include <cando/kinematics/stub.h>
#endif

#endif

#ifdef GC_INTERFACE_GC_MANAGED_TYPES
GC_MANAGED_TYPE(gctools::GCArray_moveable<Vector3>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<double>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<float>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<gctools::smart_ptr<core::T_O>>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<int>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<long>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<short>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<signed char>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<unsigned char>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<unsigned int>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<unsigned long>);
GC_MANAGED_TYPE(gctools::GCArray_moveable<unsigned short>);
GC_MANAGED_TYPE(gctools::GCBitUnitArray_moveable<1,unsigned int,int>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::AtomInfo>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::AtomPdbRec>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::BeyondThresholdFixedNonbondRestraint>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyAnchorRestraint>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyAngle>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyAtom>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyChiralRestraint>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyDihedral>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyImproperRestraint>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyNonbond>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyRigidBodyStaple>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::EnergyStretch>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::FixedNonbondRestraint>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::KeyEntry>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::ResidueOut>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<chem::RigidBodyAtomInfo>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::AuxArgument>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::CacheRecord>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::DynamicBinding>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::ExceptionEntry>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::KeywordArgument>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::OptionalArgument>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::RequiredArgument>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::SymbolClassPair>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::SymbolStorage>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<core::T_O *>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<double>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<float>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::SmallOrderedSet<gctools::smart_ptr<chem::Atom_O>>>);

#ifdef USE_KINEMATICS
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<kinematics::AtomTemplate_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<kinematics::MonomerNode_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<kinematics::ChainNode_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<kinematics::Rotamer_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<kinematics::ChiDihedral_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<kinematics::BackboneDependentRotamerSet_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>, gctools::smart_ptr<kinematics::MonomerNode_O> > > );
#endif


GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::StereoisomerAtoms_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<adapt::QDomNode_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::AGEdge_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Atom_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::BondList_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Bond_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ChemInfo_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::CDBond_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::CDFragment_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Isoname_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::AntechamberRoot_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ComplexRestraint_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ConformationCollectionEntry_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ConformationExplorerEntryStage_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ConformationExplorerEntry_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ConstitutionAtom_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ConstitutionBond_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Coupling_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::FFNonbond_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::FFVdw_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Mate_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Matter_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Molecule_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::MonomerCoordinates_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Monomer_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::OneTypeRule_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::PathMessage_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ReportBase_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Residue_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Rotamer_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::StereoConfiguration_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::StereoisomerAtom_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Stereoisomer_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Structure_Old_ListEntry_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Topology_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::TrajectoryFrame_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::Twister_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ZMatrixEntry_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::ZMatrixInternal_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<chem::OctNode_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<clbind::ClassRep_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::Cons_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::List_V>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::Package_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::SequenceStepper_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::SimpleBitVector_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::SingleDispatchMethod_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::SourceFileInfo_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::Symbol_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<gctools::smart_ptr<core::T_O>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<int>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<chem::Atom_O>,gctools::smart_ptr<chem::CDNode_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<chem::Atom_O>,int>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<chem::Atom_O>,unsigned int>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,chem::EstimateStretch>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::AtomIndexer_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::CDFragment_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::Coupling_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::EntityNameSetBase_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::Entity_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::FrameRecognizer_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::Molecule_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::MonomerCoordinates_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::PdbMonomerConnectivity_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::Plug_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::SpecificContext_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::StereoisomerAtoms_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::Stereoisomer_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<chem::Topology_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<core::Array_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<core::DoubleFloat_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<core::Symbol_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::Symbol_O>,gctools::smart_ptr<core::T_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<gctools::smart_ptr<core::T_O>,gctools::smart_ptr<core::T_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<int,gctools::smart_ptr<chem::CDNode_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<int,gctools::smart_ptr<chem::Residue_O>>>);
GC_MANAGED_TYPE(gctools::GCVector_moveable<std::pair<unsigned int,gctools::smart_ptr<chem::Residue_O>>>);
#endif

#ifdef GC_INTERFACE_FORWARD
class Vector3;
 namespace geom {
    class OVector2_O;
    class SimpleMDArrayCoordinate_O;
    class BoundingBox_O;
    class OMatrix_O;
    class GeomExposer_O;
    class OVector3_O;
    class SimpleVectorCoordinate_O;
    class Color_O;
    class MDArrayCoordinate_O;
 };
 namespace units {
    class Unit_O;
    class NamedUnit_O;
    class UnitsExposer_O;
    class Quantity_O;
    class Dimension_O;
 };
 namespace adapt {
    class StringSet_O;
    class AdaptExposer_O;
    class SymbolSet_O;
    class StringList_O;
    class ObjectSet_O;
    class IndexedObjectBag_O;
    class SymbolList_O;
    class IterateCons_O;
    class QDomNode_O;
 };
 namespace chem {
    class MessageReport_O;
    class EntityNameSetBase_O;
    class AtomTable_O;
    class EnergyRigidBodyStaple_O;
    class Command_O;
    class IterateResidues_O;
    class TagSet_O;
    class DirectionalCoupling_O;
    class EnergyChiralRestraint_O;
    class EnergyImproperRestraint_O;
    class Root_O;
    class ResidueOut;
    class Monomer_O;
    class EnergyImproperRestraint;
    class MonomerContext_O;
    class FFStretch_O;
    class FFNonbondCrossTermTable_O;
    class CDFragment_O;
    class BondList_O;
    class ResidueList_O;
    class Structure_Old_ListEntry_O;
    class OutPlug_O;
    class OligomerPart_Monomer_O;
    class ReadAmberParameters_O;
    class AtomIndexer_O;
    class Residue_O;
    class RestrainedExoCyclicAtom_O;
    class EnergyChiralRestraint;
    class Structure_Old_List_O;
    class ChemInfoMatch_O;
    class SpanningLoop_O;
    class CDBond_O;
    class OligomerPart_Link_O;
    class ChemExposer_O;
    class CalculatePositionUsingInternals_O;
    class EnergyAtom;
    class RestrainedPiBond_O;
    class RingClosingMate_O;
    class ForceField_O;
    class ReportBase_O;
    class FFAngle_O;
    class FrameRecognizer_O;
    class PathMessage_O;
    class StereoisomerAtom_O;
    class MonomerCoordinates_O;
    class SearchStatistics_O;
    class RestraintAngle_O;
    class ZMatrixBondInternal_O;
    class FFVdw_O;
    class FFStretchDb_O;
    class EnergyAnchorRestraint_O;
    class AtomId_O;
    class TypeAssignmentRules_O;
    class SmartsRoot_O;
    class EnergyNonbond;
    class WildElementDict_O;
    class PlugWithMates_O;
    class EnergyAnchorRestraint;
    class EnergyAngle_O;
    class ResidueTest_O;
    class SpecificContext_O;
    class EnergyRigidBodyStaple;
    class AtomOrBondMatchNode_O;
    class AntechamberFocusAtomMatch_O;
    class CDText_O;
    class RestraintAnchor_O;
    class PdbMonomerConnectivity_O;
    class IterateBonds_O;
    class Plug_O;
    class IterateAtoms_O;
    class SparseLargeSquareMatrix_O;
    class RepresentedEntityNameSet_O;
    class SuperposeSelectedAtoms_O;
    class MinimizerLog_O;
    class FFTypesDb_O;
    class Branch_O;
    class ChemInfo_O;
    class ZMatrixEntry_O;
    class RestraintFixedNonbond_O;
    class TopologyAtomInfo_O;
    class StereoInformation_O;
    class MacroModelFile_O;
    class RestraintChiral_O;
    class PdbWriter_O;
    class ChemDraw_O;
    class ConformationExplorerMatch_O;
    class AGEdge_O;
    class ImproperTorsion_O;
    class ConstitutionAtoms_O;
    class Aggregate_O;
    class Atom_O;
    class CoordinateSystem_O;
    class JumpPlug_O;
    class ConformationCollectionEntry_O;
    class Restraint_O;
    class ChemInfoNode_O;
    class FFNonbond_O;
    class FFPtorDb_O;
    class FFPtor_O;
    class AtomInfo;
    class Molecule_O;
    class RingTest_O;
    class EnergyDihedral_O;
    class RootMatchNode_O;
    class Minimizer_O;
    class Oligomer_O;
    class ConstitutionBond_O;
    class Entity_O;
    class RigidBodyAtomInfo;
    class StereoisomerAtoms_O;
    class EntityNameSetWithCap_O;
    class EnergyRigidBodyNonbond_O;
    class ConformationExplorerEntry_O;
    class AtomReferenceBase_O;
    class CDNode_O;
    class FixedNonbondRestraint;
    class Chain_O;
    class TrajectoryFrame_O;
    class EnergyComponent_O;
    class PdbMonomerDatabase_O;
    class Coupling_O;
    class KeyEntry;
    class PdbReader_O;
    class EnergyStretch_O;
    class ZMatrixInternal_O;
    class CipPrioritizer_O;
    class Matter_O;
    class Mate_O;
    class EstimateStretch;
    class Bond_O;
    class RestraintDihedral_O;
    class AGVertex_O;
    class FullLargeSquareMatrix_O;
    class FFItorDb_O;
    class BondTest_O;
    class Alias_O;
    class CandoDatabase_O;
    class AtomReference_O;
    class ProperTorsion_O;
    class EnergyFunction_O;
    class ComplexRestraint_O;
    class CalculatePositionRelativeToOrigin_O;
    class IterateRestraints_O;
    class AntechamberRoot_O;
    class ZMatrix_O;
    class Stereoisomer_O;
    class AntechamberBondTest_O;
    class IterateMatter_O;
    class BondMatchNode_O;
    class EnergyRigidBodyComponent_O;
    class StatusTracker_O;
    class AtomPdbRec;
    class OneTypeRule_O;
    class SuperposableConformationCollection_O;
    class ZMatrixAngleInternal_O;
    class EnergyNonbond_O;
    class RingClosingPlug_O;
    class OligomerPart_Base_O;
    class FFAngleDb_O;
    class CalculatePositionAlongBond_O;
    class EnergyAngle;
    class EntityNameSet_O;
    class ConformationCollection_O;
    class InfoDb_O;
    class ConformationExplorer_O;
    class RigidBodyEnergyFunction_O;
    class ConstitutionVirtualAtom_O;
    class ConstitutionAtom_O;
    class NumericalFunction_O;
    class AbstractLargeSquareMatrix_O;
    class TwisterDriver_O;
    class ZMatrixDihedralInternal_O;
    class FFNonbondDb_O;
    class VirtualAtom_O;
    class EnergyDihedral;
    class VirtualSphere_O;
    class Rotamer_O;
    class FFItor_O;
    class Topology_O;
    class StereoConfiguration_O;
    class FFParameterBaseDb_O;
    class BondListMatchNode_O;
    class Trajectory_O;
    class MapOfMonomerNamesToAtomIndexers_O;
    class ScoringFunction_O;
    class RingCoupling_O;
    class AfterMatchBondTest_O;
    class RingFinder_O;
    class Twister_O;
    class FFBaseDb_O;
    class InPlug_O;
    class MonomerPack_O;
    class SuperposeEngine_O;
    class CalculatePosition_O;
    class ForceMatchReport_O;
    class StereoisomerVirtualAtom_O;
    class Constitution_O;
    class FFParameter_O;
    class EnergyStretch;
    class RestraintDistance_O;
    class SpecificContextSet_O;
    class ConformationExplorerEntryStage_O;
    class Angle_O;
    class EnergyFixedNonbondRestraint_O;
    class BeyondThresholdFixedNonbondRestraint;
    class StructureComparer_O;
    class AtomTest_O;
    class Logical_O;
    class StepReport_O;
    class AtomIdToAtomMap_O;
    class FFVdwDb_O;
    class OctNode_O;
    class Octree_O;
 };
#endif // GC_INTERFACE_FORWARD
