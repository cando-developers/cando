// start
// define cpp macros: SET_SYMBOL, CREATE_CLASS, SET_CLASS, DEFINE_BASE_CLASSES, DEFINE_CLASS_NAMES, EXPOSE_TO_CANDO 
// define cpp macro: ALL_STAGES to get the effect of defining all of the macros above
// define cpp macro: EXPOSE_PYTHON to expose python
// Associating namespace(chemkw) with package(ChemKwPkg)
// Associating namespace(chem) with package(ChemPkg)
 // class  Aggregate_O : public  Matter_O
 // class Alias_O : public core::T_O
 // class Angle_O : public core::T_O
 // class Atom_O : public Matter_O
 // class AtomId_O : public core::T_O
 // class AtomIdToAtomMap_O : public core::T_O
 // class AtomIndexer_O : public core::T_O
 // class MapOfMonomerNamesToAtomIndexers_O : public core::T_O
 // class AtomReferenceBase_O : public core::T_O
 // class AtomReference_O : public AtomReferenceBase_O
 // class Bond_O : public core::T_O
 // class BondList_O : public core::T_O
 // class CalculatePosition_O : public core::T_O
 // class CalculatePositionUsingInternals_O : public CalculatePosition_O
 // class CalculatePositionAlongBond_O : public CalculatePosition_O
 // class CalculatePositionRelativeToOrigin_O : public CalculatePosition_O
 // class CandoDatabase_O : public core::T_O
 // class OligomerPart_Base_O : public core::T_O
 // class OligomerPart_Monomer_O : public OligomerPart_Base_O
 // class OligomerPart_Link_O : public OligomerPart_Base_O
 // class ChemInfoMatch_O : public core::T_O
 // class WildElementDict_O : public core::T_O
 // class ChemInfoNode_O : public core::T_O
 // class ResidueList_O : public core::T_O
 // class RootMatchNode_O : public ChemInfoNode_O
 // class BondMatchNode_O : public ChemInfoNode_O
 // class AtomOrBondMatchNode_O : public BondMatchNode_O
 // class BondListMatchNode_O : public ChemInfoNode_O
 // class Logical_O : public AtomOrBondMatchNode_O
 // class TagSet_O : public AtomOrBondMatchNode_O
 // class RingTest_O : public AtomOrBondMatchNode_O
 // class ResidueTest_O : public AtomOrBondMatchNode_O
 // class BondTest_O : public BondMatchNode_O
 // class AntechamberBondTest_O : public BondMatchNode_O
 // class AtomTest_O : public AtomOrBondMatchNode_O
 // class AntechamberFocusAtomMatch_O : public AtomOrBondMatchNode_O
 // class Chain_O : public BondListMatchNode_O
 // class Branch_O : public BondListMatchNode_O
 // class AfterMatchBondTest_O : public RootMatchNode_O
 // class Root_O : public AtomOrBondMatchNode_O
 // class SmartsRoot_O : public Root_O
 // class AntechamberRoot_O : public Root_O
 // class ChemInfo_O : public core::T_O
// Associating namespace(chemkw) with package(ChemKwPkg)
// Associating namespace(chem) with package(ChemPkg)
 // class CDNode_O : public core::T_O
 // class CDBond_O : public core::T_O
 // class CDFragment_O : public core::T_O
 // class CDText_O : public core::T_O
 // class ChemDraw_O : public core::T_O
 // class CipPrioritizer_O : public core::T_O
 // class Command_O : public core::T_O
 // class ComplexRestraint_O : public core::T_O
 // class RestrainedExoCyclicAtom_O : public ComplexRestraint_O
 // class RestrainedPiBond_O : public ComplexRestraint_O
 // class ConformationCollectionEntry_O : public core::T_O
 // class ConformationCollection_O : public core::T_O
 // class ConformationExplorerEntryStage_O : public core::T_O
 // class ConformationExplorerEntry_O : public core::T_O
 // class ConformationExplorer_O : public core::T_O
 // class ConformationExplorerMatch_O : public core::T_O
 // class Constitution_O : public Entity_O
 // class ConstitutionBond_O : public core::T_O
 // class ConstitutionAtom_O : public core::T_O
 // class ConstitutionVirtualAtom_O : public ConstitutionAtom_O
 // class ConstitutionAtoms_O : public core::T_O
 // class CoordinateSystem_O : public core::T_O
 // class Coupling_O : public core::T_O //public core::Model_O
 // class DirectionalCoupling_O : public Coupling_O
 // class RingCoupling_O : public Coupling_O
 // class EnergyAnchorRestraint_O : public EnergyComponent_O
 // class EnergyAngle_O : public EnergyComponent_O
 // class AtomTable_O : public core::T_O
 // class EnergyChiralRestraint_O : public EnergyComponent_O
 // class EnergyComponent_O : public core::T_O
 // class EnergyDihedral_O : public EnergyComponent_O
 // class EnergyFixedNonbondRestraint_O : public EnergyComponent_O
 // class EnergyFunction_O : public core::T_O
 // class EnergyImproperRestraint_O : public EnergyComponent_O
 // class EnergyNonbond_O : public EnergyComponent_O
 // class EnergyStretch_O : public EnergyComponent_O
 // class Entity_O : public core::T_O
 // class EntityNameSetBase_O : public Entity_O
 // class EntityNameSet_O : public EntityNameSetBase_O
 // class EntityNameSetWithCap_O : public EntityNameSetBase_O
 // class FFAngle_O : public FFParameter_O
 // class FFAngleDb_O : public FFBaseDb_O
 // class FFParameter_O : public core::T_O
 // class FFBaseDb_O : public core::T_O
 // class FFItor_O : public FFParameter_O
 // class FFItorDb_O : public FFBaseDb_O
 // class FFNonbondCrossTermTable_O : public core::T_O
 // class FFNonbond_O : public FFParameter_O
 // class FFNonbondDb_O : public FFBaseDb_O
 // class FFPtor_O : public FFParameter_O
 // class FFPtorDb_O : public FFBaseDb_O
 // class FFStretch_O : public FFParameter_O
 // class FFStretchDb_O : public FFBaseDb_O
 // class FFTypesDb_O : public FFBaseDb_O
 // class FFVdw_O : public FFParameter_O
 // class FFVdwDb_O : public FFBaseDb_O
 // class InfoDb_O : public core::T_O
 // class ForceField_O : public core::T_O
 // class FrameRecognizer_O : public core::T_O
 // class ImproperTorsion_O : public core::T_O
 // class IterateMatter_O : public core::Iterator_O
 // class IterateAtoms_O : public IterateMatter_O
 // class IterateResidues_O : public IterateMatter_O
 // class IterateBonds_O : public IterateMatter_O
 // class IterateRestraints_O : public core::Iterator_O
 // class AbstractLargeSquareMatrix_O : public core::T_O
 // class FullLargeSquareMatrix_O : public  AbstractLargeSquareMatrix_O
 // class SparseLargeSquareMatrix_O : public AbstractLargeSquareMatrix_O
 // class MacroModelFile_O : public core::T_O
 // class Matter_O : public core::T_O
 // class Minimizer_O : public core::T_O
 // class ReportBase_O : public core::T_O
 // class MessageReport_O : public ReportBase_O
 // class ForceMatchReport_O : public ReportBase_O
 // class StepReport_O : public ReportBase_O
 // class MinimizerLog_O : public core::T_O
 // class Molecule_O : public Matter_O
 // class Monomer_O : public Entity_O
 // class OneMonomer_O : public core::T_O
 // class MultiMonomer_O : public Monomer_O
 // class MonomerContext_O : public core::T_O
 // class MonomerCoordinates_O : public core::T_O
 // class MonomerPack_O : public EntityNameSet_O
 // class NVector_O : public core::T_O
 // class NumericalFunction_O : public core::T_O
 // class Oligomer_O : public core::T_O
 // LISP_OTHER_PACKAGE(CorePkg,geom::HighlightedObject_O,core::T_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::Iterator_O,core::T_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::ConditionBase_O,core::T_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::ConditionWithStackTrace_O,core::ConditionBase_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::CandoException_O,core::ConditionWithStackTrace_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::LispError_O,core::CandoException_O)
 // LISP_OTHER_PACKAGE(CorePkg,geom::GraphicsObject_O,core::T_O)
 // LISP_OTHER_PACKAGE(CorePkg,geom::Render_O,geom::GraphicsObject_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::SequenceBaseVirtual_O,core::T_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::SequenceBase_O,core::SequenceBaseVirtual_O)
 // LISP_OTHER_PACKAGE(CorePkg,core::Sequence_O,core::SequenceBase_O)
 // class PdbReader_O : public core::T_O
 // class PdbWriter_O : public core::T_O
 // class PdbMonomerConnectivity_O : public core::T_O
 // class PdbMonomerDatabase_O : public core::T_O
 // class Mate_O : public EntityNameSetWithCap_O
 // class RingClosingMate_O : public EntityNameSetWithCap_O
 // class Plug_O : public core::T_O
 // class PlugWithMates_O : public Plug_O
 // class OutPlug_O : public PlugWithMates_O
 // class InPlug_O : public PlugWithMates_O
 // class JumpPlug_O : public Plug_O
 // class RingClosingPlug_O : public OutPlug_O
 // class ProperTorsion_O : public core::T_O
 // class ReadAmberParameters_O : public core::T_O
 // class RepresentedEntityNameSet_O : public EntityNameSet_O
 // class Residue_O : public Matter_O
 // class Restraint_O : public core::T_O
 // class RestraintAnchor_O : public Restraint_O
 // class RestraintChiral_O : public Restraint_O
 // class RestraintDistance_O : public Restraint_O
 // class RestraintAngle_O : public Restraint_O
 // class RestraintDihedral_O : public Restraint_O
 // class RestraintFixedNonbond_O : public Restraint_O
 // class RestraintList_O : public core::Vector_O
 // class PathMessage_O : public core::T_O
 // class AGVertex_O : public core::T_O
 // class AGEdge_O : public core::T_O
 // class RingFinder_O : public core::T_O
 // class Rotamer_O : public core::T_O
 // class SearchStatistics_O : public core::T_O
 // class SpanningLoop_O : public core::T_O
 // class SpecificContext_O : public core::T_O
 // class SpecificContextSet_O : public core::T_O
 // class StatusTracker_O : public core::T_O
 // class StereoConfiguration_O : public core::T_O
 // class Stereoisomer_O : public Entity_O
 // class StereoInformation_O : public core::T_O
 // class StereoisomerAtom_O : public core::T_O
 // class StereoisomerVirtualAtom_O : public StereoisomerAtom_O
 // class StereoisomerAtoms_O : public core::T_O
 // class StructureComparer_O : public core::T_O
 // class Structure_Old_ListEntry_O : public core::T_O
 // class Structure_Old_List_O : public core::T_O
 // class SuperposableConformationCollection_O : public ConformationCollection_O
 // class SuperposeEngine_O : public core::T_O
 // class SuperposeSelectedAtoms_O : public SuperposeEngine_O
 // class Topology_O : public core::T_O
 // class TrajectoryFrame_O : public core::T_O
 // class Trajectory_O : public core::T_O
 // class Twister_O : public core::T_O
 // class TwisterDriver_O : public core::T_O
 // class OneTypeRule_O : public core::T_O
 // class TypeAssignmentRules_O : public core::T_O
 // class VirtualAtom_O : public Atom_O
 // class VirtualSphere_O : public core::T_O
 // class ZMatrixInternal_O : public core::T_O
 // class ZMatrixBondInternal_O : public ZMatrixInternal_O
 // class ZMatrixAngleInternal_O : public ZMatrixInternal_O
 // class ZMatrixDihedralInternal_O : public ZMatrixInternal_O
 // class ZMatrixEntry_O : public core::T_O
 // class ZMatrix_O : public core::T_O
#ifdef HEADER_INCLUDES
#include "include/ringFinder.h"
#include "include/largeSquareMatrix.h"
#include "include/alias.h"
#include "include/angle.h"
#include "include/atomIdMap.h"
#include "include/atomId.h"
#include "include/atomIndexer.h"
#include "include/atomReference.h"
#include "include/energyAtomTable.h"
#include "include/bond.h"
#include "include/chemdraw.h"
#include "include/calculatePosition.h"
#include "include/candoDatabase.h"
#include "include/chemInfo.h"
#include "include/cipPrioritizer.h"
#include "include/command.h"
#include "include/complexRestraints.h"
#include "include/conformationCollection.h"
#include "include/conformationExplorer.h"
#include "include/constitutionAtoms.h"
#include "include/coordSys.h"
#include "include/coupling.h"
#include "include/energyComponent.h"
#include "include/energyFunction.h"
#include "include/entity.h"
#include "include/ffBaseDb.h"
#include "include/ffNonbondDb.h"
#include "include/forceField.h"
#include "include/frameRecognizer.h"
#include "include/improperTorsion.h"
#include "include/macroModel.h"
#include "include/matter.h"
#include "include/minimizerLog.h"
#include "include/minimizer.h"
#include "include/monomerContext.h"
#include "include/monomerCoordinates.h"
#include "include/nVector.h"
#include "include/numericalFunction.h"
#include "include/candoScript.h"
#include "include/oligomer.h"
#include "include/monomer.h"
#include "include/typeAssignmentRules.h"
#include "include/pdbMonomer.h"
#include "include/pdb.h"
#include "include/plug.h"
#include "include/properTorsion.h"
#include "include/readAmberParameters.h"
#include "include/restraint.h"
#include "include/rotamer.h"
#include "include/searchStatistics.h"
#include "include/spanningLoop.h"
#include "include/specificContext.h"
#include "include/statusTracker.h"
#include "include/stereochemistry.h"
#include "include/stereoisomerAtoms.h"
#include "include/structureComparer.h"
#include "include/structureList.h"
#include "include/superposeEngine.h"
#include "include/topology.h"
#include "include/trajectory.h"
#include "include/twister.h"
#include "include/virtualSphere.h"
#include "include/zMatrix.h"
#include "include/aggregate.h"
#include "include/atom.h"
#include "include/constitution.h"
#include "include/energyAnchorRestraint.h"
#include "include/energyAngle.h"
#include "include/energyChiralRestraint.h"
#include "include/energyDihedral.h"
#include "include/energyFixedNonbond.h"
#include "include/energyImproperRestraint.h"
#include "include/energyNonbond.h"
#include "include/energyStretch.h"
#include "include/entityNameSet.h"
#include "include/ffAngleDb.h"
#include "include/ffItorDb.h"
#include "include/ffPtorDb.h"
#include "include/ffStretchDb.h"
#include "include/ffTypesDb.h"
#include "include/ffVdwDb.h"
#include "include/iterateMatter.h"
#include "include/iterateRestraints.h"
#include "include/molecule.h"
#include "include/residue.h"
#include "include/superposableConformationCollection.h"
#include "include/virtualAtom.h"
#include "include/monomerPack.h"
#include "include/representedEntityNameSet.h"
#endif // HEADER_INCLUDES
#undef HEADER_INCLUDES
#if defined(SET_SYMBOL) || defined(ALL_STAGES)
// requires LOOKUP_SYMBOL(pkg,symbolName) be defined
chem::AGEdge_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AGEdge_O::static_packageName(),chem::AGEdge_O::static_className()));
chem::AGVertex_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AGVertex_O::static_packageName(),chem::AGVertex_O::static_className()));
chem::AbstractLargeSquareMatrix_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AbstractLargeSquareMatrix_O::static_packageName(),chem::AbstractLargeSquareMatrix_O::static_className()));
chem::Alias_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Alias_O::static_packageName(),chem::Alias_O::static_className()));
chem::Angle_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Angle_O::static_packageName(),chem::Angle_O::static_className()));
chem::AtomIdToAtomMap_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomIdToAtomMap_O::static_packageName(),chem::AtomIdToAtomMap_O::static_className()));
chem::AtomId_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomId_O::static_packageName(),chem::AtomId_O::static_className()));
chem::AtomIndexer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomIndexer_O::static_packageName(),chem::AtomIndexer_O::static_className()));
chem::AtomReferenceBase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomReferenceBase_O::static_packageName(),chem::AtomReferenceBase_O::static_className()));
chem::AtomTable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomTable_O::static_packageName(),chem::AtomTable_O::static_className()));
chem::BondList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::BondList_O::static_packageName(),chem::BondList_O::static_className()));
chem::Bond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Bond_O::static_packageName(),chem::Bond_O::static_className()));
chem::CDBond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CDBond_O::static_packageName(),chem::CDBond_O::static_className()));
chem::CDFragment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CDFragment_O::static_packageName(),chem::CDFragment_O::static_className()));
chem::CDNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CDNode_O::static_packageName(),chem::CDNode_O::static_className()));
chem::CDText_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CDText_O::static_packageName(),chem::CDText_O::static_className()));
chem::CalculatePosition_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CalculatePosition_O::static_packageName(),chem::CalculatePosition_O::static_className()));
chem::CandoDatabase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CandoDatabase_O::static_packageName(),chem::CandoDatabase_O::static_className()));
chem::ChemDraw_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ChemDraw_O::static_packageName(),chem::ChemDraw_O::static_className()));
chem::ChemInfoMatch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ChemInfoMatch_O::static_packageName(),chem::ChemInfoMatch_O::static_className()));
chem::ChemInfoNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ChemInfoNode_O::static_packageName(),chem::ChemInfoNode_O::static_className()));
chem::ChemInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ChemInfo_O::static_packageName(),chem::ChemInfo_O::static_className()));
chem::CipPrioritizer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CipPrioritizer_O::static_packageName(),chem::CipPrioritizer_O::static_className()));
chem::Command_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Command_O::static_packageName(),chem::Command_O::static_className()));
chem::ComplexRestraint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ComplexRestraint_O::static_packageName(),chem::ComplexRestraint_O::static_className()));
chem::ConformationCollectionEntry_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConformationCollectionEntry_O::static_packageName(),chem::ConformationCollectionEntry_O::static_className()));
chem::ConformationCollection_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConformationCollection_O::static_packageName(),chem::ConformationCollection_O::static_className()));
chem::ConformationExplorerEntryStage_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConformationExplorerEntryStage_O::static_packageName(),chem::ConformationExplorerEntryStage_O::static_className()));
chem::ConformationExplorerEntry_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConformationExplorerEntry_O::static_packageName(),chem::ConformationExplorerEntry_O::static_className()));
chem::ConformationExplorerMatch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConformationExplorerMatch_O::static_packageName(),chem::ConformationExplorerMatch_O::static_className()));
chem::ConformationExplorer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConformationExplorer_O::static_packageName(),chem::ConformationExplorer_O::static_className()));
chem::ConstitutionAtom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConstitutionAtom_O::static_packageName(),chem::ConstitutionAtom_O::static_className()));
chem::ConstitutionAtoms_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConstitutionAtoms_O::static_packageName(),chem::ConstitutionAtoms_O::static_className()));
chem::ConstitutionBond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConstitutionBond_O::static_packageName(),chem::ConstitutionBond_O::static_className()));
chem::CoordinateSystem_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CoordinateSystem_O::static_packageName(),chem::CoordinateSystem_O::static_className()));
chem::Coupling_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Coupling_O::static_packageName(),chem::Coupling_O::static_className()));
chem::EnergyComponent_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyComponent_O::static_packageName(),chem::EnergyComponent_O::static_className()));
chem::EnergyFunction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyFunction_O::static_packageName(),chem::EnergyFunction_O::static_className()));
chem::Entity_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Entity_O::static_packageName(),chem::Entity_O::static_className()));
chem::FFBaseDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFBaseDb_O::static_packageName(),chem::FFBaseDb_O::static_className()));
chem::FFNonbondCrossTermTable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFNonbondCrossTermTable_O::static_packageName(),chem::FFNonbondCrossTermTable_O::static_className()));
chem::FFParameter_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFParameter_O::static_packageName(),chem::FFParameter_O::static_className()));
chem::ForceField_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ForceField_O::static_packageName(),chem::ForceField_O::static_className()));
chem::FrameRecognizer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FrameRecognizer_O::static_packageName(),chem::FrameRecognizer_O::static_className()));
chem::ImproperTorsion_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ImproperTorsion_O::static_packageName(),chem::ImproperTorsion_O::static_className()));
chem::InfoDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::InfoDb_O::static_packageName(),chem::InfoDb_O::static_className()));
chem::MacroModelFile_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MacroModelFile_O::static_packageName(),chem::MacroModelFile_O::static_className()));
chem::MapOfMonomerNamesToAtomIndexers_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MapOfMonomerNamesToAtomIndexers_O::static_packageName(),chem::MapOfMonomerNamesToAtomIndexers_O::static_className()));
chem::Matter_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Matter_O::static_packageName(),chem::Matter_O::static_className()));
chem::MinimizerLog_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MinimizerLog_O::static_packageName(),chem::MinimizerLog_O::static_className()));
chem::Minimizer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Minimizer_O::static_packageName(),chem::Minimizer_O::static_className()));
chem::MonomerContext_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MonomerContext_O::static_packageName(),chem::MonomerContext_O::static_className()));
chem::MonomerCoordinates_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MonomerCoordinates_O::static_packageName(),chem::MonomerCoordinates_O::static_className()));
chem::NVector_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::NVector_O::static_packageName(),chem::NVector_O::static_className()));
chem::NumericalFunction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::NumericalFunction_O::static_packageName(),chem::NumericalFunction_O::static_className()));
chem::OligomerPart_Base_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::OligomerPart_Base_O::static_packageName(),chem::OligomerPart_Base_O::static_className()));
chem::Oligomer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Oligomer_O::static_packageName(),chem::Oligomer_O::static_className()));
chem::OneMonomer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::OneMonomer_O::static_packageName(),chem::OneMonomer_O::static_className()));
chem::OneTypeRule_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::OneTypeRule_O::static_packageName(),chem::OneTypeRule_O::static_className()));
chem::PathMessage_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::PathMessage_O::static_packageName(),chem::PathMessage_O::static_className()));
chem::PdbMonomerConnectivity_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::PdbMonomerConnectivity_O::static_packageName(),chem::PdbMonomerConnectivity_O::static_className()));
chem::PdbMonomerDatabase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::PdbMonomerDatabase_O::static_packageName(),chem::PdbMonomerDatabase_O::static_className()));
chem::PdbReader_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::PdbReader_O::static_packageName(),chem::PdbReader_O::static_className()));
chem::PdbWriter_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::PdbWriter_O::static_packageName(),chem::PdbWriter_O::static_className()));
chem::Plug_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Plug_O::static_packageName(),chem::Plug_O::static_className()));
chem::ProperTorsion_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ProperTorsion_O::static_packageName(),chem::ProperTorsion_O::static_className()));
chem::ReadAmberParameters_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ReadAmberParameters_O::static_packageName(),chem::ReadAmberParameters_O::static_className()));
chem::ReportBase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ReportBase_O::static_packageName(),chem::ReportBase_O::static_className()));
chem::ResidueList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ResidueList_O::static_packageName(),chem::ResidueList_O::static_className()));
chem::RestraintList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintList_O::static_packageName(),chem::RestraintList_O::static_className()));
chem::Restraint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Restraint_O::static_packageName(),chem::Restraint_O::static_className()));
chem::RingFinder_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RingFinder_O::static_packageName(),chem::RingFinder_O::static_className()));
chem::Rotamer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Rotamer_O::static_packageName(),chem::Rotamer_O::static_className()));
chem::SearchStatistics_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SearchStatistics_O::static_packageName(),chem::SearchStatistics_O::static_className()));
chem::SpanningLoop_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SpanningLoop_O::static_packageName(),chem::SpanningLoop_O::static_className()));
chem::SpecificContextSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SpecificContextSet_O::static_packageName(),chem::SpecificContextSet_O::static_className()));
chem::SpecificContext_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SpecificContext_O::static_packageName(),chem::SpecificContext_O::static_className()));
chem::StatusTracker_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StatusTracker_O::static_packageName(),chem::StatusTracker_O::static_className()));
chem::StereoConfiguration_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StereoConfiguration_O::static_packageName(),chem::StereoConfiguration_O::static_className()));
chem::StereoInformation_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StereoInformation_O::static_packageName(),chem::StereoInformation_O::static_className()));
chem::StereoisomerAtom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StereoisomerAtom_O::static_packageName(),chem::StereoisomerAtom_O::static_className()));
chem::StereoisomerAtoms_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StereoisomerAtoms_O::static_packageName(),chem::StereoisomerAtoms_O::static_className()));
chem::StructureComparer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StructureComparer_O::static_packageName(),chem::StructureComparer_O::static_className()));
chem::Structure_Old_ListEntry_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Structure_Old_ListEntry_O::static_packageName(),chem::Structure_Old_ListEntry_O::static_className()));
chem::Structure_Old_List_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Structure_Old_List_O::static_packageName(),chem::Structure_Old_List_O::static_className()));
chem::SuperposeEngine_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SuperposeEngine_O::static_packageName(),chem::SuperposeEngine_O::static_className()));
chem::Topology_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Topology_O::static_packageName(),chem::Topology_O::static_className()));
chem::TrajectoryFrame_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::TrajectoryFrame_O::static_packageName(),chem::TrajectoryFrame_O::static_className()));
chem::Trajectory_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Trajectory_O::static_packageName(),chem::Trajectory_O::static_className()));
chem::TwisterDriver_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::TwisterDriver_O::static_packageName(),chem::TwisterDriver_O::static_className()));
chem::Twister_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Twister_O::static_packageName(),chem::Twister_O::static_className()));
chem::TypeAssignmentRules_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::TypeAssignmentRules_O::static_packageName(),chem::TypeAssignmentRules_O::static_className()));
chem::VirtualSphere_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::VirtualSphere_O::static_packageName(),chem::VirtualSphere_O::static_className()));
chem::WildElementDict_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::WildElementDict_O::static_packageName(),chem::WildElementDict_O::static_className()));
chem::ZMatrixEntry_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ZMatrixEntry_O::static_packageName(),chem::ZMatrixEntry_O::static_className()));
chem::ZMatrixInternal_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ZMatrixInternal_O::static_packageName(),chem::ZMatrixInternal_O::static_className()));
chem::ZMatrix_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ZMatrix_O::static_packageName(),chem::ZMatrix_O::static_className()));
chem::Aggregate_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Aggregate_O::static_packageName(),chem::Aggregate_O::static_className()));
chem::AtomReference_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomReference_O::static_packageName(),chem::AtomReference_O::static_className()));
chem::Atom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Atom_O::static_packageName(),chem::Atom_O::static_className()));
chem::BondListMatchNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::BondListMatchNode_O::static_packageName(),chem::BondListMatchNode_O::static_className()));
chem::BondMatchNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::BondMatchNode_O::static_packageName(),chem::BondMatchNode_O::static_className()));
chem::CalculatePositionAlongBond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CalculatePositionAlongBond_O::static_packageName(),chem::CalculatePositionAlongBond_O::static_className()));
chem::CalculatePositionRelativeToOrigin_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CalculatePositionRelativeToOrigin_O::static_packageName(),chem::CalculatePositionRelativeToOrigin_O::static_className()));
chem::CalculatePositionUsingInternals_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::CalculatePositionUsingInternals_O::static_packageName(),chem::CalculatePositionUsingInternals_O::static_className()));
chem::ConstitutionVirtualAtom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ConstitutionVirtualAtom_O::static_packageName(),chem::ConstitutionVirtualAtom_O::static_className()));
chem::Constitution_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Constitution_O::static_packageName(),chem::Constitution_O::static_className()));
chem::DirectionalCoupling_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::DirectionalCoupling_O::static_packageName(),chem::DirectionalCoupling_O::static_className()));
chem::EnergyAnchorRestraint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyAnchorRestraint_O::static_packageName(),chem::EnergyAnchorRestraint_O::static_className()));
chem::EnergyAngle_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyAngle_O::static_packageName(),chem::EnergyAngle_O::static_className()));
chem::EnergyChiralRestraint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyChiralRestraint_O::static_packageName(),chem::EnergyChiralRestraint_O::static_className()));
chem::EnergyDihedral_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyDihedral_O::static_packageName(),chem::EnergyDihedral_O::static_className()));
chem::EnergyFixedNonbondRestraint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyFixedNonbondRestraint_O::static_packageName(),chem::EnergyFixedNonbondRestraint_O::static_className()));
chem::EnergyImproperRestraint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyImproperRestraint_O::static_packageName(),chem::EnergyImproperRestraint_O::static_className()));
chem::EnergyNonbond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyNonbond_O::static_packageName(),chem::EnergyNonbond_O::static_className()));
chem::EnergyStretch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EnergyStretch_O::static_packageName(),chem::EnergyStretch_O::static_className()));
chem::EntityNameSetBase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EntityNameSetBase_O::static_packageName(),chem::EntityNameSetBase_O::static_className()));
chem::FFAngleDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFAngleDb_O::static_packageName(),chem::FFAngleDb_O::static_className()));
chem::FFAngle_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFAngle_O::static_packageName(),chem::FFAngle_O::static_className()));
chem::FFItorDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFItorDb_O::static_packageName(),chem::FFItorDb_O::static_className()));
chem::FFItor_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFItor_O::static_packageName(),chem::FFItor_O::static_className()));
chem::FFNonbondDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFNonbondDb_O::static_packageName(),chem::FFNonbondDb_O::static_className()));
chem::FFNonbond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFNonbond_O::static_packageName(),chem::FFNonbond_O::static_className()));
chem::FFPtorDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFPtorDb_O::static_packageName(),chem::FFPtorDb_O::static_className()));
chem::FFPtor_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFPtor_O::static_packageName(),chem::FFPtor_O::static_className()));
chem::FFStretchDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFStretchDb_O::static_packageName(),chem::FFStretchDb_O::static_className()));
chem::FFStretch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFStretch_O::static_packageName(),chem::FFStretch_O::static_className()));
chem::FFTypesDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFTypesDb_O::static_packageName(),chem::FFTypesDb_O::static_className()));
chem::FFVdwDb_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFVdwDb_O::static_packageName(),chem::FFVdwDb_O::static_className()));
chem::FFVdw_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FFVdw_O::static_packageName(),chem::FFVdw_O::static_className()));
chem::ForceMatchReport_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ForceMatchReport_O::static_packageName(),chem::ForceMatchReport_O::static_className()));
chem::FullLargeSquareMatrix_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::FullLargeSquareMatrix_O::static_packageName(),chem::FullLargeSquareMatrix_O::static_className()));
chem::IterateMatter_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::IterateMatter_O::static_packageName(),chem::IterateMatter_O::static_className()));
chem::IterateRestraints_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::IterateRestraints_O::static_packageName(),chem::IterateRestraints_O::static_className()));
chem::JumpPlug_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::JumpPlug_O::static_packageName(),chem::JumpPlug_O::static_className()));
chem::MessageReport_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MessageReport_O::static_packageName(),chem::MessageReport_O::static_className()));
chem::Molecule_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Molecule_O::static_packageName(),chem::Molecule_O::static_className()));
chem::Monomer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Monomer_O::static_packageName(),chem::Monomer_O::static_className()));
chem::OligomerPart_Link_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::OligomerPart_Link_O::static_packageName(),chem::OligomerPart_Link_O::static_className()));
chem::OligomerPart_Monomer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::OligomerPart_Monomer_O::static_packageName(),chem::OligomerPart_Monomer_O::static_className()));
chem::PlugWithMates_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::PlugWithMates_O::static_packageName(),chem::PlugWithMates_O::static_className()));
chem::Residue_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Residue_O::static_packageName(),chem::Residue_O::static_className()));
chem::RestrainedExoCyclicAtom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestrainedExoCyclicAtom_O::static_packageName(),chem::RestrainedExoCyclicAtom_O::static_className()));
chem::RestrainedPiBond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestrainedPiBond_O::static_packageName(),chem::RestrainedPiBond_O::static_className()));
chem::RestraintAnchor_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintAnchor_O::static_packageName(),chem::RestraintAnchor_O::static_className()));
chem::RestraintAngle_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintAngle_O::static_packageName(),chem::RestraintAngle_O::static_className()));
chem::RestraintChiral_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintChiral_O::static_packageName(),chem::RestraintChiral_O::static_className()));
chem::RestraintDihedral_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintDihedral_O::static_packageName(),chem::RestraintDihedral_O::static_className()));
chem::RestraintDistance_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintDistance_O::static_packageName(),chem::RestraintDistance_O::static_className()));
chem::RestraintFixedNonbond_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RestraintFixedNonbond_O::static_packageName(),chem::RestraintFixedNonbond_O::static_className()));
chem::RingCoupling_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RingCoupling_O::static_packageName(),chem::RingCoupling_O::static_className()));
chem::RootMatchNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RootMatchNode_O::static_packageName(),chem::RootMatchNode_O::static_className()));
chem::SparseLargeSquareMatrix_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SparseLargeSquareMatrix_O::static_packageName(),chem::SparseLargeSquareMatrix_O::static_className()));
chem::StepReport_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StepReport_O::static_packageName(),chem::StepReport_O::static_className()));
chem::StereoisomerVirtualAtom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::StereoisomerVirtualAtom_O::static_packageName(),chem::StereoisomerVirtualAtom_O::static_className()));
chem::Stereoisomer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Stereoisomer_O::static_packageName(),chem::Stereoisomer_O::static_className()));
chem::SuperposableConformationCollection_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SuperposableConformationCollection_O::static_packageName(),chem::SuperposableConformationCollection_O::static_className()));
chem::SuperposeSelectedAtoms_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SuperposeSelectedAtoms_O::static_packageName(),chem::SuperposeSelectedAtoms_O::static_className()));
chem::ZMatrixAngleInternal_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ZMatrixAngleInternal_O::static_packageName(),chem::ZMatrixAngleInternal_O::static_className()));
chem::ZMatrixBondInternal_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ZMatrixBondInternal_O::static_packageName(),chem::ZMatrixBondInternal_O::static_className()));
chem::ZMatrixDihedralInternal_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ZMatrixDihedralInternal_O::static_packageName(),chem::ZMatrixDihedralInternal_O::static_className()));
chem::AfterMatchBondTest_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AfterMatchBondTest_O::static_packageName(),chem::AfterMatchBondTest_O::static_className()));
chem::AntechamberBondTest_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AntechamberBondTest_O::static_packageName(),chem::AntechamberBondTest_O::static_className()));
chem::AtomOrBondMatchNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomOrBondMatchNode_O::static_packageName(),chem::AtomOrBondMatchNode_O::static_className()));
chem::BondTest_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::BondTest_O::static_packageName(),chem::BondTest_O::static_className()));
chem::Branch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Branch_O::static_packageName(),chem::Branch_O::static_className()));
chem::Chain_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Chain_O::static_packageName(),chem::Chain_O::static_className()));
chem::EntityNameSetWithCap_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EntityNameSetWithCap_O::static_packageName(),chem::EntityNameSetWithCap_O::static_className()));
chem::EntityNameSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::EntityNameSet_O::static_packageName(),chem::EntityNameSet_O::static_className()));
chem::InPlug_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::InPlug_O::static_packageName(),chem::InPlug_O::static_className()));
chem::IterateAtoms_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::IterateAtoms_O::static_packageName(),chem::IterateAtoms_O::static_className()));
chem::IterateBonds_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::IterateBonds_O::static_packageName(),chem::IterateBonds_O::static_className()));
chem::IterateResidues_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::IterateResidues_O::static_packageName(),chem::IterateResidues_O::static_className()));
chem::MultiMonomer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MultiMonomer_O::static_packageName(),chem::MultiMonomer_O::static_className()));
chem::OutPlug_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::OutPlug_O::static_packageName(),chem::OutPlug_O::static_className()));
chem::VirtualAtom_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::VirtualAtom_O::static_packageName(),chem::VirtualAtom_O::static_className()));
chem::AntechamberFocusAtomMatch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AntechamberFocusAtomMatch_O::static_packageName(),chem::AntechamberFocusAtomMatch_O::static_className()));
chem::AtomTest_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AtomTest_O::static_packageName(),chem::AtomTest_O::static_className()));
chem::Logical_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Logical_O::static_packageName(),chem::Logical_O::static_className()));
chem::Mate_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Mate_O::static_packageName(),chem::Mate_O::static_className()));
chem::MonomerPack_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::MonomerPack_O::static_packageName(),chem::MonomerPack_O::static_className()));
chem::RepresentedEntityNameSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RepresentedEntityNameSet_O::static_packageName(),chem::RepresentedEntityNameSet_O::static_className()));
chem::ResidueTest_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::ResidueTest_O::static_packageName(),chem::ResidueTest_O::static_className()));
chem::RingClosingMate_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RingClosingMate_O::static_packageName(),chem::RingClosingMate_O::static_className()));
chem::RingClosingPlug_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RingClosingPlug_O::static_packageName(),chem::RingClosingPlug_O::static_className()));
chem::RingTest_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::RingTest_O::static_packageName(),chem::RingTest_O::static_className()));
chem::Root_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::Root_O::static_packageName(),chem::Root_O::static_className()));
chem::TagSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::TagSet_O::static_packageName(),chem::TagSet_O::static_className()));
chem::AntechamberRoot_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::AntechamberRoot_O::static_packageName(),chem::AntechamberRoot_O::static_className()));
chem::SmartsRoot_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(chem::SmartsRoot_O::static_packageName(),chem::SmartsRoot_O::static_className()));
#endif // SET_SYMBOL
#undef SET_SYMBOL
#if defined(CREATE_CLASS) || defined(ALL_STAGES)

    LOG(BF("Creating class[classchem__AGEdge_Oval]"));
    core::BuiltInClass_sp classchem__AGEdge_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AGEdge_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AGEdge_Oval,_lisp,chem::AGEdge_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AGEdge_O>::id,chem::AGEdge_O::static_classSymbol());
    chem::AGEdge_O::___staticClass = classchem__AGEdge_Oval;
#ifdef USE_MPS
    chem::AGEdge_O::static_Kind = gctools::GCKind<chem::AGEdge_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AGEdge_Oval,chem::AGEdge_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AGEdge_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AGEdge_O>>::allocateClass();
        chem::AGEdge_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AGEdge_O::static_className() % (void*)(chem::AGEdge_O::static_allocator) );
    classchem__AGEdge_Oval->setCreator(chem::AGEdge_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AGEdge_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AGEdge_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AGVertex_Oval]"));
    core::BuiltInClass_sp classchem__AGVertex_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AGVertex_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AGVertex_Oval,_lisp,chem::AGVertex_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AGVertex_O>::id,chem::AGVertex_O::static_classSymbol());
    chem::AGVertex_O::___staticClass = classchem__AGVertex_Oval;
#ifdef USE_MPS
    chem::AGVertex_O::static_Kind = gctools::GCKind<chem::AGVertex_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AGVertex_Oval,chem::AGVertex_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AGVertex_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AGVertex_O>>::allocateClass();
        chem::AGVertex_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AGVertex_O::static_className() % (void*)(chem::AGVertex_O::static_allocator) );
    classchem__AGVertex_Oval->setCreator(chem::AGVertex_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AGVertex_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AGVertex_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AbstractLargeSquareMatrix_Oval]"));
    core::BuiltInClass_sp classchem__AbstractLargeSquareMatrix_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AbstractLargeSquareMatrix_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AbstractLargeSquareMatrix_Oval,_lisp,chem::AbstractLargeSquareMatrix_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AbstractLargeSquareMatrix_O>::id,chem::AbstractLargeSquareMatrix_O::static_classSymbol());
    chem::AbstractLargeSquareMatrix_O::___staticClass = classchem__AbstractLargeSquareMatrix_Oval;
#ifdef USE_MPS
    chem::AbstractLargeSquareMatrix_O::static_Kind = gctools::GCKind<chem::AbstractLargeSquareMatrix_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AbstractLargeSquareMatrix_Oval,chem::AbstractLargeSquareMatrix_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AbstractLargeSquareMatrix_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AbstractLargeSquareMatrix_O>>::allocateClass();
        chem::AbstractLargeSquareMatrix_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AbstractLargeSquareMatrix_O::static_className() % (void*)(chem::AbstractLargeSquareMatrix_O::static_allocator) );
    classchem__AbstractLargeSquareMatrix_Oval->setCreator(chem::AbstractLargeSquareMatrix_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AbstractLargeSquareMatrix_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AbstractLargeSquareMatrix_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Alias_Oval]"));
    core::BuiltInClass_sp classchem__Alias_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Alias_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Alias_Oval,_lisp,chem::Alias_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Alias_O>::id,chem::Alias_O::static_classSymbol());
    chem::Alias_O::___staticClass = classchem__Alias_Oval;
#ifdef USE_MPS
    chem::Alias_O::static_Kind = gctools::GCKind<chem::Alias_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Alias_Oval,chem::Alias_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Alias_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Alias_O>>::allocateClass();
        chem::Alias_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Alias_O::static_className() % (void*)(chem::Alias_O::static_allocator) );
    classchem__Alias_Oval->setCreator(chem::Alias_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Alias_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Alias_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Angle_Oval]"));
    core::BuiltInClass_sp classchem__Angle_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Angle_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Angle_Oval,_lisp,chem::Angle_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Angle_O>::id,chem::Angle_O::static_classSymbol());
    chem::Angle_O::___staticClass = classchem__Angle_Oval;
#ifdef USE_MPS
    chem::Angle_O::static_Kind = gctools::GCKind<chem::Angle_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Angle_Oval,chem::Angle_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Angle_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Angle_O>>::allocateClass();
        chem::Angle_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Angle_O::static_className() % (void*)(chem::Angle_O::static_allocator) );
    classchem__Angle_Oval->setCreator(chem::Angle_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Angle_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Angle_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomIdToAtomMap_Oval]"));
    core::BuiltInClass_sp classchem__AtomIdToAtomMap_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomIdToAtomMap_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomIdToAtomMap_Oval,_lisp,chem::AtomIdToAtomMap_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomIdToAtomMap_O>::id,chem::AtomIdToAtomMap_O::static_classSymbol());
    chem::AtomIdToAtomMap_O::___staticClass = classchem__AtomIdToAtomMap_Oval;
#ifdef USE_MPS
    chem::AtomIdToAtomMap_O::static_Kind = gctools::GCKind<chem::AtomIdToAtomMap_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomIdToAtomMap_Oval,chem::AtomIdToAtomMap_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomIdToAtomMap_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomIdToAtomMap_O>>::allocateClass();
        chem::AtomIdToAtomMap_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomIdToAtomMap_O::static_className() % (void*)(chem::AtomIdToAtomMap_O::static_allocator) );
    classchem__AtomIdToAtomMap_Oval->setCreator(chem::AtomIdToAtomMap_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomIdToAtomMap_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomIdToAtomMap_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomId_Oval]"));
    core::BuiltInClass_sp classchem__AtomId_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomId_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomId_Oval,_lisp,chem::AtomId_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomId_O>::id,chem::AtomId_O::static_classSymbol());
    chem::AtomId_O::___staticClass = classchem__AtomId_Oval;
#ifdef USE_MPS
    chem::AtomId_O::static_Kind = gctools::GCKind<chem::AtomId_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomId_Oval,chem::AtomId_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomId_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomId_O>>::allocateClass();
        chem::AtomId_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomId_O::static_className() % (void*)(chem::AtomId_O::static_allocator) );
    classchem__AtomId_Oval->setCreator(chem::AtomId_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomId_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomId_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomIndexer_Oval]"));
    core::BuiltInClass_sp classchem__AtomIndexer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomIndexer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomIndexer_Oval,_lisp,chem::AtomIndexer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomIndexer_O>::id,chem::AtomIndexer_O::static_classSymbol());
    chem::AtomIndexer_O::___staticClass = classchem__AtomIndexer_Oval;
#ifdef USE_MPS
    chem::AtomIndexer_O::static_Kind = gctools::GCKind<chem::AtomIndexer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomIndexer_Oval,chem::AtomIndexer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomIndexer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomIndexer_O>>::allocateClass();
        chem::AtomIndexer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomIndexer_O::static_className() % (void*)(chem::AtomIndexer_O::static_allocator) );
    classchem__AtomIndexer_Oval->setCreator(chem::AtomIndexer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomIndexer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomIndexer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomReferenceBase_Oval]"));
    core::BuiltInClass_sp classchem__AtomReferenceBase_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomReferenceBase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomReferenceBase_Oval,_lisp,chem::AtomReferenceBase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomReferenceBase_O>::id,chem::AtomReferenceBase_O::static_classSymbol());
    chem::AtomReferenceBase_O::___staticClass = classchem__AtomReferenceBase_Oval;
#ifdef USE_MPS
    chem::AtomReferenceBase_O::static_Kind = gctools::GCKind<chem::AtomReferenceBase_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomReferenceBase_Oval,chem::AtomReferenceBase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomReferenceBase_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomReferenceBase_O>>::allocateClass();
        chem::AtomReferenceBase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomReferenceBase_O::static_className() % (void*)(chem::AtomReferenceBase_O::static_allocator) );
    classchem__AtomReferenceBase_Oval->setCreator(chem::AtomReferenceBase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomReferenceBase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomReferenceBase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomTable_Oval]"));
    core::BuiltInClass_sp classchem__AtomTable_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomTable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomTable_Oval,_lisp,chem::AtomTable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomTable_O>::id,chem::AtomTable_O::static_classSymbol());
    chem::AtomTable_O::___staticClass = classchem__AtomTable_Oval;
#ifdef USE_MPS
    chem::AtomTable_O::static_Kind = gctools::GCKind<chem::AtomTable_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomTable_Oval,chem::AtomTable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomTable_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomTable_O>>::allocateClass();
        chem::AtomTable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomTable_O::static_className() % (void*)(chem::AtomTable_O::static_allocator) );
    classchem__AtomTable_Oval->setCreator(chem::AtomTable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomTable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomTable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__BondList_Oval]"));
    core::BuiltInClass_sp classchem__BondList_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__BondList_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__BondList_Oval,_lisp,chem::BondList_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::BondList_O>::id,chem::BondList_O::static_classSymbol());
    chem::BondList_O::___staticClass = classchem__BondList_Oval;
#ifdef USE_MPS
    chem::BondList_O::static_Kind = gctools::GCKind<chem::BondList_O>::Kind;
#endif
    core::af_setf_findClass(classchem__BondList_Oval,chem::BondList_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::BondList_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::BondList_O>>::allocateClass();
        chem::BondList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::BondList_O::static_className() % (void*)(chem::BondList_O::static_allocator) );
    classchem__BondList_Oval->setCreator(chem::BondList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::BondList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__BondList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Bond_Oval]"));
    core::BuiltInClass_sp classchem__Bond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Bond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Bond_Oval,_lisp,chem::Bond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Bond_O>::id,chem::Bond_O::static_classSymbol());
    chem::Bond_O::___staticClass = classchem__Bond_Oval;
#ifdef USE_MPS
    chem::Bond_O::static_Kind = gctools::GCKind<chem::Bond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Bond_Oval,chem::Bond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Bond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Bond_O>>::allocateClass();
        chem::Bond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Bond_O::static_className() % (void*)(chem::Bond_O::static_allocator) );
    classchem__Bond_Oval->setCreator(chem::Bond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Bond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Bond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CDBond_Oval]"));
    core::BuiltInClass_sp classchem__CDBond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CDBond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CDBond_Oval,_lisp,chem::CDBond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CDBond_O>::id,chem::CDBond_O::static_classSymbol());
    chem::CDBond_O::___staticClass = classchem__CDBond_Oval;
#ifdef USE_MPS
    chem::CDBond_O::static_Kind = gctools::GCKind<chem::CDBond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CDBond_Oval,chem::CDBond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CDBond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CDBond_O>>::allocateClass();
        chem::CDBond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CDBond_O::static_className() % (void*)(chem::CDBond_O::static_allocator) );
    classchem__CDBond_Oval->setCreator(chem::CDBond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CDBond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CDBond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CDFragment_Oval]"));
    core::BuiltInClass_sp classchem__CDFragment_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CDFragment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CDFragment_Oval,_lisp,chem::CDFragment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CDFragment_O>::id,chem::CDFragment_O::static_classSymbol());
    chem::CDFragment_O::___staticClass = classchem__CDFragment_Oval;
#ifdef USE_MPS
    chem::CDFragment_O::static_Kind = gctools::GCKind<chem::CDFragment_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CDFragment_Oval,chem::CDFragment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CDFragment_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CDFragment_O>>::allocateClass();
        chem::CDFragment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CDFragment_O::static_className() % (void*)(chem::CDFragment_O::static_allocator) );
    classchem__CDFragment_Oval->setCreator(chem::CDFragment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CDFragment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CDFragment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CDNode_Oval]"));
    core::BuiltInClass_sp classchem__CDNode_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CDNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CDNode_Oval,_lisp,chem::CDNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CDNode_O>::id,chem::CDNode_O::static_classSymbol());
    chem::CDNode_O::___staticClass = classchem__CDNode_Oval;
#ifdef USE_MPS
    chem::CDNode_O::static_Kind = gctools::GCKind<chem::CDNode_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CDNode_Oval,chem::CDNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CDNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CDNode_O>>::allocateClass();
        chem::CDNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CDNode_O::static_className() % (void*)(chem::CDNode_O::static_allocator) );
    classchem__CDNode_Oval->setCreator(chem::CDNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CDNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CDNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CDText_Oval]"));
    core::BuiltInClass_sp classchem__CDText_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CDText_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CDText_Oval,_lisp,chem::CDText_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CDText_O>::id,chem::CDText_O::static_classSymbol());
    chem::CDText_O::___staticClass = classchem__CDText_Oval;
#ifdef USE_MPS
    chem::CDText_O::static_Kind = gctools::GCKind<chem::CDText_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CDText_Oval,chem::CDText_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CDText_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CDText_O>>::allocateClass();
        chem::CDText_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CDText_O::static_className() % (void*)(chem::CDText_O::static_allocator) );
    classchem__CDText_Oval->setCreator(chem::CDText_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CDText_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CDText_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CalculatePosition_Oval]"));
    core::BuiltInClass_sp classchem__CalculatePosition_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CalculatePosition_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CalculatePosition_Oval,_lisp,chem::CalculatePosition_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CalculatePosition_O>::id,chem::CalculatePosition_O::static_classSymbol());
    chem::CalculatePosition_O::___staticClass = classchem__CalculatePosition_Oval;
#ifdef USE_MPS
    chem::CalculatePosition_O::static_Kind = gctools::GCKind<chem::CalculatePosition_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CalculatePosition_Oval,chem::CalculatePosition_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CalculatePosition_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CalculatePosition_O>>::allocateClass();
        chem::CalculatePosition_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CalculatePosition_O::static_className() % (void*)(chem::CalculatePosition_O::static_allocator) );
    classchem__CalculatePosition_Oval->setCreator(chem::CalculatePosition_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CalculatePosition_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CalculatePosition_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CandoDatabase_Oval]"));
    core::BuiltInClass_sp classchem__CandoDatabase_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CandoDatabase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CandoDatabase_Oval,_lisp,chem::CandoDatabase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CandoDatabase_O>::id,chem::CandoDatabase_O::static_classSymbol());
    chem::CandoDatabase_O::___staticClass = classchem__CandoDatabase_Oval;
#ifdef USE_MPS
    chem::CandoDatabase_O::static_Kind = gctools::GCKind<chem::CandoDatabase_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CandoDatabase_Oval,chem::CandoDatabase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CandoDatabase_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CandoDatabase_O>>::allocateClass();
        chem::CandoDatabase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CandoDatabase_O::static_className() % (void*)(chem::CandoDatabase_O::static_allocator) );
    classchem__CandoDatabase_Oval->setCreator(chem::CandoDatabase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CandoDatabase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CandoDatabase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ChemDraw_Oval]"));
    core::BuiltInClass_sp classchem__ChemDraw_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ChemDraw_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ChemDraw_Oval,_lisp,chem::ChemDraw_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ChemDraw_O>::id,chem::ChemDraw_O::static_classSymbol());
    chem::ChemDraw_O::___staticClass = classchem__ChemDraw_Oval;
#ifdef USE_MPS
    chem::ChemDraw_O::static_Kind = gctools::GCKind<chem::ChemDraw_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ChemDraw_Oval,chem::ChemDraw_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ChemDraw_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ChemDraw_O>>::allocateClass();
        chem::ChemDraw_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ChemDraw_O::static_className() % (void*)(chem::ChemDraw_O::static_allocator) );
    classchem__ChemDraw_Oval->setCreator(chem::ChemDraw_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ChemDraw_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ChemDraw_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ChemInfoMatch_Oval]"));
    core::BuiltInClass_sp classchem__ChemInfoMatch_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ChemInfoMatch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ChemInfoMatch_Oval,_lisp,chem::ChemInfoMatch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ChemInfoMatch_O>::id,chem::ChemInfoMatch_O::static_classSymbol());
    chem::ChemInfoMatch_O::___staticClass = classchem__ChemInfoMatch_Oval;
#ifdef USE_MPS
    chem::ChemInfoMatch_O::static_Kind = gctools::GCKind<chem::ChemInfoMatch_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ChemInfoMatch_Oval,chem::ChemInfoMatch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ChemInfoMatch_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ChemInfoMatch_O>>::allocateClass();
        chem::ChemInfoMatch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ChemInfoMatch_O::static_className() % (void*)(chem::ChemInfoMatch_O::static_allocator) );
    classchem__ChemInfoMatch_Oval->setCreator(chem::ChemInfoMatch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ChemInfoMatch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ChemInfoMatch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ChemInfoNode_Oval]"));
    core::BuiltInClass_sp classchem__ChemInfoNode_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ChemInfoNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ChemInfoNode_Oval,_lisp,chem::ChemInfoNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ChemInfoNode_O>::id,chem::ChemInfoNode_O::static_classSymbol());
    chem::ChemInfoNode_O::___staticClass = classchem__ChemInfoNode_Oval;
#ifdef USE_MPS
    chem::ChemInfoNode_O::static_Kind = gctools::GCKind<chem::ChemInfoNode_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ChemInfoNode_Oval,chem::ChemInfoNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ChemInfoNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ChemInfoNode_O>>::allocateClass();
        chem::ChemInfoNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ChemInfoNode_O::static_className() % (void*)(chem::ChemInfoNode_O::static_allocator) );
    classchem__ChemInfoNode_Oval->setCreator(chem::ChemInfoNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ChemInfoNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ChemInfoNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ChemInfo_Oval]"));
    core::BuiltInClass_sp classchem__ChemInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ChemInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ChemInfo_Oval,_lisp,chem::ChemInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ChemInfo_O>::id,chem::ChemInfo_O::static_classSymbol());
    chem::ChemInfo_O::___staticClass = classchem__ChemInfo_Oval;
#ifdef USE_MPS
    chem::ChemInfo_O::static_Kind = gctools::GCKind<chem::ChemInfo_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ChemInfo_Oval,chem::ChemInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ChemInfo_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ChemInfo_O>>::allocateClass();
        chem::ChemInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ChemInfo_O::static_className() % (void*)(chem::ChemInfo_O::static_allocator) );
    classchem__ChemInfo_Oval->setCreator(chem::ChemInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ChemInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ChemInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CipPrioritizer_Oval]"));
    core::BuiltInClass_sp classchem__CipPrioritizer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CipPrioritizer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CipPrioritizer_Oval,_lisp,chem::CipPrioritizer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CipPrioritizer_O>::id,chem::CipPrioritizer_O::static_classSymbol());
    chem::CipPrioritizer_O::___staticClass = classchem__CipPrioritizer_Oval;
#ifdef USE_MPS
    chem::CipPrioritizer_O::static_Kind = gctools::GCKind<chem::CipPrioritizer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CipPrioritizer_Oval,chem::CipPrioritizer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CipPrioritizer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CipPrioritizer_O>>::allocateClass();
        chem::CipPrioritizer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CipPrioritizer_O::static_className() % (void*)(chem::CipPrioritizer_O::static_allocator) );
    classchem__CipPrioritizer_Oval->setCreator(chem::CipPrioritizer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CipPrioritizer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CipPrioritizer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Command_Oval]"));
    core::BuiltInClass_sp classchem__Command_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Command_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Command_Oval,_lisp,chem::Command_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Command_O>::id,chem::Command_O::static_classSymbol());
    chem::Command_O::___staticClass = classchem__Command_Oval;
#ifdef USE_MPS
    chem::Command_O::static_Kind = gctools::GCKind<chem::Command_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Command_Oval,chem::Command_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Command_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Command_O>>::allocateClass();
        chem::Command_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Command_O::static_className() % (void*)(chem::Command_O::static_allocator) );
    classchem__Command_Oval->setCreator(chem::Command_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Command_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Command_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ComplexRestraint_Oval]"));
    core::BuiltInClass_sp classchem__ComplexRestraint_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ComplexRestraint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ComplexRestraint_Oval,_lisp,chem::ComplexRestraint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ComplexRestraint_O>::id,chem::ComplexRestraint_O::static_classSymbol());
    chem::ComplexRestraint_O::___staticClass = classchem__ComplexRestraint_Oval;
#ifdef USE_MPS
    chem::ComplexRestraint_O::static_Kind = gctools::GCKind<chem::ComplexRestraint_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ComplexRestraint_Oval,chem::ComplexRestraint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ComplexRestraint_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ComplexRestraint_O>>::allocateClass();
        chem::ComplexRestraint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ComplexRestraint_O::static_className() % (void*)(chem::ComplexRestraint_O::static_allocator) );
    classchem__ComplexRestraint_Oval->setCreator(chem::ComplexRestraint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ComplexRestraint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ComplexRestraint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConformationCollectionEntry_Oval]"));
    core::BuiltInClass_sp classchem__ConformationCollectionEntry_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConformationCollectionEntry_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConformationCollectionEntry_Oval,_lisp,chem::ConformationCollectionEntry_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConformationCollectionEntry_O>::id,chem::ConformationCollectionEntry_O::static_classSymbol());
    chem::ConformationCollectionEntry_O::___staticClass = classchem__ConformationCollectionEntry_Oval;
#ifdef USE_MPS
    chem::ConformationCollectionEntry_O::static_Kind = gctools::GCKind<chem::ConformationCollectionEntry_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConformationCollectionEntry_Oval,chem::ConformationCollectionEntry_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConformationCollectionEntry_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConformationCollectionEntry_O>>::allocateClass();
        chem::ConformationCollectionEntry_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConformationCollectionEntry_O::static_className() % (void*)(chem::ConformationCollectionEntry_O::static_allocator) );
    classchem__ConformationCollectionEntry_Oval->setCreator(chem::ConformationCollectionEntry_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConformationCollectionEntry_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConformationCollectionEntry_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConformationCollection_Oval]"));
    core::BuiltInClass_sp classchem__ConformationCollection_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConformationCollection_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConformationCollection_Oval,_lisp,chem::ConformationCollection_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConformationCollection_O>::id,chem::ConformationCollection_O::static_classSymbol());
    chem::ConformationCollection_O::___staticClass = classchem__ConformationCollection_Oval;
#ifdef USE_MPS
    chem::ConformationCollection_O::static_Kind = gctools::GCKind<chem::ConformationCollection_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConformationCollection_Oval,chem::ConformationCollection_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConformationCollection_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConformationCollection_O>>::allocateClass();
        chem::ConformationCollection_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConformationCollection_O::static_className() % (void*)(chem::ConformationCollection_O::static_allocator) );
    classchem__ConformationCollection_Oval->setCreator(chem::ConformationCollection_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConformationCollection_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConformationCollection_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConformationExplorerEntryStage_Oval]"));
    core::BuiltInClass_sp classchem__ConformationExplorerEntryStage_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConformationExplorerEntryStage_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConformationExplorerEntryStage_Oval,_lisp,chem::ConformationExplorerEntryStage_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConformationExplorerEntryStage_O>::id,chem::ConformationExplorerEntryStage_O::static_classSymbol());
    chem::ConformationExplorerEntryStage_O::___staticClass = classchem__ConformationExplorerEntryStage_Oval;
#ifdef USE_MPS
    chem::ConformationExplorerEntryStage_O::static_Kind = gctools::GCKind<chem::ConformationExplorerEntryStage_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConformationExplorerEntryStage_Oval,chem::ConformationExplorerEntryStage_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConformationExplorerEntryStage_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConformationExplorerEntryStage_O>>::allocateClass();
        chem::ConformationExplorerEntryStage_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConformationExplorerEntryStage_O::static_className() % (void*)(chem::ConformationExplorerEntryStage_O::static_allocator) );
    classchem__ConformationExplorerEntryStage_Oval->setCreator(chem::ConformationExplorerEntryStage_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConformationExplorerEntryStage_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConformationExplorerEntryStage_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConformationExplorerEntry_Oval]"));
    core::BuiltInClass_sp classchem__ConformationExplorerEntry_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConformationExplorerEntry_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConformationExplorerEntry_Oval,_lisp,chem::ConformationExplorerEntry_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConformationExplorerEntry_O>::id,chem::ConformationExplorerEntry_O::static_classSymbol());
    chem::ConformationExplorerEntry_O::___staticClass = classchem__ConformationExplorerEntry_Oval;
#ifdef USE_MPS
    chem::ConformationExplorerEntry_O::static_Kind = gctools::GCKind<chem::ConformationExplorerEntry_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConformationExplorerEntry_Oval,chem::ConformationExplorerEntry_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConformationExplorerEntry_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConformationExplorerEntry_O>>::allocateClass();
        chem::ConformationExplorerEntry_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConformationExplorerEntry_O::static_className() % (void*)(chem::ConformationExplorerEntry_O::static_allocator) );
    classchem__ConformationExplorerEntry_Oval->setCreator(chem::ConformationExplorerEntry_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConformationExplorerEntry_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConformationExplorerEntry_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConformationExplorerMatch_Oval]"));
    core::BuiltInClass_sp classchem__ConformationExplorerMatch_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConformationExplorerMatch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConformationExplorerMatch_Oval,_lisp,chem::ConformationExplorerMatch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConformationExplorerMatch_O>::id,chem::ConformationExplorerMatch_O::static_classSymbol());
    chem::ConformationExplorerMatch_O::___staticClass = classchem__ConformationExplorerMatch_Oval;
#ifdef USE_MPS
    chem::ConformationExplorerMatch_O::static_Kind = gctools::GCKind<chem::ConformationExplorerMatch_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConformationExplorerMatch_Oval,chem::ConformationExplorerMatch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConformationExplorerMatch_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConformationExplorerMatch_O>>::allocateClass();
        chem::ConformationExplorerMatch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConformationExplorerMatch_O::static_className() % (void*)(chem::ConformationExplorerMatch_O::static_allocator) );
    classchem__ConformationExplorerMatch_Oval->setCreator(chem::ConformationExplorerMatch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConformationExplorerMatch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConformationExplorerMatch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConformationExplorer_Oval]"));
    core::BuiltInClass_sp classchem__ConformationExplorer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConformationExplorer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConformationExplorer_Oval,_lisp,chem::ConformationExplorer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConformationExplorer_O>::id,chem::ConformationExplorer_O::static_classSymbol());
    chem::ConformationExplorer_O::___staticClass = classchem__ConformationExplorer_Oval;
#ifdef USE_MPS
    chem::ConformationExplorer_O::static_Kind = gctools::GCKind<chem::ConformationExplorer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConformationExplorer_Oval,chem::ConformationExplorer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConformationExplorer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConformationExplorer_O>>::allocateClass();
        chem::ConformationExplorer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConformationExplorer_O::static_className() % (void*)(chem::ConformationExplorer_O::static_allocator) );
    classchem__ConformationExplorer_Oval->setCreator(chem::ConformationExplorer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConformationExplorer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConformationExplorer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConstitutionAtom_Oval]"));
    core::BuiltInClass_sp classchem__ConstitutionAtom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConstitutionAtom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConstitutionAtom_Oval,_lisp,chem::ConstitutionAtom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConstitutionAtom_O>::id,chem::ConstitutionAtom_O::static_classSymbol());
    chem::ConstitutionAtom_O::___staticClass = classchem__ConstitutionAtom_Oval;
#ifdef USE_MPS
    chem::ConstitutionAtom_O::static_Kind = gctools::GCKind<chem::ConstitutionAtom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConstitutionAtom_Oval,chem::ConstitutionAtom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConstitutionAtom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConstitutionAtom_O>>::allocateClass();
        chem::ConstitutionAtom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConstitutionAtom_O::static_className() % (void*)(chem::ConstitutionAtom_O::static_allocator) );
    classchem__ConstitutionAtom_Oval->setCreator(chem::ConstitutionAtom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConstitutionAtom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConstitutionAtom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConstitutionAtoms_Oval]"));
    core::BuiltInClass_sp classchem__ConstitutionAtoms_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConstitutionAtoms_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConstitutionAtoms_Oval,_lisp,chem::ConstitutionAtoms_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConstitutionAtoms_O>::id,chem::ConstitutionAtoms_O::static_classSymbol());
    chem::ConstitutionAtoms_O::___staticClass = classchem__ConstitutionAtoms_Oval;
#ifdef USE_MPS
    chem::ConstitutionAtoms_O::static_Kind = gctools::GCKind<chem::ConstitutionAtoms_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConstitutionAtoms_Oval,chem::ConstitutionAtoms_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConstitutionAtoms_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConstitutionAtoms_O>>::allocateClass();
        chem::ConstitutionAtoms_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConstitutionAtoms_O::static_className() % (void*)(chem::ConstitutionAtoms_O::static_allocator) );
    classchem__ConstitutionAtoms_Oval->setCreator(chem::ConstitutionAtoms_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConstitutionAtoms_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConstitutionAtoms_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConstitutionBond_Oval]"));
    core::BuiltInClass_sp classchem__ConstitutionBond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConstitutionBond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConstitutionBond_Oval,_lisp,chem::ConstitutionBond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConstitutionBond_O>::id,chem::ConstitutionBond_O::static_classSymbol());
    chem::ConstitutionBond_O::___staticClass = classchem__ConstitutionBond_Oval;
#ifdef USE_MPS
    chem::ConstitutionBond_O::static_Kind = gctools::GCKind<chem::ConstitutionBond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConstitutionBond_Oval,chem::ConstitutionBond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConstitutionBond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConstitutionBond_O>>::allocateClass();
        chem::ConstitutionBond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConstitutionBond_O::static_className() % (void*)(chem::ConstitutionBond_O::static_allocator) );
    classchem__ConstitutionBond_Oval->setCreator(chem::ConstitutionBond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConstitutionBond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConstitutionBond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CoordinateSystem_Oval]"));
    core::BuiltInClass_sp classchem__CoordinateSystem_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CoordinateSystem_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CoordinateSystem_Oval,_lisp,chem::CoordinateSystem_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CoordinateSystem_O>::id,chem::CoordinateSystem_O::static_classSymbol());
    chem::CoordinateSystem_O::___staticClass = classchem__CoordinateSystem_Oval;
#ifdef USE_MPS
    chem::CoordinateSystem_O::static_Kind = gctools::GCKind<chem::CoordinateSystem_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CoordinateSystem_Oval,chem::CoordinateSystem_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CoordinateSystem_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CoordinateSystem_O>>::allocateClass();
        chem::CoordinateSystem_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CoordinateSystem_O::static_className() % (void*)(chem::CoordinateSystem_O::static_allocator) );
    classchem__CoordinateSystem_Oval->setCreator(chem::CoordinateSystem_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CoordinateSystem_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CoordinateSystem_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Coupling_Oval]"));
    core::BuiltInClass_sp classchem__Coupling_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Coupling_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Coupling_Oval,_lisp,chem::Coupling_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Coupling_O>::id,chem::Coupling_O::static_classSymbol());
    chem::Coupling_O::___staticClass = classchem__Coupling_Oval;
#ifdef USE_MPS
    chem::Coupling_O::static_Kind = gctools::GCKind<chem::Coupling_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Coupling_Oval,chem::Coupling_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Coupling_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Coupling_O>>::allocateClass();
        chem::Coupling_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Coupling_O::static_className() % (void*)(chem::Coupling_O::static_allocator) );
    classchem__Coupling_Oval->setCreator(chem::Coupling_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Coupling_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Coupling_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyComponent_Oval]"));
    core::BuiltInClass_sp classchem__EnergyComponent_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyComponent_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyComponent_Oval,_lisp,chem::EnergyComponent_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyComponent_O>::id,chem::EnergyComponent_O::static_classSymbol());
    chem::EnergyComponent_O::___staticClass = classchem__EnergyComponent_Oval;
#ifdef USE_MPS
    chem::EnergyComponent_O::static_Kind = gctools::GCKind<chem::EnergyComponent_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyComponent_Oval,chem::EnergyComponent_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyComponent_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyComponent_O>>::allocateClass();
        chem::EnergyComponent_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyComponent_O::static_className() % (void*)(chem::EnergyComponent_O::static_allocator) );
    classchem__EnergyComponent_Oval->setCreator(chem::EnergyComponent_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyComponent_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyComponent_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyFunction_Oval]"));
    core::BuiltInClass_sp classchem__EnergyFunction_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyFunction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyFunction_Oval,_lisp,chem::EnergyFunction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyFunction_O>::id,chem::EnergyFunction_O::static_classSymbol());
    chem::EnergyFunction_O::___staticClass = classchem__EnergyFunction_Oval;
#ifdef USE_MPS
    chem::EnergyFunction_O::static_Kind = gctools::GCKind<chem::EnergyFunction_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyFunction_Oval,chem::EnergyFunction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyFunction_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyFunction_O>>::allocateClass();
        chem::EnergyFunction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyFunction_O::static_className() % (void*)(chem::EnergyFunction_O::static_allocator) );
    classchem__EnergyFunction_Oval->setCreator(chem::EnergyFunction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyFunction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyFunction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Entity_Oval]"));
    core::BuiltInClass_sp classchem__Entity_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Entity_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Entity_Oval,_lisp,chem::Entity_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Entity_O>::id,chem::Entity_O::static_classSymbol());
    chem::Entity_O::___staticClass = classchem__Entity_Oval;
#ifdef USE_MPS
    chem::Entity_O::static_Kind = gctools::GCKind<chem::Entity_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Entity_Oval,chem::Entity_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Entity_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Entity_O>>::allocateClass();
        chem::Entity_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Entity_O::static_className() % (void*)(chem::Entity_O::static_allocator) );
    classchem__Entity_Oval->setCreator(chem::Entity_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Entity_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Entity_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFBaseDb_Oval]"));
    core::BuiltInClass_sp classchem__FFBaseDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFBaseDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFBaseDb_Oval,_lisp,chem::FFBaseDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFBaseDb_O>::id,chem::FFBaseDb_O::static_classSymbol());
    chem::FFBaseDb_O::___staticClass = classchem__FFBaseDb_Oval;
#ifdef USE_MPS
    chem::FFBaseDb_O::static_Kind = gctools::GCKind<chem::FFBaseDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFBaseDb_Oval,chem::FFBaseDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFBaseDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFBaseDb_O>>::allocateClass();
        chem::FFBaseDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFBaseDb_O::static_className() % (void*)(chem::FFBaseDb_O::static_allocator) );
    classchem__FFBaseDb_Oval->setCreator(chem::FFBaseDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFBaseDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFBaseDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFNonbondCrossTermTable_Oval]"));
    core::BuiltInClass_sp classchem__FFNonbondCrossTermTable_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFNonbondCrossTermTable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFNonbondCrossTermTable_Oval,_lisp,chem::FFNonbondCrossTermTable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFNonbondCrossTermTable_O>::id,chem::FFNonbondCrossTermTable_O::static_classSymbol());
    chem::FFNonbondCrossTermTable_O::___staticClass = classchem__FFNonbondCrossTermTable_Oval;
#ifdef USE_MPS
    chem::FFNonbondCrossTermTable_O::static_Kind = gctools::GCKind<chem::FFNonbondCrossTermTable_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFNonbondCrossTermTable_Oval,chem::FFNonbondCrossTermTable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFNonbondCrossTermTable_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFNonbondCrossTermTable_O>>::allocateClass();
        chem::FFNonbondCrossTermTable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFNonbondCrossTermTable_O::static_className() % (void*)(chem::FFNonbondCrossTermTable_O::static_allocator) );
    classchem__FFNonbondCrossTermTable_Oval->setCreator(chem::FFNonbondCrossTermTable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFNonbondCrossTermTable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFNonbondCrossTermTable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFParameter_Oval]"));
    core::BuiltInClass_sp classchem__FFParameter_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFParameter_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFParameter_Oval,_lisp,chem::FFParameter_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFParameter_O>::id,chem::FFParameter_O::static_classSymbol());
    chem::FFParameter_O::___staticClass = classchem__FFParameter_Oval;
#ifdef USE_MPS
    chem::FFParameter_O::static_Kind = gctools::GCKind<chem::FFParameter_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFParameter_Oval,chem::FFParameter_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFParameter_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFParameter_O>>::allocateClass();
        chem::FFParameter_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFParameter_O::static_className() % (void*)(chem::FFParameter_O::static_allocator) );
    classchem__FFParameter_Oval->setCreator(chem::FFParameter_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFParameter_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFParameter_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ForceField_Oval]"));
    core::BuiltInClass_sp classchem__ForceField_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ForceField_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ForceField_Oval,_lisp,chem::ForceField_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ForceField_O>::id,chem::ForceField_O::static_classSymbol());
    chem::ForceField_O::___staticClass = classchem__ForceField_Oval;
#ifdef USE_MPS
    chem::ForceField_O::static_Kind = gctools::GCKind<chem::ForceField_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ForceField_Oval,chem::ForceField_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ForceField_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ForceField_O>>::allocateClass();
        chem::ForceField_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ForceField_O::static_className() % (void*)(chem::ForceField_O::static_allocator) );
    classchem__ForceField_Oval->setCreator(chem::ForceField_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ForceField_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ForceField_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FrameRecognizer_Oval]"));
    core::BuiltInClass_sp classchem__FrameRecognizer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FrameRecognizer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FrameRecognizer_Oval,_lisp,chem::FrameRecognizer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FrameRecognizer_O>::id,chem::FrameRecognizer_O::static_classSymbol());
    chem::FrameRecognizer_O::___staticClass = classchem__FrameRecognizer_Oval;
#ifdef USE_MPS
    chem::FrameRecognizer_O::static_Kind = gctools::GCKind<chem::FrameRecognizer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FrameRecognizer_Oval,chem::FrameRecognizer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FrameRecognizer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FrameRecognizer_O>>::allocateClass();
        chem::FrameRecognizer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FrameRecognizer_O::static_className() % (void*)(chem::FrameRecognizer_O::static_allocator) );
    classchem__FrameRecognizer_Oval->setCreator(chem::FrameRecognizer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FrameRecognizer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FrameRecognizer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ImproperTorsion_Oval]"));
    core::BuiltInClass_sp classchem__ImproperTorsion_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ImproperTorsion_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ImproperTorsion_Oval,_lisp,chem::ImproperTorsion_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ImproperTorsion_O>::id,chem::ImproperTorsion_O::static_classSymbol());
    chem::ImproperTorsion_O::___staticClass = classchem__ImproperTorsion_Oval;
#ifdef USE_MPS
    chem::ImproperTorsion_O::static_Kind = gctools::GCKind<chem::ImproperTorsion_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ImproperTorsion_Oval,chem::ImproperTorsion_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ImproperTorsion_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ImproperTorsion_O>>::allocateClass();
        chem::ImproperTorsion_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ImproperTorsion_O::static_className() % (void*)(chem::ImproperTorsion_O::static_allocator) );
    classchem__ImproperTorsion_Oval->setCreator(chem::ImproperTorsion_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ImproperTorsion_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ImproperTorsion_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__InfoDb_Oval]"));
    core::BuiltInClass_sp classchem__InfoDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__InfoDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__InfoDb_Oval,_lisp,chem::InfoDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::InfoDb_O>::id,chem::InfoDb_O::static_classSymbol());
    chem::InfoDb_O::___staticClass = classchem__InfoDb_Oval;
#ifdef USE_MPS
    chem::InfoDb_O::static_Kind = gctools::GCKind<chem::InfoDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__InfoDb_Oval,chem::InfoDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::InfoDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::InfoDb_O>>::allocateClass();
        chem::InfoDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::InfoDb_O::static_className() % (void*)(chem::InfoDb_O::static_allocator) );
    classchem__InfoDb_Oval->setCreator(chem::InfoDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::InfoDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__InfoDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MacroModelFile_Oval]"));
    core::BuiltInClass_sp classchem__MacroModelFile_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MacroModelFile_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MacroModelFile_Oval,_lisp,chem::MacroModelFile_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MacroModelFile_O>::id,chem::MacroModelFile_O::static_classSymbol());
    chem::MacroModelFile_O::___staticClass = classchem__MacroModelFile_Oval;
#ifdef USE_MPS
    chem::MacroModelFile_O::static_Kind = gctools::GCKind<chem::MacroModelFile_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MacroModelFile_Oval,chem::MacroModelFile_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MacroModelFile_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MacroModelFile_O>>::allocateClass();
        chem::MacroModelFile_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MacroModelFile_O::static_className() % (void*)(chem::MacroModelFile_O::static_allocator) );
    classchem__MacroModelFile_Oval->setCreator(chem::MacroModelFile_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MacroModelFile_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MacroModelFile_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MapOfMonomerNamesToAtomIndexers_Oval]"));
    core::BuiltInClass_sp classchem__MapOfMonomerNamesToAtomIndexers_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MapOfMonomerNamesToAtomIndexers_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MapOfMonomerNamesToAtomIndexers_Oval,_lisp,chem::MapOfMonomerNamesToAtomIndexers_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MapOfMonomerNamesToAtomIndexers_O>::id,chem::MapOfMonomerNamesToAtomIndexers_O::static_classSymbol());
    chem::MapOfMonomerNamesToAtomIndexers_O::___staticClass = classchem__MapOfMonomerNamesToAtomIndexers_Oval;
#ifdef USE_MPS
    chem::MapOfMonomerNamesToAtomIndexers_O::static_Kind = gctools::GCKind<chem::MapOfMonomerNamesToAtomIndexers_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MapOfMonomerNamesToAtomIndexers_Oval,chem::MapOfMonomerNamesToAtomIndexers_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MapOfMonomerNamesToAtomIndexers_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MapOfMonomerNamesToAtomIndexers_O>>::allocateClass();
        chem::MapOfMonomerNamesToAtomIndexers_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MapOfMonomerNamesToAtomIndexers_O::static_className() % (void*)(chem::MapOfMonomerNamesToAtomIndexers_O::static_allocator) );
    classchem__MapOfMonomerNamesToAtomIndexers_Oval->setCreator(chem::MapOfMonomerNamesToAtomIndexers_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MapOfMonomerNamesToAtomIndexers_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MapOfMonomerNamesToAtomIndexers_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Matter_Oval]"));
    core::BuiltInClass_sp classchem__Matter_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Matter_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Matter_Oval,_lisp,chem::Matter_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Matter_O>::id,chem::Matter_O::static_classSymbol());
    chem::Matter_O::___staticClass = classchem__Matter_Oval;
#ifdef USE_MPS
    chem::Matter_O::static_Kind = gctools::GCKind<chem::Matter_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Matter_Oval,chem::Matter_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Matter_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Matter_O>>::allocateClass();
        chem::Matter_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Matter_O::static_className() % (void*)(chem::Matter_O::static_allocator) );
    classchem__Matter_Oval->setCreator(chem::Matter_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Matter_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Matter_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MinimizerLog_Oval]"));
    core::BuiltInClass_sp classchem__MinimizerLog_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MinimizerLog_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MinimizerLog_Oval,_lisp,chem::MinimizerLog_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MinimizerLog_O>::id,chem::MinimizerLog_O::static_classSymbol());
    chem::MinimizerLog_O::___staticClass = classchem__MinimizerLog_Oval;
#ifdef USE_MPS
    chem::MinimizerLog_O::static_Kind = gctools::GCKind<chem::MinimizerLog_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MinimizerLog_Oval,chem::MinimizerLog_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MinimizerLog_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MinimizerLog_O>>::allocateClass();
        chem::MinimizerLog_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MinimizerLog_O::static_className() % (void*)(chem::MinimizerLog_O::static_allocator) );
    classchem__MinimizerLog_Oval->setCreator(chem::MinimizerLog_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MinimizerLog_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MinimizerLog_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Minimizer_Oval]"));
    core::BuiltInClass_sp classchem__Minimizer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Minimizer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Minimizer_Oval,_lisp,chem::Minimizer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Minimizer_O>::id,chem::Minimizer_O::static_classSymbol());
    chem::Minimizer_O::___staticClass = classchem__Minimizer_Oval;
#ifdef USE_MPS
    chem::Minimizer_O::static_Kind = gctools::GCKind<chem::Minimizer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Minimizer_Oval,chem::Minimizer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Minimizer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Minimizer_O>>::allocateClass();
        chem::Minimizer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Minimizer_O::static_className() % (void*)(chem::Minimizer_O::static_allocator) );
    classchem__Minimizer_Oval->setCreator(chem::Minimizer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Minimizer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Minimizer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MonomerContext_Oval]"));
    core::BuiltInClass_sp classchem__MonomerContext_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MonomerContext_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MonomerContext_Oval,_lisp,chem::MonomerContext_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MonomerContext_O>::id,chem::MonomerContext_O::static_classSymbol());
    chem::MonomerContext_O::___staticClass = classchem__MonomerContext_Oval;
#ifdef USE_MPS
    chem::MonomerContext_O::static_Kind = gctools::GCKind<chem::MonomerContext_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MonomerContext_Oval,chem::MonomerContext_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MonomerContext_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MonomerContext_O>>::allocateClass();
        chem::MonomerContext_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MonomerContext_O::static_className() % (void*)(chem::MonomerContext_O::static_allocator) );
    classchem__MonomerContext_Oval->setCreator(chem::MonomerContext_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MonomerContext_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MonomerContext_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MonomerCoordinates_Oval]"));
    core::BuiltInClass_sp classchem__MonomerCoordinates_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MonomerCoordinates_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MonomerCoordinates_Oval,_lisp,chem::MonomerCoordinates_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MonomerCoordinates_O>::id,chem::MonomerCoordinates_O::static_classSymbol());
    chem::MonomerCoordinates_O::___staticClass = classchem__MonomerCoordinates_Oval;
#ifdef USE_MPS
    chem::MonomerCoordinates_O::static_Kind = gctools::GCKind<chem::MonomerCoordinates_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MonomerCoordinates_Oval,chem::MonomerCoordinates_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MonomerCoordinates_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MonomerCoordinates_O>>::allocateClass();
        chem::MonomerCoordinates_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MonomerCoordinates_O::static_className() % (void*)(chem::MonomerCoordinates_O::static_allocator) );
    classchem__MonomerCoordinates_Oval->setCreator(chem::MonomerCoordinates_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MonomerCoordinates_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MonomerCoordinates_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__NVector_Oval]"));
    core::BuiltInClass_sp classchem__NVector_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__NVector_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__NVector_Oval,_lisp,chem::NVector_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::NVector_O>::id,chem::NVector_O::static_classSymbol());
    chem::NVector_O::___staticClass = classchem__NVector_Oval;
#ifdef USE_MPS
    chem::NVector_O::static_Kind = gctools::GCKind<chem::NVector_O>::Kind;
#endif
    core::af_setf_findClass(classchem__NVector_Oval,chem::NVector_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::NVector_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::NVector_O>>::allocateClass();
        chem::NVector_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::NVector_O::static_className() % (void*)(chem::NVector_O::static_allocator) );
    classchem__NVector_Oval->setCreator(chem::NVector_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::NVector_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__NVector_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__NumericalFunction_Oval]"));
    core::BuiltInClass_sp classchem__NumericalFunction_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__NumericalFunction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__NumericalFunction_Oval,_lisp,chem::NumericalFunction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::NumericalFunction_O>::id,chem::NumericalFunction_O::static_classSymbol());
    chem::NumericalFunction_O::___staticClass = classchem__NumericalFunction_Oval;
#ifdef USE_MPS
    chem::NumericalFunction_O::static_Kind = gctools::GCKind<chem::NumericalFunction_O>::Kind;
#endif
    core::af_setf_findClass(classchem__NumericalFunction_Oval,chem::NumericalFunction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::NumericalFunction_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::NumericalFunction_O>>::allocateClass();
        chem::NumericalFunction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::NumericalFunction_O::static_className() % (void*)(chem::NumericalFunction_O::static_allocator) );
    classchem__NumericalFunction_Oval->setCreator(chem::NumericalFunction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::NumericalFunction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__NumericalFunction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__OligomerPart_Base_Oval]"));
    core::BuiltInClass_sp classchem__OligomerPart_Base_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__OligomerPart_Base_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__OligomerPart_Base_Oval,_lisp,chem::OligomerPart_Base_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::OligomerPart_Base_O>::id,chem::OligomerPart_Base_O::static_classSymbol());
    chem::OligomerPart_Base_O::___staticClass = classchem__OligomerPart_Base_Oval;
#ifdef USE_MPS
    chem::OligomerPart_Base_O::static_Kind = gctools::GCKind<chem::OligomerPart_Base_O>::Kind;
#endif
    core::af_setf_findClass(classchem__OligomerPart_Base_Oval,chem::OligomerPart_Base_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::OligomerPart_Base_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::OligomerPart_Base_O>>::allocateClass();
        chem::OligomerPart_Base_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::OligomerPart_Base_O::static_className() % (void*)(chem::OligomerPart_Base_O::static_allocator) );
    classchem__OligomerPart_Base_Oval->setCreator(chem::OligomerPart_Base_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::OligomerPart_Base_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__OligomerPart_Base_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Oligomer_Oval]"));
    core::BuiltInClass_sp classchem__Oligomer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Oligomer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Oligomer_Oval,_lisp,chem::Oligomer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Oligomer_O>::id,chem::Oligomer_O::static_classSymbol());
    chem::Oligomer_O::___staticClass = classchem__Oligomer_Oval;
#ifdef USE_MPS
    chem::Oligomer_O::static_Kind = gctools::GCKind<chem::Oligomer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Oligomer_Oval,chem::Oligomer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Oligomer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Oligomer_O>>::allocateClass();
        chem::Oligomer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Oligomer_O::static_className() % (void*)(chem::Oligomer_O::static_allocator) );
    classchem__Oligomer_Oval->setCreator(chem::Oligomer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Oligomer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Oligomer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__OneMonomer_Oval]"));
    core::BuiltInClass_sp classchem__OneMonomer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__OneMonomer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__OneMonomer_Oval,_lisp,chem::OneMonomer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::OneMonomer_O>::id,chem::OneMonomer_O::static_classSymbol());
    chem::OneMonomer_O::___staticClass = classchem__OneMonomer_Oval;
#ifdef USE_MPS
    chem::OneMonomer_O::static_Kind = gctools::GCKind<chem::OneMonomer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__OneMonomer_Oval,chem::OneMonomer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::OneMonomer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::OneMonomer_O>>::allocateClass();
        chem::OneMonomer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::OneMonomer_O::static_className() % (void*)(chem::OneMonomer_O::static_allocator) );
    classchem__OneMonomer_Oval->setCreator(chem::OneMonomer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::OneMonomer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__OneMonomer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__OneTypeRule_Oval]"));
    core::BuiltInClass_sp classchem__OneTypeRule_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__OneTypeRule_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__OneTypeRule_Oval,_lisp,chem::OneTypeRule_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::OneTypeRule_O>::id,chem::OneTypeRule_O::static_classSymbol());
    chem::OneTypeRule_O::___staticClass = classchem__OneTypeRule_Oval;
#ifdef USE_MPS
    chem::OneTypeRule_O::static_Kind = gctools::GCKind<chem::OneTypeRule_O>::Kind;
#endif
    core::af_setf_findClass(classchem__OneTypeRule_Oval,chem::OneTypeRule_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::OneTypeRule_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::OneTypeRule_O>>::allocateClass();
        chem::OneTypeRule_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::OneTypeRule_O::static_className() % (void*)(chem::OneTypeRule_O::static_allocator) );
    classchem__OneTypeRule_Oval->setCreator(chem::OneTypeRule_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::OneTypeRule_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__OneTypeRule_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__PathMessage_Oval]"));
    core::BuiltInClass_sp classchem__PathMessage_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__PathMessage_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__PathMessage_Oval,_lisp,chem::PathMessage_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::PathMessage_O>::id,chem::PathMessage_O::static_classSymbol());
    chem::PathMessage_O::___staticClass = classchem__PathMessage_Oval;
#ifdef USE_MPS
    chem::PathMessage_O::static_Kind = gctools::GCKind<chem::PathMessage_O>::Kind;
#endif
    core::af_setf_findClass(classchem__PathMessage_Oval,chem::PathMessage_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::PathMessage_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::PathMessage_O>>::allocateClass();
        chem::PathMessage_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::PathMessage_O::static_className() % (void*)(chem::PathMessage_O::static_allocator) );
    classchem__PathMessage_Oval->setCreator(chem::PathMessage_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::PathMessage_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__PathMessage_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__PdbMonomerConnectivity_Oval]"));
    core::BuiltInClass_sp classchem__PdbMonomerConnectivity_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__PdbMonomerConnectivity_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__PdbMonomerConnectivity_Oval,_lisp,chem::PdbMonomerConnectivity_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::PdbMonomerConnectivity_O>::id,chem::PdbMonomerConnectivity_O::static_classSymbol());
    chem::PdbMonomerConnectivity_O::___staticClass = classchem__PdbMonomerConnectivity_Oval;
#ifdef USE_MPS
    chem::PdbMonomerConnectivity_O::static_Kind = gctools::GCKind<chem::PdbMonomerConnectivity_O>::Kind;
#endif
    core::af_setf_findClass(classchem__PdbMonomerConnectivity_Oval,chem::PdbMonomerConnectivity_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::PdbMonomerConnectivity_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::PdbMonomerConnectivity_O>>::allocateClass();
        chem::PdbMonomerConnectivity_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::PdbMonomerConnectivity_O::static_className() % (void*)(chem::PdbMonomerConnectivity_O::static_allocator) );
    classchem__PdbMonomerConnectivity_Oval->setCreator(chem::PdbMonomerConnectivity_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::PdbMonomerConnectivity_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__PdbMonomerConnectivity_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__PdbMonomerDatabase_Oval]"));
    core::BuiltInClass_sp classchem__PdbMonomerDatabase_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__PdbMonomerDatabase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__PdbMonomerDatabase_Oval,_lisp,chem::PdbMonomerDatabase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::PdbMonomerDatabase_O>::id,chem::PdbMonomerDatabase_O::static_classSymbol());
    chem::PdbMonomerDatabase_O::___staticClass = classchem__PdbMonomerDatabase_Oval;
#ifdef USE_MPS
    chem::PdbMonomerDatabase_O::static_Kind = gctools::GCKind<chem::PdbMonomerDatabase_O>::Kind;
#endif
    core::af_setf_findClass(classchem__PdbMonomerDatabase_Oval,chem::PdbMonomerDatabase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::PdbMonomerDatabase_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::PdbMonomerDatabase_O>>::allocateClass();
        chem::PdbMonomerDatabase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::PdbMonomerDatabase_O::static_className() % (void*)(chem::PdbMonomerDatabase_O::static_allocator) );
    classchem__PdbMonomerDatabase_Oval->setCreator(chem::PdbMonomerDatabase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::PdbMonomerDatabase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__PdbMonomerDatabase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__PdbReader_Oval]"));
    core::BuiltInClass_sp classchem__PdbReader_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__PdbReader_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__PdbReader_Oval,_lisp,chem::PdbReader_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::PdbReader_O>::id,chem::PdbReader_O::static_classSymbol());
    chem::PdbReader_O::___staticClass = classchem__PdbReader_Oval;
#ifdef USE_MPS
    chem::PdbReader_O::static_Kind = gctools::GCKind<chem::PdbReader_O>::Kind;
#endif
    core::af_setf_findClass(classchem__PdbReader_Oval,chem::PdbReader_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::PdbReader_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::PdbReader_O>>::allocateClass();
        chem::PdbReader_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::PdbReader_O::static_className() % (void*)(chem::PdbReader_O::static_allocator) );
    classchem__PdbReader_Oval->setCreator(chem::PdbReader_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::PdbReader_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__PdbReader_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__PdbWriter_Oval]"));
    core::BuiltInClass_sp classchem__PdbWriter_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__PdbWriter_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__PdbWriter_Oval,_lisp,chem::PdbWriter_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::PdbWriter_O>::id,chem::PdbWriter_O::static_classSymbol());
    chem::PdbWriter_O::___staticClass = classchem__PdbWriter_Oval;
#ifdef USE_MPS
    chem::PdbWriter_O::static_Kind = gctools::GCKind<chem::PdbWriter_O>::Kind;
#endif
    core::af_setf_findClass(classchem__PdbWriter_Oval,chem::PdbWriter_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::PdbWriter_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::PdbWriter_O>>::allocateClass();
        chem::PdbWriter_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::PdbWriter_O::static_className() % (void*)(chem::PdbWriter_O::static_allocator) );
    classchem__PdbWriter_Oval->setCreator(chem::PdbWriter_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::PdbWriter_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__PdbWriter_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Plug_Oval]"));
    core::BuiltInClass_sp classchem__Plug_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Plug_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Plug_Oval,_lisp,chem::Plug_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Plug_O>::id,chem::Plug_O::static_classSymbol());
    chem::Plug_O::___staticClass = classchem__Plug_Oval;
#ifdef USE_MPS
    chem::Plug_O::static_Kind = gctools::GCKind<chem::Plug_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Plug_Oval,chem::Plug_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Plug_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Plug_O>>::allocateClass();
        chem::Plug_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Plug_O::static_className() % (void*)(chem::Plug_O::static_allocator) );
    classchem__Plug_Oval->setCreator(chem::Plug_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Plug_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Plug_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ProperTorsion_Oval]"));
    core::BuiltInClass_sp classchem__ProperTorsion_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ProperTorsion_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ProperTorsion_Oval,_lisp,chem::ProperTorsion_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ProperTorsion_O>::id,chem::ProperTorsion_O::static_classSymbol());
    chem::ProperTorsion_O::___staticClass = classchem__ProperTorsion_Oval;
#ifdef USE_MPS
    chem::ProperTorsion_O::static_Kind = gctools::GCKind<chem::ProperTorsion_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ProperTorsion_Oval,chem::ProperTorsion_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ProperTorsion_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ProperTorsion_O>>::allocateClass();
        chem::ProperTorsion_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ProperTorsion_O::static_className() % (void*)(chem::ProperTorsion_O::static_allocator) );
    classchem__ProperTorsion_Oval->setCreator(chem::ProperTorsion_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ProperTorsion_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ProperTorsion_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ReadAmberParameters_Oval]"));
    core::BuiltInClass_sp classchem__ReadAmberParameters_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ReadAmberParameters_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ReadAmberParameters_Oval,_lisp,chem::ReadAmberParameters_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ReadAmberParameters_O>::id,chem::ReadAmberParameters_O::static_classSymbol());
    chem::ReadAmberParameters_O::___staticClass = classchem__ReadAmberParameters_Oval;
#ifdef USE_MPS
    chem::ReadAmberParameters_O::static_Kind = gctools::GCKind<chem::ReadAmberParameters_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ReadAmberParameters_Oval,chem::ReadAmberParameters_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ReadAmberParameters_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ReadAmberParameters_O>>::allocateClass();
        chem::ReadAmberParameters_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ReadAmberParameters_O::static_className() % (void*)(chem::ReadAmberParameters_O::static_allocator) );
    classchem__ReadAmberParameters_Oval->setCreator(chem::ReadAmberParameters_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ReadAmberParameters_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ReadAmberParameters_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ReportBase_Oval]"));
    core::BuiltInClass_sp classchem__ReportBase_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ReportBase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ReportBase_Oval,_lisp,chem::ReportBase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ReportBase_O>::id,chem::ReportBase_O::static_classSymbol());
    chem::ReportBase_O::___staticClass = classchem__ReportBase_Oval;
#ifdef USE_MPS
    chem::ReportBase_O::static_Kind = gctools::GCKind<chem::ReportBase_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ReportBase_Oval,chem::ReportBase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ReportBase_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ReportBase_O>>::allocateClass();
        chem::ReportBase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ReportBase_O::static_className() % (void*)(chem::ReportBase_O::static_allocator) );
    classchem__ReportBase_Oval->setCreator(chem::ReportBase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ReportBase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ReportBase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ResidueList_Oval]"));
    core::BuiltInClass_sp classchem__ResidueList_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ResidueList_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ResidueList_Oval,_lisp,chem::ResidueList_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ResidueList_O>::id,chem::ResidueList_O::static_classSymbol());
    chem::ResidueList_O::___staticClass = classchem__ResidueList_Oval;
#ifdef USE_MPS
    chem::ResidueList_O::static_Kind = gctools::GCKind<chem::ResidueList_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ResidueList_Oval,chem::ResidueList_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ResidueList_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ResidueList_O>>::allocateClass();
        chem::ResidueList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ResidueList_O::static_className() % (void*)(chem::ResidueList_O::static_allocator) );
    classchem__ResidueList_Oval->setCreator(chem::ResidueList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ResidueList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ResidueList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintList_Oval]"));
    core::BuiltInClass_sp classchem__RestraintList_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintList_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintList_Oval,_lisp,chem::RestraintList_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintList_O>::id,chem::RestraintList_O::static_classSymbol());
    chem::RestraintList_O::___staticClass = classchem__RestraintList_Oval;
#ifdef USE_MPS
    chem::RestraintList_O::static_Kind = gctools::GCKind<chem::RestraintList_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintList_Oval,chem::RestraintList_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintList_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintList_O>>::allocateClass();
        chem::RestraintList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintList_O::static_className() % (void*)(chem::RestraintList_O::static_allocator) );
    classchem__RestraintList_Oval->setCreator(chem::RestraintList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Restraint_Oval]"));
    core::BuiltInClass_sp classchem__Restraint_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Restraint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Restraint_Oval,_lisp,chem::Restraint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Restraint_O>::id,chem::Restraint_O::static_classSymbol());
    chem::Restraint_O::___staticClass = classchem__Restraint_Oval;
#ifdef USE_MPS
    chem::Restraint_O::static_Kind = gctools::GCKind<chem::Restraint_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Restraint_Oval,chem::Restraint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Restraint_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Restraint_O>>::allocateClass();
        chem::Restraint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Restraint_O::static_className() % (void*)(chem::Restraint_O::static_allocator) );
    classchem__Restraint_Oval->setCreator(chem::Restraint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Restraint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Restraint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RingFinder_Oval]"));
    core::BuiltInClass_sp classchem__RingFinder_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RingFinder_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RingFinder_Oval,_lisp,chem::RingFinder_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RingFinder_O>::id,chem::RingFinder_O::static_classSymbol());
    chem::RingFinder_O::___staticClass = classchem__RingFinder_Oval;
#ifdef USE_MPS
    chem::RingFinder_O::static_Kind = gctools::GCKind<chem::RingFinder_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RingFinder_Oval,chem::RingFinder_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RingFinder_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RingFinder_O>>::allocateClass();
        chem::RingFinder_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RingFinder_O::static_className() % (void*)(chem::RingFinder_O::static_allocator) );
    classchem__RingFinder_Oval->setCreator(chem::RingFinder_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RingFinder_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RingFinder_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Rotamer_Oval]"));
    core::BuiltInClass_sp classchem__Rotamer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Rotamer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Rotamer_Oval,_lisp,chem::Rotamer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Rotamer_O>::id,chem::Rotamer_O::static_classSymbol());
    chem::Rotamer_O::___staticClass = classchem__Rotamer_Oval;
#ifdef USE_MPS
    chem::Rotamer_O::static_Kind = gctools::GCKind<chem::Rotamer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Rotamer_Oval,chem::Rotamer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Rotamer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Rotamer_O>>::allocateClass();
        chem::Rotamer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Rotamer_O::static_className() % (void*)(chem::Rotamer_O::static_allocator) );
    classchem__Rotamer_Oval->setCreator(chem::Rotamer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Rotamer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Rotamer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SearchStatistics_Oval]"));
    core::BuiltInClass_sp classchem__SearchStatistics_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SearchStatistics_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SearchStatistics_Oval,_lisp,chem::SearchStatistics_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SearchStatistics_O>::id,chem::SearchStatistics_O::static_classSymbol());
    chem::SearchStatistics_O::___staticClass = classchem__SearchStatistics_Oval;
#ifdef USE_MPS
    chem::SearchStatistics_O::static_Kind = gctools::GCKind<chem::SearchStatistics_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SearchStatistics_Oval,chem::SearchStatistics_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SearchStatistics_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SearchStatistics_O>>::allocateClass();
        chem::SearchStatistics_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SearchStatistics_O::static_className() % (void*)(chem::SearchStatistics_O::static_allocator) );
    classchem__SearchStatistics_Oval->setCreator(chem::SearchStatistics_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SearchStatistics_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SearchStatistics_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SpanningLoop_Oval]"));
    core::BuiltInClass_sp classchem__SpanningLoop_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SpanningLoop_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SpanningLoop_Oval,_lisp,chem::SpanningLoop_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SpanningLoop_O>::id,chem::SpanningLoop_O::static_classSymbol());
    chem::SpanningLoop_O::___staticClass = classchem__SpanningLoop_Oval;
#ifdef USE_MPS
    chem::SpanningLoop_O::static_Kind = gctools::GCKind<chem::SpanningLoop_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SpanningLoop_Oval,chem::SpanningLoop_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SpanningLoop_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SpanningLoop_O>>::allocateClass();
        chem::SpanningLoop_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SpanningLoop_O::static_className() % (void*)(chem::SpanningLoop_O::static_allocator) );
    classchem__SpanningLoop_Oval->setCreator(chem::SpanningLoop_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SpanningLoop_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SpanningLoop_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SpecificContextSet_Oval]"));
    core::BuiltInClass_sp classchem__SpecificContextSet_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SpecificContextSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SpecificContextSet_Oval,_lisp,chem::SpecificContextSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SpecificContextSet_O>::id,chem::SpecificContextSet_O::static_classSymbol());
    chem::SpecificContextSet_O::___staticClass = classchem__SpecificContextSet_Oval;
#ifdef USE_MPS
    chem::SpecificContextSet_O::static_Kind = gctools::GCKind<chem::SpecificContextSet_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SpecificContextSet_Oval,chem::SpecificContextSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SpecificContextSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SpecificContextSet_O>>::allocateClass();
        chem::SpecificContextSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SpecificContextSet_O::static_className() % (void*)(chem::SpecificContextSet_O::static_allocator) );
    classchem__SpecificContextSet_Oval->setCreator(chem::SpecificContextSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SpecificContextSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SpecificContextSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SpecificContext_Oval]"));
    core::BuiltInClass_sp classchem__SpecificContext_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SpecificContext_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SpecificContext_Oval,_lisp,chem::SpecificContext_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SpecificContext_O>::id,chem::SpecificContext_O::static_classSymbol());
    chem::SpecificContext_O::___staticClass = classchem__SpecificContext_Oval;
#ifdef USE_MPS
    chem::SpecificContext_O::static_Kind = gctools::GCKind<chem::SpecificContext_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SpecificContext_Oval,chem::SpecificContext_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SpecificContext_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SpecificContext_O>>::allocateClass();
        chem::SpecificContext_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SpecificContext_O::static_className() % (void*)(chem::SpecificContext_O::static_allocator) );
    classchem__SpecificContext_Oval->setCreator(chem::SpecificContext_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SpecificContext_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SpecificContext_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StatusTracker_Oval]"));
    core::BuiltInClass_sp classchem__StatusTracker_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StatusTracker_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StatusTracker_Oval,_lisp,chem::StatusTracker_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StatusTracker_O>::id,chem::StatusTracker_O::static_classSymbol());
    chem::StatusTracker_O::___staticClass = classchem__StatusTracker_Oval;
#ifdef USE_MPS
    chem::StatusTracker_O::static_Kind = gctools::GCKind<chem::StatusTracker_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StatusTracker_Oval,chem::StatusTracker_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StatusTracker_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StatusTracker_O>>::allocateClass();
        chem::StatusTracker_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StatusTracker_O::static_className() % (void*)(chem::StatusTracker_O::static_allocator) );
    classchem__StatusTracker_Oval->setCreator(chem::StatusTracker_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StatusTracker_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StatusTracker_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StereoConfiguration_Oval]"));
    core::BuiltInClass_sp classchem__StereoConfiguration_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StereoConfiguration_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StereoConfiguration_Oval,_lisp,chem::StereoConfiguration_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StereoConfiguration_O>::id,chem::StereoConfiguration_O::static_classSymbol());
    chem::StereoConfiguration_O::___staticClass = classchem__StereoConfiguration_Oval;
#ifdef USE_MPS
    chem::StereoConfiguration_O::static_Kind = gctools::GCKind<chem::StereoConfiguration_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StereoConfiguration_Oval,chem::StereoConfiguration_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StereoConfiguration_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StereoConfiguration_O>>::allocateClass();
        chem::StereoConfiguration_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StereoConfiguration_O::static_className() % (void*)(chem::StereoConfiguration_O::static_allocator) );
    classchem__StereoConfiguration_Oval->setCreator(chem::StereoConfiguration_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StereoConfiguration_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StereoConfiguration_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StereoInformation_Oval]"));
    core::BuiltInClass_sp classchem__StereoInformation_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StereoInformation_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StereoInformation_Oval,_lisp,chem::StereoInformation_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StereoInformation_O>::id,chem::StereoInformation_O::static_classSymbol());
    chem::StereoInformation_O::___staticClass = classchem__StereoInformation_Oval;
#ifdef USE_MPS
    chem::StereoInformation_O::static_Kind = gctools::GCKind<chem::StereoInformation_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StereoInformation_Oval,chem::StereoInformation_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StereoInformation_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StereoInformation_O>>::allocateClass();
        chem::StereoInformation_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StereoInformation_O::static_className() % (void*)(chem::StereoInformation_O::static_allocator) );
    classchem__StereoInformation_Oval->setCreator(chem::StereoInformation_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StereoInformation_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StereoInformation_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StereoisomerAtom_Oval]"));
    core::BuiltInClass_sp classchem__StereoisomerAtom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StereoisomerAtom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StereoisomerAtom_Oval,_lisp,chem::StereoisomerAtom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StereoisomerAtom_O>::id,chem::StereoisomerAtom_O::static_classSymbol());
    chem::StereoisomerAtom_O::___staticClass = classchem__StereoisomerAtom_Oval;
#ifdef USE_MPS
    chem::StereoisomerAtom_O::static_Kind = gctools::GCKind<chem::StereoisomerAtom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StereoisomerAtom_Oval,chem::StereoisomerAtom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StereoisomerAtom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StereoisomerAtom_O>>::allocateClass();
        chem::StereoisomerAtom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StereoisomerAtom_O::static_className() % (void*)(chem::StereoisomerAtom_O::static_allocator) );
    classchem__StereoisomerAtom_Oval->setCreator(chem::StereoisomerAtom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StereoisomerAtom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StereoisomerAtom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StereoisomerAtoms_Oval]"));
    core::BuiltInClass_sp classchem__StereoisomerAtoms_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StereoisomerAtoms_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StereoisomerAtoms_Oval,_lisp,chem::StereoisomerAtoms_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StereoisomerAtoms_O>::id,chem::StereoisomerAtoms_O::static_classSymbol());
    chem::StereoisomerAtoms_O::___staticClass = classchem__StereoisomerAtoms_Oval;
#ifdef USE_MPS
    chem::StereoisomerAtoms_O::static_Kind = gctools::GCKind<chem::StereoisomerAtoms_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StereoisomerAtoms_Oval,chem::StereoisomerAtoms_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StereoisomerAtoms_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StereoisomerAtoms_O>>::allocateClass();
        chem::StereoisomerAtoms_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StereoisomerAtoms_O::static_className() % (void*)(chem::StereoisomerAtoms_O::static_allocator) );
    classchem__StereoisomerAtoms_Oval->setCreator(chem::StereoisomerAtoms_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StereoisomerAtoms_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StereoisomerAtoms_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StructureComparer_Oval]"));
    core::BuiltInClass_sp classchem__StructureComparer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StructureComparer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StructureComparer_Oval,_lisp,chem::StructureComparer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StructureComparer_O>::id,chem::StructureComparer_O::static_classSymbol());
    chem::StructureComparer_O::___staticClass = classchem__StructureComparer_Oval;
#ifdef USE_MPS
    chem::StructureComparer_O::static_Kind = gctools::GCKind<chem::StructureComparer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StructureComparer_Oval,chem::StructureComparer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StructureComparer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StructureComparer_O>>::allocateClass();
        chem::StructureComparer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StructureComparer_O::static_className() % (void*)(chem::StructureComparer_O::static_allocator) );
    classchem__StructureComparer_Oval->setCreator(chem::StructureComparer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StructureComparer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StructureComparer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Structure_Old_ListEntry_Oval]"));
    core::BuiltInClass_sp classchem__Structure_Old_ListEntry_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Structure_Old_ListEntry_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Structure_Old_ListEntry_Oval,_lisp,chem::Structure_Old_ListEntry_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Structure_Old_ListEntry_O>::id,chem::Structure_Old_ListEntry_O::static_classSymbol());
    chem::Structure_Old_ListEntry_O::___staticClass = classchem__Structure_Old_ListEntry_Oval;
#ifdef USE_MPS
    chem::Structure_Old_ListEntry_O::static_Kind = gctools::GCKind<chem::Structure_Old_ListEntry_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Structure_Old_ListEntry_Oval,chem::Structure_Old_ListEntry_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Structure_Old_ListEntry_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Structure_Old_ListEntry_O>>::allocateClass();
        chem::Structure_Old_ListEntry_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Structure_Old_ListEntry_O::static_className() % (void*)(chem::Structure_Old_ListEntry_O::static_allocator) );
    classchem__Structure_Old_ListEntry_Oval->setCreator(chem::Structure_Old_ListEntry_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Structure_Old_ListEntry_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Structure_Old_ListEntry_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Structure_Old_List_Oval]"));
    core::BuiltInClass_sp classchem__Structure_Old_List_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Structure_Old_List_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Structure_Old_List_Oval,_lisp,chem::Structure_Old_List_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Structure_Old_List_O>::id,chem::Structure_Old_List_O::static_classSymbol());
    chem::Structure_Old_List_O::___staticClass = classchem__Structure_Old_List_Oval;
#ifdef USE_MPS
    chem::Structure_Old_List_O::static_Kind = gctools::GCKind<chem::Structure_Old_List_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Structure_Old_List_Oval,chem::Structure_Old_List_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Structure_Old_List_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Structure_Old_List_O>>::allocateClass();
        chem::Structure_Old_List_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Structure_Old_List_O::static_className() % (void*)(chem::Structure_Old_List_O::static_allocator) );
    classchem__Structure_Old_List_Oval->setCreator(chem::Structure_Old_List_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Structure_Old_List_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Structure_Old_List_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SuperposeEngine_Oval]"));
    core::BuiltInClass_sp classchem__SuperposeEngine_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SuperposeEngine_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SuperposeEngine_Oval,_lisp,chem::SuperposeEngine_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SuperposeEngine_O>::id,chem::SuperposeEngine_O::static_classSymbol());
    chem::SuperposeEngine_O::___staticClass = classchem__SuperposeEngine_Oval;
#ifdef USE_MPS
    chem::SuperposeEngine_O::static_Kind = gctools::GCKind<chem::SuperposeEngine_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SuperposeEngine_Oval,chem::SuperposeEngine_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SuperposeEngine_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SuperposeEngine_O>>::allocateClass();
        chem::SuperposeEngine_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SuperposeEngine_O::static_className() % (void*)(chem::SuperposeEngine_O::static_allocator) );
    classchem__SuperposeEngine_Oval->setCreator(chem::SuperposeEngine_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SuperposeEngine_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SuperposeEngine_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Topology_Oval]"));
    core::BuiltInClass_sp classchem__Topology_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Topology_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Topology_Oval,_lisp,chem::Topology_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Topology_O>::id,chem::Topology_O::static_classSymbol());
    chem::Topology_O::___staticClass = classchem__Topology_Oval;
#ifdef USE_MPS
    chem::Topology_O::static_Kind = gctools::GCKind<chem::Topology_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Topology_Oval,chem::Topology_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Topology_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Topology_O>>::allocateClass();
        chem::Topology_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Topology_O::static_className() % (void*)(chem::Topology_O::static_allocator) );
    classchem__Topology_Oval->setCreator(chem::Topology_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Topology_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Topology_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__TrajectoryFrame_Oval]"));
    core::BuiltInClass_sp classchem__TrajectoryFrame_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__TrajectoryFrame_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__TrajectoryFrame_Oval,_lisp,chem::TrajectoryFrame_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::TrajectoryFrame_O>::id,chem::TrajectoryFrame_O::static_classSymbol());
    chem::TrajectoryFrame_O::___staticClass = classchem__TrajectoryFrame_Oval;
#ifdef USE_MPS
    chem::TrajectoryFrame_O::static_Kind = gctools::GCKind<chem::TrajectoryFrame_O>::Kind;
#endif
    core::af_setf_findClass(classchem__TrajectoryFrame_Oval,chem::TrajectoryFrame_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::TrajectoryFrame_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::TrajectoryFrame_O>>::allocateClass();
        chem::TrajectoryFrame_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::TrajectoryFrame_O::static_className() % (void*)(chem::TrajectoryFrame_O::static_allocator) );
    classchem__TrajectoryFrame_Oval->setCreator(chem::TrajectoryFrame_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::TrajectoryFrame_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__TrajectoryFrame_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Trajectory_Oval]"));
    core::BuiltInClass_sp classchem__Trajectory_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Trajectory_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Trajectory_Oval,_lisp,chem::Trajectory_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Trajectory_O>::id,chem::Trajectory_O::static_classSymbol());
    chem::Trajectory_O::___staticClass = classchem__Trajectory_Oval;
#ifdef USE_MPS
    chem::Trajectory_O::static_Kind = gctools::GCKind<chem::Trajectory_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Trajectory_Oval,chem::Trajectory_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Trajectory_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Trajectory_O>>::allocateClass();
        chem::Trajectory_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Trajectory_O::static_className() % (void*)(chem::Trajectory_O::static_allocator) );
    classchem__Trajectory_Oval->setCreator(chem::Trajectory_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Trajectory_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Trajectory_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__TwisterDriver_Oval]"));
    core::BuiltInClass_sp classchem__TwisterDriver_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__TwisterDriver_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__TwisterDriver_Oval,_lisp,chem::TwisterDriver_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::TwisterDriver_O>::id,chem::TwisterDriver_O::static_classSymbol());
    chem::TwisterDriver_O::___staticClass = classchem__TwisterDriver_Oval;
#ifdef USE_MPS
    chem::TwisterDriver_O::static_Kind = gctools::GCKind<chem::TwisterDriver_O>::Kind;
#endif
    core::af_setf_findClass(classchem__TwisterDriver_Oval,chem::TwisterDriver_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::TwisterDriver_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::TwisterDriver_O>>::allocateClass();
        chem::TwisterDriver_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::TwisterDriver_O::static_className() % (void*)(chem::TwisterDriver_O::static_allocator) );
    classchem__TwisterDriver_Oval->setCreator(chem::TwisterDriver_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::TwisterDriver_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__TwisterDriver_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Twister_Oval]"));
    core::BuiltInClass_sp classchem__Twister_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Twister_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Twister_Oval,_lisp,chem::Twister_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Twister_O>::id,chem::Twister_O::static_classSymbol());
    chem::Twister_O::___staticClass = classchem__Twister_Oval;
#ifdef USE_MPS
    chem::Twister_O::static_Kind = gctools::GCKind<chem::Twister_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Twister_Oval,chem::Twister_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Twister_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Twister_O>>::allocateClass();
        chem::Twister_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Twister_O::static_className() % (void*)(chem::Twister_O::static_allocator) );
    classchem__Twister_Oval->setCreator(chem::Twister_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Twister_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Twister_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__TypeAssignmentRules_Oval]"));
    core::BuiltInClass_sp classchem__TypeAssignmentRules_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__TypeAssignmentRules_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__TypeAssignmentRules_Oval,_lisp,chem::TypeAssignmentRules_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::TypeAssignmentRules_O>::id,chem::TypeAssignmentRules_O::static_classSymbol());
    chem::TypeAssignmentRules_O::___staticClass = classchem__TypeAssignmentRules_Oval;
#ifdef USE_MPS
    chem::TypeAssignmentRules_O::static_Kind = gctools::GCKind<chem::TypeAssignmentRules_O>::Kind;
#endif
    core::af_setf_findClass(classchem__TypeAssignmentRules_Oval,chem::TypeAssignmentRules_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::TypeAssignmentRules_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::TypeAssignmentRules_O>>::allocateClass();
        chem::TypeAssignmentRules_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::TypeAssignmentRules_O::static_className() % (void*)(chem::TypeAssignmentRules_O::static_allocator) );
    classchem__TypeAssignmentRules_Oval->setCreator(chem::TypeAssignmentRules_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::TypeAssignmentRules_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__TypeAssignmentRules_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__VirtualSphere_Oval]"));
    core::BuiltInClass_sp classchem__VirtualSphere_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__VirtualSphere_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__VirtualSphere_Oval,_lisp,chem::VirtualSphere_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::VirtualSphere_O>::id,chem::VirtualSphere_O::static_classSymbol());
    chem::VirtualSphere_O::___staticClass = classchem__VirtualSphere_Oval;
#ifdef USE_MPS
    chem::VirtualSphere_O::static_Kind = gctools::GCKind<chem::VirtualSphere_O>::Kind;
#endif
    core::af_setf_findClass(classchem__VirtualSphere_Oval,chem::VirtualSphere_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::VirtualSphere_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::VirtualSphere_O>>::allocateClass();
        chem::VirtualSphere_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::VirtualSphere_O::static_className() % (void*)(chem::VirtualSphere_O::static_allocator) );
    classchem__VirtualSphere_Oval->setCreator(chem::VirtualSphere_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::VirtualSphere_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__VirtualSphere_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__WildElementDict_Oval]"));
    core::BuiltInClass_sp classchem__WildElementDict_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__WildElementDict_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__WildElementDict_Oval,_lisp,chem::WildElementDict_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::WildElementDict_O>::id,chem::WildElementDict_O::static_classSymbol());
    chem::WildElementDict_O::___staticClass = classchem__WildElementDict_Oval;
#ifdef USE_MPS
    chem::WildElementDict_O::static_Kind = gctools::GCKind<chem::WildElementDict_O>::Kind;
#endif
    core::af_setf_findClass(classchem__WildElementDict_Oval,chem::WildElementDict_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::WildElementDict_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::WildElementDict_O>>::allocateClass();
        chem::WildElementDict_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::WildElementDict_O::static_className() % (void*)(chem::WildElementDict_O::static_allocator) );
    classchem__WildElementDict_Oval->setCreator(chem::WildElementDict_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::WildElementDict_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__WildElementDict_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ZMatrixEntry_Oval]"));
    core::BuiltInClass_sp classchem__ZMatrixEntry_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ZMatrixEntry_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ZMatrixEntry_Oval,_lisp,chem::ZMatrixEntry_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ZMatrixEntry_O>::id,chem::ZMatrixEntry_O::static_classSymbol());
    chem::ZMatrixEntry_O::___staticClass = classchem__ZMatrixEntry_Oval;
#ifdef USE_MPS
    chem::ZMatrixEntry_O::static_Kind = gctools::GCKind<chem::ZMatrixEntry_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ZMatrixEntry_Oval,chem::ZMatrixEntry_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ZMatrixEntry_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ZMatrixEntry_O>>::allocateClass();
        chem::ZMatrixEntry_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ZMatrixEntry_O::static_className() % (void*)(chem::ZMatrixEntry_O::static_allocator) );
    classchem__ZMatrixEntry_Oval->setCreator(chem::ZMatrixEntry_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ZMatrixEntry_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ZMatrixEntry_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ZMatrixInternal_Oval]"));
    core::BuiltInClass_sp classchem__ZMatrixInternal_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ZMatrixInternal_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ZMatrixInternal_Oval,_lisp,chem::ZMatrixInternal_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ZMatrixInternal_O>::id,chem::ZMatrixInternal_O::static_classSymbol());
    chem::ZMatrixInternal_O::___staticClass = classchem__ZMatrixInternal_Oval;
#ifdef USE_MPS
    chem::ZMatrixInternal_O::static_Kind = gctools::GCKind<chem::ZMatrixInternal_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ZMatrixInternal_Oval,chem::ZMatrixInternal_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ZMatrixInternal_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ZMatrixInternal_O>>::allocateClass();
        chem::ZMatrixInternal_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ZMatrixInternal_O::static_className() % (void*)(chem::ZMatrixInternal_O::static_allocator) );
    classchem__ZMatrixInternal_Oval->setCreator(chem::ZMatrixInternal_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ZMatrixInternal_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ZMatrixInternal_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ZMatrix_Oval]"));
    core::BuiltInClass_sp classchem__ZMatrix_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ZMatrix_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ZMatrix_Oval,_lisp,chem::ZMatrix_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ZMatrix_O>::id,chem::ZMatrix_O::static_classSymbol());
    chem::ZMatrix_O::___staticClass = classchem__ZMatrix_Oval;
#ifdef USE_MPS
    chem::ZMatrix_O::static_Kind = gctools::GCKind<chem::ZMatrix_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ZMatrix_Oval,chem::ZMatrix_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ZMatrix_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ZMatrix_O>>::allocateClass();
        chem::ZMatrix_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ZMatrix_O::static_className() % (void*)(chem::ZMatrix_O::static_allocator) );
    classchem__ZMatrix_Oval->setCreator(chem::ZMatrix_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ZMatrix_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ZMatrix_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Aggregate_Oval]"));
    core::BuiltInClass_sp classchem__Aggregate_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Aggregate_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Aggregate_Oval,_lisp,chem::Aggregate_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Aggregate_O>::id,chem::Aggregate_O::static_classSymbol());
    chem::Aggregate_O::___staticClass = classchem__Aggregate_Oval;
#ifdef USE_MPS
    chem::Aggregate_O::static_Kind = gctools::GCKind<chem::Aggregate_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Aggregate_Oval,chem::Aggregate_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Aggregate_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Aggregate_O>>::allocateClass();
        chem::Aggregate_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Aggregate_O::static_className() % (void*)(chem::Aggregate_O::static_allocator) );
    classchem__Aggregate_Oval->setCreator(chem::Aggregate_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Aggregate_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Aggregate_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomReference_Oval]"));
    core::BuiltInClass_sp classchem__AtomReference_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomReference_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomReference_Oval,_lisp,chem::AtomReference_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomReference_O>::id,chem::AtomReference_O::static_classSymbol());
    chem::AtomReference_O::___staticClass = classchem__AtomReference_Oval;
#ifdef USE_MPS
    chem::AtomReference_O::static_Kind = gctools::GCKind<chem::AtomReference_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomReference_Oval,chem::AtomReference_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomReference_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomReference_O>>::allocateClass();
        chem::AtomReference_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomReference_O::static_className() % (void*)(chem::AtomReference_O::static_allocator) );
    classchem__AtomReference_Oval->setCreator(chem::AtomReference_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomReference_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomReference_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Atom_Oval]"));
    core::BuiltInClass_sp classchem__Atom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Atom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Atom_Oval,_lisp,chem::Atom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Atom_O>::id,chem::Atom_O::static_classSymbol());
    chem::Atom_O::___staticClass = classchem__Atom_Oval;
#ifdef USE_MPS
    chem::Atom_O::static_Kind = gctools::GCKind<chem::Atom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Atom_Oval,chem::Atom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Atom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Atom_O>>::allocateClass();
        chem::Atom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Atom_O::static_className() % (void*)(chem::Atom_O::static_allocator) );
    classchem__Atom_Oval->setCreator(chem::Atom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Atom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Atom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__BondListMatchNode_Oval]"));
    core::BuiltInClass_sp classchem__BondListMatchNode_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__BondListMatchNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__BondListMatchNode_Oval,_lisp,chem::BondListMatchNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::BondListMatchNode_O>::id,chem::BondListMatchNode_O::static_classSymbol());
    chem::BondListMatchNode_O::___staticClass = classchem__BondListMatchNode_Oval;
#ifdef USE_MPS
    chem::BondListMatchNode_O::static_Kind = gctools::GCKind<chem::BondListMatchNode_O>::Kind;
#endif
    core::af_setf_findClass(classchem__BondListMatchNode_Oval,chem::BondListMatchNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::BondListMatchNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::BondListMatchNode_O>>::allocateClass();
        chem::BondListMatchNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::BondListMatchNode_O::static_className() % (void*)(chem::BondListMatchNode_O::static_allocator) );
    classchem__BondListMatchNode_Oval->setCreator(chem::BondListMatchNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::BondListMatchNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__BondListMatchNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__BondMatchNode_Oval]"));
    core::BuiltInClass_sp classchem__BondMatchNode_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__BondMatchNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__BondMatchNode_Oval,_lisp,chem::BondMatchNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::BondMatchNode_O>::id,chem::BondMatchNode_O::static_classSymbol());
    chem::BondMatchNode_O::___staticClass = classchem__BondMatchNode_Oval;
#ifdef USE_MPS
    chem::BondMatchNode_O::static_Kind = gctools::GCKind<chem::BondMatchNode_O>::Kind;
#endif
    core::af_setf_findClass(classchem__BondMatchNode_Oval,chem::BondMatchNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::BondMatchNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::BondMatchNode_O>>::allocateClass();
        chem::BondMatchNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::BondMatchNode_O::static_className() % (void*)(chem::BondMatchNode_O::static_allocator) );
    classchem__BondMatchNode_Oval->setCreator(chem::BondMatchNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::BondMatchNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__BondMatchNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CalculatePositionAlongBond_Oval]"));
    core::BuiltInClass_sp classchem__CalculatePositionAlongBond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CalculatePositionAlongBond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CalculatePositionAlongBond_Oval,_lisp,chem::CalculatePositionAlongBond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CalculatePositionAlongBond_O>::id,chem::CalculatePositionAlongBond_O::static_classSymbol());
    chem::CalculatePositionAlongBond_O::___staticClass = classchem__CalculatePositionAlongBond_Oval;
#ifdef USE_MPS
    chem::CalculatePositionAlongBond_O::static_Kind = gctools::GCKind<chem::CalculatePositionAlongBond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CalculatePositionAlongBond_Oval,chem::CalculatePositionAlongBond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CalculatePositionAlongBond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CalculatePositionAlongBond_O>>::allocateClass();
        chem::CalculatePositionAlongBond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CalculatePositionAlongBond_O::static_className() % (void*)(chem::CalculatePositionAlongBond_O::static_allocator) );
    classchem__CalculatePositionAlongBond_Oval->setCreator(chem::CalculatePositionAlongBond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CalculatePositionAlongBond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CalculatePositionAlongBond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CalculatePositionRelativeToOrigin_Oval]"));
    core::BuiltInClass_sp classchem__CalculatePositionRelativeToOrigin_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CalculatePositionRelativeToOrigin_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CalculatePositionRelativeToOrigin_Oval,_lisp,chem::CalculatePositionRelativeToOrigin_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CalculatePositionRelativeToOrigin_O>::id,chem::CalculatePositionRelativeToOrigin_O::static_classSymbol());
    chem::CalculatePositionRelativeToOrigin_O::___staticClass = classchem__CalculatePositionRelativeToOrigin_Oval;
#ifdef USE_MPS
    chem::CalculatePositionRelativeToOrigin_O::static_Kind = gctools::GCKind<chem::CalculatePositionRelativeToOrigin_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CalculatePositionRelativeToOrigin_Oval,chem::CalculatePositionRelativeToOrigin_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CalculatePositionRelativeToOrigin_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CalculatePositionRelativeToOrigin_O>>::allocateClass();
        chem::CalculatePositionRelativeToOrigin_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CalculatePositionRelativeToOrigin_O::static_className() % (void*)(chem::CalculatePositionRelativeToOrigin_O::static_allocator) );
    classchem__CalculatePositionRelativeToOrigin_Oval->setCreator(chem::CalculatePositionRelativeToOrigin_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CalculatePositionRelativeToOrigin_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CalculatePositionRelativeToOrigin_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__CalculatePositionUsingInternals_Oval]"));
    core::BuiltInClass_sp classchem__CalculatePositionUsingInternals_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__CalculatePositionUsingInternals_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__CalculatePositionUsingInternals_Oval,_lisp,chem::CalculatePositionUsingInternals_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::CalculatePositionUsingInternals_O>::id,chem::CalculatePositionUsingInternals_O::static_classSymbol());
    chem::CalculatePositionUsingInternals_O::___staticClass = classchem__CalculatePositionUsingInternals_Oval;
#ifdef USE_MPS
    chem::CalculatePositionUsingInternals_O::static_Kind = gctools::GCKind<chem::CalculatePositionUsingInternals_O>::Kind;
#endif
    core::af_setf_findClass(classchem__CalculatePositionUsingInternals_Oval,chem::CalculatePositionUsingInternals_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::CalculatePositionUsingInternals_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::CalculatePositionUsingInternals_O>>::allocateClass();
        chem::CalculatePositionUsingInternals_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::CalculatePositionUsingInternals_O::static_className() % (void*)(chem::CalculatePositionUsingInternals_O::static_allocator) );
    classchem__CalculatePositionUsingInternals_Oval->setCreator(chem::CalculatePositionUsingInternals_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::CalculatePositionUsingInternals_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__CalculatePositionUsingInternals_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ConstitutionVirtualAtom_Oval]"));
    core::BuiltInClass_sp classchem__ConstitutionVirtualAtom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ConstitutionVirtualAtom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ConstitutionVirtualAtom_Oval,_lisp,chem::ConstitutionVirtualAtom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ConstitutionVirtualAtom_O>::id,chem::ConstitutionVirtualAtom_O::static_classSymbol());
    chem::ConstitutionVirtualAtom_O::___staticClass = classchem__ConstitutionVirtualAtom_Oval;
#ifdef USE_MPS
    chem::ConstitutionVirtualAtom_O::static_Kind = gctools::GCKind<chem::ConstitutionVirtualAtom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ConstitutionVirtualAtom_Oval,chem::ConstitutionVirtualAtom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ConstitutionVirtualAtom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ConstitutionVirtualAtom_O>>::allocateClass();
        chem::ConstitutionVirtualAtom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ConstitutionVirtualAtom_O::static_className() % (void*)(chem::ConstitutionVirtualAtom_O::static_allocator) );
    classchem__ConstitutionVirtualAtom_Oval->setCreator(chem::ConstitutionVirtualAtom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ConstitutionVirtualAtom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ConstitutionVirtualAtom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Constitution_Oval]"));
    core::BuiltInClass_sp classchem__Constitution_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Constitution_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Constitution_Oval,_lisp,chem::Constitution_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Constitution_O>::id,chem::Constitution_O::static_classSymbol());
    chem::Constitution_O::___staticClass = classchem__Constitution_Oval;
#ifdef USE_MPS
    chem::Constitution_O::static_Kind = gctools::GCKind<chem::Constitution_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Constitution_Oval,chem::Constitution_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Constitution_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Constitution_O>>::allocateClass();
        chem::Constitution_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Constitution_O::static_className() % (void*)(chem::Constitution_O::static_allocator) );
    classchem__Constitution_Oval->setCreator(chem::Constitution_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Constitution_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Constitution_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__DirectionalCoupling_Oval]"));
    core::BuiltInClass_sp classchem__DirectionalCoupling_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__DirectionalCoupling_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__DirectionalCoupling_Oval,_lisp,chem::DirectionalCoupling_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::DirectionalCoupling_O>::id,chem::DirectionalCoupling_O::static_classSymbol());
    chem::DirectionalCoupling_O::___staticClass = classchem__DirectionalCoupling_Oval;
#ifdef USE_MPS
    chem::DirectionalCoupling_O::static_Kind = gctools::GCKind<chem::DirectionalCoupling_O>::Kind;
#endif
    core::af_setf_findClass(classchem__DirectionalCoupling_Oval,chem::DirectionalCoupling_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::DirectionalCoupling_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::DirectionalCoupling_O>>::allocateClass();
        chem::DirectionalCoupling_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::DirectionalCoupling_O::static_className() % (void*)(chem::DirectionalCoupling_O::static_allocator) );
    classchem__DirectionalCoupling_Oval->setCreator(chem::DirectionalCoupling_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::DirectionalCoupling_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__DirectionalCoupling_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyAnchorRestraint_Oval]"));
    core::BuiltInClass_sp classchem__EnergyAnchorRestraint_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyAnchorRestraint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyAnchorRestraint_Oval,_lisp,chem::EnergyAnchorRestraint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyAnchorRestraint_O>::id,chem::EnergyAnchorRestraint_O::static_classSymbol());
    chem::EnergyAnchorRestraint_O::___staticClass = classchem__EnergyAnchorRestraint_Oval;
#ifdef USE_MPS
    chem::EnergyAnchorRestraint_O::static_Kind = gctools::GCKind<chem::EnergyAnchorRestraint_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyAnchorRestraint_Oval,chem::EnergyAnchorRestraint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyAnchorRestraint_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyAnchorRestraint_O>>::allocateClass();
        chem::EnergyAnchorRestraint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyAnchorRestraint_O::static_className() % (void*)(chem::EnergyAnchorRestraint_O::static_allocator) );
    classchem__EnergyAnchorRestraint_Oval->setCreator(chem::EnergyAnchorRestraint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyAnchorRestraint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyAnchorRestraint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyAngle_Oval]"));
    core::BuiltInClass_sp classchem__EnergyAngle_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyAngle_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyAngle_Oval,_lisp,chem::EnergyAngle_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyAngle_O>::id,chem::EnergyAngle_O::static_classSymbol());
    chem::EnergyAngle_O::___staticClass = classchem__EnergyAngle_Oval;
#ifdef USE_MPS
    chem::EnergyAngle_O::static_Kind = gctools::GCKind<chem::EnergyAngle_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyAngle_Oval,chem::EnergyAngle_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyAngle_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyAngle_O>>::allocateClass();
        chem::EnergyAngle_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyAngle_O::static_className() % (void*)(chem::EnergyAngle_O::static_allocator) );
    classchem__EnergyAngle_Oval->setCreator(chem::EnergyAngle_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyAngle_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyAngle_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyChiralRestraint_Oval]"));
    core::BuiltInClass_sp classchem__EnergyChiralRestraint_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyChiralRestraint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyChiralRestraint_Oval,_lisp,chem::EnergyChiralRestraint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyChiralRestraint_O>::id,chem::EnergyChiralRestraint_O::static_classSymbol());
    chem::EnergyChiralRestraint_O::___staticClass = classchem__EnergyChiralRestraint_Oval;
#ifdef USE_MPS
    chem::EnergyChiralRestraint_O::static_Kind = gctools::GCKind<chem::EnergyChiralRestraint_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyChiralRestraint_Oval,chem::EnergyChiralRestraint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyChiralRestraint_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyChiralRestraint_O>>::allocateClass();
        chem::EnergyChiralRestraint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyChiralRestraint_O::static_className() % (void*)(chem::EnergyChiralRestraint_O::static_allocator) );
    classchem__EnergyChiralRestraint_Oval->setCreator(chem::EnergyChiralRestraint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyChiralRestraint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyChiralRestraint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyDihedral_Oval]"));
    core::BuiltInClass_sp classchem__EnergyDihedral_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyDihedral_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyDihedral_Oval,_lisp,chem::EnergyDihedral_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyDihedral_O>::id,chem::EnergyDihedral_O::static_classSymbol());
    chem::EnergyDihedral_O::___staticClass = classchem__EnergyDihedral_Oval;
#ifdef USE_MPS
    chem::EnergyDihedral_O::static_Kind = gctools::GCKind<chem::EnergyDihedral_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyDihedral_Oval,chem::EnergyDihedral_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyDihedral_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyDihedral_O>>::allocateClass();
        chem::EnergyDihedral_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyDihedral_O::static_className() % (void*)(chem::EnergyDihedral_O::static_allocator) );
    classchem__EnergyDihedral_Oval->setCreator(chem::EnergyDihedral_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyDihedral_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyDihedral_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyFixedNonbondRestraint_Oval]"));
    core::BuiltInClass_sp classchem__EnergyFixedNonbondRestraint_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyFixedNonbondRestraint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyFixedNonbondRestraint_Oval,_lisp,chem::EnergyFixedNonbondRestraint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyFixedNonbondRestraint_O>::id,chem::EnergyFixedNonbondRestraint_O::static_classSymbol());
    chem::EnergyFixedNonbondRestraint_O::___staticClass = classchem__EnergyFixedNonbondRestraint_Oval;
#ifdef USE_MPS
    chem::EnergyFixedNonbondRestraint_O::static_Kind = gctools::GCKind<chem::EnergyFixedNonbondRestraint_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyFixedNonbondRestraint_Oval,chem::EnergyFixedNonbondRestraint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyFixedNonbondRestraint_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyFixedNonbondRestraint_O>>::allocateClass();
        chem::EnergyFixedNonbondRestraint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyFixedNonbondRestraint_O::static_className() % (void*)(chem::EnergyFixedNonbondRestraint_O::static_allocator) );
    classchem__EnergyFixedNonbondRestraint_Oval->setCreator(chem::EnergyFixedNonbondRestraint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyFixedNonbondRestraint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyFixedNonbondRestraint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyImproperRestraint_Oval]"));
    core::BuiltInClass_sp classchem__EnergyImproperRestraint_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyImproperRestraint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyImproperRestraint_Oval,_lisp,chem::EnergyImproperRestraint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyImproperRestraint_O>::id,chem::EnergyImproperRestraint_O::static_classSymbol());
    chem::EnergyImproperRestraint_O::___staticClass = classchem__EnergyImproperRestraint_Oval;
#ifdef USE_MPS
    chem::EnergyImproperRestraint_O::static_Kind = gctools::GCKind<chem::EnergyImproperRestraint_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyImproperRestraint_Oval,chem::EnergyImproperRestraint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyImproperRestraint_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyImproperRestraint_O>>::allocateClass();
        chem::EnergyImproperRestraint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyImproperRestraint_O::static_className() % (void*)(chem::EnergyImproperRestraint_O::static_allocator) );
    classchem__EnergyImproperRestraint_Oval->setCreator(chem::EnergyImproperRestraint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyImproperRestraint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyImproperRestraint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyNonbond_Oval]"));
    core::BuiltInClass_sp classchem__EnergyNonbond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyNonbond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyNonbond_Oval,_lisp,chem::EnergyNonbond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyNonbond_O>::id,chem::EnergyNonbond_O::static_classSymbol());
    chem::EnergyNonbond_O::___staticClass = classchem__EnergyNonbond_Oval;
#ifdef USE_MPS
    chem::EnergyNonbond_O::static_Kind = gctools::GCKind<chem::EnergyNonbond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyNonbond_Oval,chem::EnergyNonbond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyNonbond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyNonbond_O>>::allocateClass();
        chem::EnergyNonbond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyNonbond_O::static_className() % (void*)(chem::EnergyNonbond_O::static_allocator) );
    classchem__EnergyNonbond_Oval->setCreator(chem::EnergyNonbond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyNonbond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyNonbond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EnergyStretch_Oval]"));
    core::BuiltInClass_sp classchem__EnergyStretch_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EnergyStretch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EnergyStretch_Oval,_lisp,chem::EnergyStretch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EnergyStretch_O>::id,chem::EnergyStretch_O::static_classSymbol());
    chem::EnergyStretch_O::___staticClass = classchem__EnergyStretch_Oval;
#ifdef USE_MPS
    chem::EnergyStretch_O::static_Kind = gctools::GCKind<chem::EnergyStretch_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EnergyStretch_Oval,chem::EnergyStretch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EnergyStretch_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EnergyStretch_O>>::allocateClass();
        chem::EnergyStretch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EnergyStretch_O::static_className() % (void*)(chem::EnergyStretch_O::static_allocator) );
    classchem__EnergyStretch_Oval->setCreator(chem::EnergyStretch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EnergyStretch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EnergyStretch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EntityNameSetBase_Oval]"));
    core::BuiltInClass_sp classchem__EntityNameSetBase_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EntityNameSetBase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EntityNameSetBase_Oval,_lisp,chem::EntityNameSetBase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EntityNameSetBase_O>::id,chem::EntityNameSetBase_O::static_classSymbol());
    chem::EntityNameSetBase_O::___staticClass = classchem__EntityNameSetBase_Oval;
#ifdef USE_MPS
    chem::EntityNameSetBase_O::static_Kind = gctools::GCKind<chem::EntityNameSetBase_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EntityNameSetBase_Oval,chem::EntityNameSetBase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EntityNameSetBase_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EntityNameSetBase_O>>::allocateClass();
        chem::EntityNameSetBase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EntityNameSetBase_O::static_className() % (void*)(chem::EntityNameSetBase_O::static_allocator) );
    classchem__EntityNameSetBase_Oval->setCreator(chem::EntityNameSetBase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EntityNameSetBase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EntityNameSetBase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFAngleDb_Oval]"));
    core::BuiltInClass_sp classchem__FFAngleDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFAngleDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFAngleDb_Oval,_lisp,chem::FFAngleDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFAngleDb_O>::id,chem::FFAngleDb_O::static_classSymbol());
    chem::FFAngleDb_O::___staticClass = classchem__FFAngleDb_Oval;
#ifdef USE_MPS
    chem::FFAngleDb_O::static_Kind = gctools::GCKind<chem::FFAngleDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFAngleDb_Oval,chem::FFAngleDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFAngleDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFAngleDb_O>>::allocateClass();
        chem::FFAngleDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFAngleDb_O::static_className() % (void*)(chem::FFAngleDb_O::static_allocator) );
    classchem__FFAngleDb_Oval->setCreator(chem::FFAngleDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFAngleDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFAngleDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFAngle_Oval]"));
    core::BuiltInClass_sp classchem__FFAngle_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFAngle_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFAngle_Oval,_lisp,chem::FFAngle_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFAngle_O>::id,chem::FFAngle_O::static_classSymbol());
    chem::FFAngle_O::___staticClass = classchem__FFAngle_Oval;
#ifdef USE_MPS
    chem::FFAngle_O::static_Kind = gctools::GCKind<chem::FFAngle_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFAngle_Oval,chem::FFAngle_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFAngle_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFAngle_O>>::allocateClass();
        chem::FFAngle_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFAngle_O::static_className() % (void*)(chem::FFAngle_O::static_allocator) );
    classchem__FFAngle_Oval->setCreator(chem::FFAngle_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFAngle_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFAngle_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFItorDb_Oval]"));
    core::BuiltInClass_sp classchem__FFItorDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFItorDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFItorDb_Oval,_lisp,chem::FFItorDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFItorDb_O>::id,chem::FFItorDb_O::static_classSymbol());
    chem::FFItorDb_O::___staticClass = classchem__FFItorDb_Oval;
#ifdef USE_MPS
    chem::FFItorDb_O::static_Kind = gctools::GCKind<chem::FFItorDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFItorDb_Oval,chem::FFItorDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFItorDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFItorDb_O>>::allocateClass();
        chem::FFItorDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFItorDb_O::static_className() % (void*)(chem::FFItorDb_O::static_allocator) );
    classchem__FFItorDb_Oval->setCreator(chem::FFItorDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFItorDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFItorDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFItor_Oval]"));
    core::BuiltInClass_sp classchem__FFItor_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFItor_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFItor_Oval,_lisp,chem::FFItor_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFItor_O>::id,chem::FFItor_O::static_classSymbol());
    chem::FFItor_O::___staticClass = classchem__FFItor_Oval;
#ifdef USE_MPS
    chem::FFItor_O::static_Kind = gctools::GCKind<chem::FFItor_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFItor_Oval,chem::FFItor_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFItor_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFItor_O>>::allocateClass();
        chem::FFItor_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFItor_O::static_className() % (void*)(chem::FFItor_O::static_allocator) );
    classchem__FFItor_Oval->setCreator(chem::FFItor_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFItor_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFItor_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFNonbondDb_Oval]"));
    core::BuiltInClass_sp classchem__FFNonbondDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFNonbondDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFNonbondDb_Oval,_lisp,chem::FFNonbondDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFNonbondDb_O>::id,chem::FFNonbondDb_O::static_classSymbol());
    chem::FFNonbondDb_O::___staticClass = classchem__FFNonbondDb_Oval;
#ifdef USE_MPS
    chem::FFNonbondDb_O::static_Kind = gctools::GCKind<chem::FFNonbondDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFNonbondDb_Oval,chem::FFNonbondDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFNonbondDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFNonbondDb_O>>::allocateClass();
        chem::FFNonbondDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFNonbondDb_O::static_className() % (void*)(chem::FFNonbondDb_O::static_allocator) );
    classchem__FFNonbondDb_Oval->setCreator(chem::FFNonbondDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFNonbondDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFNonbondDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFNonbond_Oval]"));
    core::BuiltInClass_sp classchem__FFNonbond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFNonbond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFNonbond_Oval,_lisp,chem::FFNonbond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFNonbond_O>::id,chem::FFNonbond_O::static_classSymbol());
    chem::FFNonbond_O::___staticClass = classchem__FFNonbond_Oval;
#ifdef USE_MPS
    chem::FFNonbond_O::static_Kind = gctools::GCKind<chem::FFNonbond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFNonbond_Oval,chem::FFNonbond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFNonbond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFNonbond_O>>::allocateClass();
        chem::FFNonbond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFNonbond_O::static_className() % (void*)(chem::FFNonbond_O::static_allocator) );
    classchem__FFNonbond_Oval->setCreator(chem::FFNonbond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFNonbond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFNonbond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFPtorDb_Oval]"));
    core::BuiltInClass_sp classchem__FFPtorDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFPtorDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFPtorDb_Oval,_lisp,chem::FFPtorDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFPtorDb_O>::id,chem::FFPtorDb_O::static_classSymbol());
    chem::FFPtorDb_O::___staticClass = classchem__FFPtorDb_Oval;
#ifdef USE_MPS
    chem::FFPtorDb_O::static_Kind = gctools::GCKind<chem::FFPtorDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFPtorDb_Oval,chem::FFPtorDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFPtorDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFPtorDb_O>>::allocateClass();
        chem::FFPtorDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFPtorDb_O::static_className() % (void*)(chem::FFPtorDb_O::static_allocator) );
    classchem__FFPtorDb_Oval->setCreator(chem::FFPtorDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFPtorDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFPtorDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFPtor_Oval]"));
    core::BuiltInClass_sp classchem__FFPtor_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFPtor_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFPtor_Oval,_lisp,chem::FFPtor_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFPtor_O>::id,chem::FFPtor_O::static_classSymbol());
    chem::FFPtor_O::___staticClass = classchem__FFPtor_Oval;
#ifdef USE_MPS
    chem::FFPtor_O::static_Kind = gctools::GCKind<chem::FFPtor_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFPtor_Oval,chem::FFPtor_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFPtor_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFPtor_O>>::allocateClass();
        chem::FFPtor_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFPtor_O::static_className() % (void*)(chem::FFPtor_O::static_allocator) );
    classchem__FFPtor_Oval->setCreator(chem::FFPtor_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFPtor_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFPtor_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFStretchDb_Oval]"));
    core::BuiltInClass_sp classchem__FFStretchDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFStretchDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFStretchDb_Oval,_lisp,chem::FFStretchDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFStretchDb_O>::id,chem::FFStretchDb_O::static_classSymbol());
    chem::FFStretchDb_O::___staticClass = classchem__FFStretchDb_Oval;
#ifdef USE_MPS
    chem::FFStretchDb_O::static_Kind = gctools::GCKind<chem::FFStretchDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFStretchDb_Oval,chem::FFStretchDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFStretchDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFStretchDb_O>>::allocateClass();
        chem::FFStretchDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFStretchDb_O::static_className() % (void*)(chem::FFStretchDb_O::static_allocator) );
    classchem__FFStretchDb_Oval->setCreator(chem::FFStretchDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFStretchDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFStretchDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFStretch_Oval]"));
    core::BuiltInClass_sp classchem__FFStretch_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFStretch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFStretch_Oval,_lisp,chem::FFStretch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFStretch_O>::id,chem::FFStretch_O::static_classSymbol());
    chem::FFStretch_O::___staticClass = classchem__FFStretch_Oval;
#ifdef USE_MPS
    chem::FFStretch_O::static_Kind = gctools::GCKind<chem::FFStretch_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFStretch_Oval,chem::FFStretch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFStretch_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFStretch_O>>::allocateClass();
        chem::FFStretch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFStretch_O::static_className() % (void*)(chem::FFStretch_O::static_allocator) );
    classchem__FFStretch_Oval->setCreator(chem::FFStretch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFStretch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFStretch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFTypesDb_Oval]"));
    core::BuiltInClass_sp classchem__FFTypesDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFTypesDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFTypesDb_Oval,_lisp,chem::FFTypesDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFTypesDb_O>::id,chem::FFTypesDb_O::static_classSymbol());
    chem::FFTypesDb_O::___staticClass = classchem__FFTypesDb_Oval;
#ifdef USE_MPS
    chem::FFTypesDb_O::static_Kind = gctools::GCKind<chem::FFTypesDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFTypesDb_Oval,chem::FFTypesDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFTypesDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFTypesDb_O>>::allocateClass();
        chem::FFTypesDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFTypesDb_O::static_className() % (void*)(chem::FFTypesDb_O::static_allocator) );
    classchem__FFTypesDb_Oval->setCreator(chem::FFTypesDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFTypesDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFTypesDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFVdwDb_Oval]"));
    core::BuiltInClass_sp classchem__FFVdwDb_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFVdwDb_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFVdwDb_Oval,_lisp,chem::FFVdwDb_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFVdwDb_O>::id,chem::FFVdwDb_O::static_classSymbol());
    chem::FFVdwDb_O::___staticClass = classchem__FFVdwDb_Oval;
#ifdef USE_MPS
    chem::FFVdwDb_O::static_Kind = gctools::GCKind<chem::FFVdwDb_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFVdwDb_Oval,chem::FFVdwDb_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFVdwDb_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFVdwDb_O>>::allocateClass();
        chem::FFVdwDb_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFVdwDb_O::static_className() % (void*)(chem::FFVdwDb_O::static_allocator) );
    classchem__FFVdwDb_Oval->setCreator(chem::FFVdwDb_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFVdwDb_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFVdwDb_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FFVdw_Oval]"));
    core::BuiltInClass_sp classchem__FFVdw_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FFVdw_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FFVdw_Oval,_lisp,chem::FFVdw_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FFVdw_O>::id,chem::FFVdw_O::static_classSymbol());
    chem::FFVdw_O::___staticClass = classchem__FFVdw_Oval;
#ifdef USE_MPS
    chem::FFVdw_O::static_Kind = gctools::GCKind<chem::FFVdw_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FFVdw_Oval,chem::FFVdw_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FFVdw_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FFVdw_O>>::allocateClass();
        chem::FFVdw_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FFVdw_O::static_className() % (void*)(chem::FFVdw_O::static_allocator) );
    classchem__FFVdw_Oval->setCreator(chem::FFVdw_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FFVdw_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FFVdw_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ForceMatchReport_Oval]"));
    core::BuiltInClass_sp classchem__ForceMatchReport_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ForceMatchReport_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ForceMatchReport_Oval,_lisp,chem::ForceMatchReport_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ForceMatchReport_O>::id,chem::ForceMatchReport_O::static_classSymbol());
    chem::ForceMatchReport_O::___staticClass = classchem__ForceMatchReport_Oval;
#ifdef USE_MPS
    chem::ForceMatchReport_O::static_Kind = gctools::GCKind<chem::ForceMatchReport_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ForceMatchReport_Oval,chem::ForceMatchReport_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ForceMatchReport_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ForceMatchReport_O>>::allocateClass();
        chem::ForceMatchReport_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ForceMatchReport_O::static_className() % (void*)(chem::ForceMatchReport_O::static_allocator) );
    classchem__ForceMatchReport_Oval->setCreator(chem::ForceMatchReport_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ForceMatchReport_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ForceMatchReport_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__FullLargeSquareMatrix_Oval]"));
    core::BuiltInClass_sp classchem__FullLargeSquareMatrix_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__FullLargeSquareMatrix_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__FullLargeSquareMatrix_Oval,_lisp,chem::FullLargeSquareMatrix_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::FullLargeSquareMatrix_O>::id,chem::FullLargeSquareMatrix_O::static_classSymbol());
    chem::FullLargeSquareMatrix_O::___staticClass = classchem__FullLargeSquareMatrix_Oval;
#ifdef USE_MPS
    chem::FullLargeSquareMatrix_O::static_Kind = gctools::GCKind<chem::FullLargeSquareMatrix_O>::Kind;
#endif
    core::af_setf_findClass(classchem__FullLargeSquareMatrix_Oval,chem::FullLargeSquareMatrix_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::FullLargeSquareMatrix_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::FullLargeSquareMatrix_O>>::allocateClass();
        chem::FullLargeSquareMatrix_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::FullLargeSquareMatrix_O::static_className() % (void*)(chem::FullLargeSquareMatrix_O::static_allocator) );
    classchem__FullLargeSquareMatrix_Oval->setCreator(chem::FullLargeSquareMatrix_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::FullLargeSquareMatrix_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__FullLargeSquareMatrix_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__IterateMatter_Oval]"));
    core::BuiltInClass_sp classchem__IterateMatter_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__IterateMatter_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__IterateMatter_Oval,_lisp,chem::IterateMatter_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::IterateMatter_O>::id,chem::IterateMatter_O::static_classSymbol());
    chem::IterateMatter_O::___staticClass = classchem__IterateMatter_Oval;
#ifdef USE_MPS
    chem::IterateMatter_O::static_Kind = gctools::GCKind<chem::IterateMatter_O>::Kind;
#endif
    core::af_setf_findClass(classchem__IterateMatter_Oval,chem::IterateMatter_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::IterateMatter_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::IterateMatter_O>>::allocateClass();
        chem::IterateMatter_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::IterateMatter_O::static_className() % (void*)(chem::IterateMatter_O::static_allocator) );
    classchem__IterateMatter_Oval->setCreator(chem::IterateMatter_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::IterateMatter_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__IterateMatter_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__IterateRestraints_Oval]"));
    core::BuiltInClass_sp classchem__IterateRestraints_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__IterateRestraints_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__IterateRestraints_Oval,_lisp,chem::IterateRestraints_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::IterateRestraints_O>::id,chem::IterateRestraints_O::static_classSymbol());
    chem::IterateRestraints_O::___staticClass = classchem__IterateRestraints_Oval;
#ifdef USE_MPS
    chem::IterateRestraints_O::static_Kind = gctools::GCKind<chem::IterateRestraints_O>::Kind;
#endif
    core::af_setf_findClass(classchem__IterateRestraints_Oval,chem::IterateRestraints_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::IterateRestraints_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::IterateRestraints_O>>::allocateClass();
        chem::IterateRestraints_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::IterateRestraints_O::static_className() % (void*)(chem::IterateRestraints_O::static_allocator) );
    classchem__IterateRestraints_Oval->setCreator(chem::IterateRestraints_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::IterateRestraints_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__IterateRestraints_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__JumpPlug_Oval]"));
    core::BuiltInClass_sp classchem__JumpPlug_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__JumpPlug_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__JumpPlug_Oval,_lisp,chem::JumpPlug_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::JumpPlug_O>::id,chem::JumpPlug_O::static_classSymbol());
    chem::JumpPlug_O::___staticClass = classchem__JumpPlug_Oval;
#ifdef USE_MPS
    chem::JumpPlug_O::static_Kind = gctools::GCKind<chem::JumpPlug_O>::Kind;
#endif
    core::af_setf_findClass(classchem__JumpPlug_Oval,chem::JumpPlug_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::JumpPlug_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::JumpPlug_O>>::allocateClass();
        chem::JumpPlug_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::JumpPlug_O::static_className() % (void*)(chem::JumpPlug_O::static_allocator) );
    classchem__JumpPlug_Oval->setCreator(chem::JumpPlug_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::JumpPlug_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__JumpPlug_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MessageReport_Oval]"));
    core::BuiltInClass_sp classchem__MessageReport_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MessageReport_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MessageReport_Oval,_lisp,chem::MessageReport_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MessageReport_O>::id,chem::MessageReport_O::static_classSymbol());
    chem::MessageReport_O::___staticClass = classchem__MessageReport_Oval;
#ifdef USE_MPS
    chem::MessageReport_O::static_Kind = gctools::GCKind<chem::MessageReport_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MessageReport_Oval,chem::MessageReport_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MessageReport_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MessageReport_O>>::allocateClass();
        chem::MessageReport_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MessageReport_O::static_className() % (void*)(chem::MessageReport_O::static_allocator) );
    classchem__MessageReport_Oval->setCreator(chem::MessageReport_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MessageReport_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MessageReport_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Molecule_Oval]"));
    core::BuiltInClass_sp classchem__Molecule_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Molecule_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Molecule_Oval,_lisp,chem::Molecule_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Molecule_O>::id,chem::Molecule_O::static_classSymbol());
    chem::Molecule_O::___staticClass = classchem__Molecule_Oval;
#ifdef USE_MPS
    chem::Molecule_O::static_Kind = gctools::GCKind<chem::Molecule_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Molecule_Oval,chem::Molecule_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Molecule_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Molecule_O>>::allocateClass();
        chem::Molecule_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Molecule_O::static_className() % (void*)(chem::Molecule_O::static_allocator) );
    classchem__Molecule_Oval->setCreator(chem::Molecule_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Molecule_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Molecule_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Monomer_Oval]"));
    core::BuiltInClass_sp classchem__Monomer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Monomer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Monomer_Oval,_lisp,chem::Monomer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Monomer_O>::id,chem::Monomer_O::static_classSymbol());
    chem::Monomer_O::___staticClass = classchem__Monomer_Oval;
#ifdef USE_MPS
    chem::Monomer_O::static_Kind = gctools::GCKind<chem::Monomer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Monomer_Oval,chem::Monomer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Monomer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Monomer_O>>::allocateClass();
        chem::Monomer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Monomer_O::static_className() % (void*)(chem::Monomer_O::static_allocator) );
    classchem__Monomer_Oval->setCreator(chem::Monomer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Monomer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Monomer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__OligomerPart_Link_Oval]"));
    core::BuiltInClass_sp classchem__OligomerPart_Link_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__OligomerPart_Link_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__OligomerPart_Link_Oval,_lisp,chem::OligomerPart_Link_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::OligomerPart_Link_O>::id,chem::OligomerPart_Link_O::static_classSymbol());
    chem::OligomerPart_Link_O::___staticClass = classchem__OligomerPart_Link_Oval;
#ifdef USE_MPS
    chem::OligomerPart_Link_O::static_Kind = gctools::GCKind<chem::OligomerPart_Link_O>::Kind;
#endif
    core::af_setf_findClass(classchem__OligomerPart_Link_Oval,chem::OligomerPart_Link_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::OligomerPart_Link_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::OligomerPart_Link_O>>::allocateClass();
        chem::OligomerPart_Link_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::OligomerPart_Link_O::static_className() % (void*)(chem::OligomerPart_Link_O::static_allocator) );
    classchem__OligomerPart_Link_Oval->setCreator(chem::OligomerPart_Link_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::OligomerPart_Link_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__OligomerPart_Link_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__OligomerPart_Monomer_Oval]"));
    core::BuiltInClass_sp classchem__OligomerPart_Monomer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__OligomerPart_Monomer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__OligomerPart_Monomer_Oval,_lisp,chem::OligomerPart_Monomer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::OligomerPart_Monomer_O>::id,chem::OligomerPart_Monomer_O::static_classSymbol());
    chem::OligomerPart_Monomer_O::___staticClass = classchem__OligomerPart_Monomer_Oval;
#ifdef USE_MPS
    chem::OligomerPart_Monomer_O::static_Kind = gctools::GCKind<chem::OligomerPart_Monomer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__OligomerPart_Monomer_Oval,chem::OligomerPart_Monomer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::OligomerPart_Monomer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::OligomerPart_Monomer_O>>::allocateClass();
        chem::OligomerPart_Monomer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::OligomerPart_Monomer_O::static_className() % (void*)(chem::OligomerPart_Monomer_O::static_allocator) );
    classchem__OligomerPart_Monomer_Oval->setCreator(chem::OligomerPart_Monomer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::OligomerPart_Monomer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__OligomerPart_Monomer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__PlugWithMates_Oval]"));
    core::BuiltInClass_sp classchem__PlugWithMates_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__PlugWithMates_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__PlugWithMates_Oval,_lisp,chem::PlugWithMates_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::PlugWithMates_O>::id,chem::PlugWithMates_O::static_classSymbol());
    chem::PlugWithMates_O::___staticClass = classchem__PlugWithMates_Oval;
#ifdef USE_MPS
    chem::PlugWithMates_O::static_Kind = gctools::GCKind<chem::PlugWithMates_O>::Kind;
#endif
    core::af_setf_findClass(classchem__PlugWithMates_Oval,chem::PlugWithMates_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::PlugWithMates_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::PlugWithMates_O>>::allocateClass();
        chem::PlugWithMates_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::PlugWithMates_O::static_className() % (void*)(chem::PlugWithMates_O::static_allocator) );
    classchem__PlugWithMates_Oval->setCreator(chem::PlugWithMates_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::PlugWithMates_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__PlugWithMates_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Residue_Oval]"));
    core::BuiltInClass_sp classchem__Residue_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Residue_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Residue_Oval,_lisp,chem::Residue_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Residue_O>::id,chem::Residue_O::static_classSymbol());
    chem::Residue_O::___staticClass = classchem__Residue_Oval;
#ifdef USE_MPS
    chem::Residue_O::static_Kind = gctools::GCKind<chem::Residue_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Residue_Oval,chem::Residue_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Residue_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Residue_O>>::allocateClass();
        chem::Residue_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Residue_O::static_className() % (void*)(chem::Residue_O::static_allocator) );
    classchem__Residue_Oval->setCreator(chem::Residue_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Residue_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Residue_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestrainedExoCyclicAtom_Oval]"));
    core::BuiltInClass_sp classchem__RestrainedExoCyclicAtom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestrainedExoCyclicAtom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestrainedExoCyclicAtom_Oval,_lisp,chem::RestrainedExoCyclicAtom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestrainedExoCyclicAtom_O>::id,chem::RestrainedExoCyclicAtom_O::static_classSymbol());
    chem::RestrainedExoCyclicAtom_O::___staticClass = classchem__RestrainedExoCyclicAtom_Oval;
#ifdef USE_MPS
    chem::RestrainedExoCyclicAtom_O::static_Kind = gctools::GCKind<chem::RestrainedExoCyclicAtom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestrainedExoCyclicAtom_Oval,chem::RestrainedExoCyclicAtom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestrainedExoCyclicAtom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestrainedExoCyclicAtom_O>>::allocateClass();
        chem::RestrainedExoCyclicAtom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestrainedExoCyclicAtom_O::static_className() % (void*)(chem::RestrainedExoCyclicAtom_O::static_allocator) );
    classchem__RestrainedExoCyclicAtom_Oval->setCreator(chem::RestrainedExoCyclicAtom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestrainedExoCyclicAtom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestrainedExoCyclicAtom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestrainedPiBond_Oval]"));
    core::BuiltInClass_sp classchem__RestrainedPiBond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestrainedPiBond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestrainedPiBond_Oval,_lisp,chem::RestrainedPiBond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestrainedPiBond_O>::id,chem::RestrainedPiBond_O::static_classSymbol());
    chem::RestrainedPiBond_O::___staticClass = classchem__RestrainedPiBond_Oval;
#ifdef USE_MPS
    chem::RestrainedPiBond_O::static_Kind = gctools::GCKind<chem::RestrainedPiBond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestrainedPiBond_Oval,chem::RestrainedPiBond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestrainedPiBond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestrainedPiBond_O>>::allocateClass();
        chem::RestrainedPiBond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestrainedPiBond_O::static_className() % (void*)(chem::RestrainedPiBond_O::static_allocator) );
    classchem__RestrainedPiBond_Oval->setCreator(chem::RestrainedPiBond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestrainedPiBond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestrainedPiBond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintAnchor_Oval]"));
    core::BuiltInClass_sp classchem__RestraintAnchor_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintAnchor_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintAnchor_Oval,_lisp,chem::RestraintAnchor_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintAnchor_O>::id,chem::RestraintAnchor_O::static_classSymbol());
    chem::RestraintAnchor_O::___staticClass = classchem__RestraintAnchor_Oval;
#ifdef USE_MPS
    chem::RestraintAnchor_O::static_Kind = gctools::GCKind<chem::RestraintAnchor_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintAnchor_Oval,chem::RestraintAnchor_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintAnchor_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintAnchor_O>>::allocateClass();
        chem::RestraintAnchor_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintAnchor_O::static_className() % (void*)(chem::RestraintAnchor_O::static_allocator) );
    classchem__RestraintAnchor_Oval->setCreator(chem::RestraintAnchor_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintAnchor_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintAnchor_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintAngle_Oval]"));
    core::BuiltInClass_sp classchem__RestraintAngle_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintAngle_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintAngle_Oval,_lisp,chem::RestraintAngle_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintAngle_O>::id,chem::RestraintAngle_O::static_classSymbol());
    chem::RestraintAngle_O::___staticClass = classchem__RestraintAngle_Oval;
#ifdef USE_MPS
    chem::RestraintAngle_O::static_Kind = gctools::GCKind<chem::RestraintAngle_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintAngle_Oval,chem::RestraintAngle_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintAngle_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintAngle_O>>::allocateClass();
        chem::RestraintAngle_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintAngle_O::static_className() % (void*)(chem::RestraintAngle_O::static_allocator) );
    classchem__RestraintAngle_Oval->setCreator(chem::RestraintAngle_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintAngle_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintAngle_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintChiral_Oval]"));
    core::BuiltInClass_sp classchem__RestraintChiral_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintChiral_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintChiral_Oval,_lisp,chem::RestraintChiral_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintChiral_O>::id,chem::RestraintChiral_O::static_classSymbol());
    chem::RestraintChiral_O::___staticClass = classchem__RestraintChiral_Oval;
#ifdef USE_MPS
    chem::RestraintChiral_O::static_Kind = gctools::GCKind<chem::RestraintChiral_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintChiral_Oval,chem::RestraintChiral_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintChiral_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintChiral_O>>::allocateClass();
        chem::RestraintChiral_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintChiral_O::static_className() % (void*)(chem::RestraintChiral_O::static_allocator) );
    classchem__RestraintChiral_Oval->setCreator(chem::RestraintChiral_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintChiral_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintChiral_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintDihedral_Oval]"));
    core::BuiltInClass_sp classchem__RestraintDihedral_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintDihedral_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintDihedral_Oval,_lisp,chem::RestraintDihedral_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintDihedral_O>::id,chem::RestraintDihedral_O::static_classSymbol());
    chem::RestraintDihedral_O::___staticClass = classchem__RestraintDihedral_Oval;
#ifdef USE_MPS
    chem::RestraintDihedral_O::static_Kind = gctools::GCKind<chem::RestraintDihedral_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintDihedral_Oval,chem::RestraintDihedral_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintDihedral_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintDihedral_O>>::allocateClass();
        chem::RestraintDihedral_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintDihedral_O::static_className() % (void*)(chem::RestraintDihedral_O::static_allocator) );
    classchem__RestraintDihedral_Oval->setCreator(chem::RestraintDihedral_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintDihedral_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintDihedral_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintDistance_Oval]"));
    core::BuiltInClass_sp classchem__RestraintDistance_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintDistance_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintDistance_Oval,_lisp,chem::RestraintDistance_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintDistance_O>::id,chem::RestraintDistance_O::static_classSymbol());
    chem::RestraintDistance_O::___staticClass = classchem__RestraintDistance_Oval;
#ifdef USE_MPS
    chem::RestraintDistance_O::static_Kind = gctools::GCKind<chem::RestraintDistance_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintDistance_Oval,chem::RestraintDistance_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintDistance_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintDistance_O>>::allocateClass();
        chem::RestraintDistance_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintDistance_O::static_className() % (void*)(chem::RestraintDistance_O::static_allocator) );
    classchem__RestraintDistance_Oval->setCreator(chem::RestraintDistance_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintDistance_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintDistance_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RestraintFixedNonbond_Oval]"));
    core::BuiltInClass_sp classchem__RestraintFixedNonbond_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RestraintFixedNonbond_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RestraintFixedNonbond_Oval,_lisp,chem::RestraintFixedNonbond_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RestraintFixedNonbond_O>::id,chem::RestraintFixedNonbond_O::static_classSymbol());
    chem::RestraintFixedNonbond_O::___staticClass = classchem__RestraintFixedNonbond_Oval;
#ifdef USE_MPS
    chem::RestraintFixedNonbond_O::static_Kind = gctools::GCKind<chem::RestraintFixedNonbond_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RestraintFixedNonbond_Oval,chem::RestraintFixedNonbond_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RestraintFixedNonbond_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RestraintFixedNonbond_O>>::allocateClass();
        chem::RestraintFixedNonbond_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RestraintFixedNonbond_O::static_className() % (void*)(chem::RestraintFixedNonbond_O::static_allocator) );
    classchem__RestraintFixedNonbond_Oval->setCreator(chem::RestraintFixedNonbond_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RestraintFixedNonbond_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RestraintFixedNonbond_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RingCoupling_Oval]"));
    core::BuiltInClass_sp classchem__RingCoupling_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RingCoupling_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RingCoupling_Oval,_lisp,chem::RingCoupling_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RingCoupling_O>::id,chem::RingCoupling_O::static_classSymbol());
    chem::RingCoupling_O::___staticClass = classchem__RingCoupling_Oval;
#ifdef USE_MPS
    chem::RingCoupling_O::static_Kind = gctools::GCKind<chem::RingCoupling_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RingCoupling_Oval,chem::RingCoupling_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RingCoupling_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RingCoupling_O>>::allocateClass();
        chem::RingCoupling_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RingCoupling_O::static_className() % (void*)(chem::RingCoupling_O::static_allocator) );
    classchem__RingCoupling_Oval->setCreator(chem::RingCoupling_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RingCoupling_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RingCoupling_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RootMatchNode_Oval]"));
    core::BuiltInClass_sp classchem__RootMatchNode_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RootMatchNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RootMatchNode_Oval,_lisp,chem::RootMatchNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RootMatchNode_O>::id,chem::RootMatchNode_O::static_classSymbol());
    chem::RootMatchNode_O::___staticClass = classchem__RootMatchNode_Oval;
#ifdef USE_MPS
    chem::RootMatchNode_O::static_Kind = gctools::GCKind<chem::RootMatchNode_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RootMatchNode_Oval,chem::RootMatchNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RootMatchNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RootMatchNode_O>>::allocateClass();
        chem::RootMatchNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RootMatchNode_O::static_className() % (void*)(chem::RootMatchNode_O::static_allocator) );
    classchem__RootMatchNode_Oval->setCreator(chem::RootMatchNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RootMatchNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RootMatchNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SparseLargeSquareMatrix_Oval]"));
    core::BuiltInClass_sp classchem__SparseLargeSquareMatrix_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SparseLargeSquareMatrix_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SparseLargeSquareMatrix_Oval,_lisp,chem::SparseLargeSquareMatrix_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SparseLargeSquareMatrix_O>::id,chem::SparseLargeSquareMatrix_O::static_classSymbol());
    chem::SparseLargeSquareMatrix_O::___staticClass = classchem__SparseLargeSquareMatrix_Oval;
#ifdef USE_MPS
    chem::SparseLargeSquareMatrix_O::static_Kind = gctools::GCKind<chem::SparseLargeSquareMatrix_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SparseLargeSquareMatrix_Oval,chem::SparseLargeSquareMatrix_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SparseLargeSquareMatrix_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SparseLargeSquareMatrix_O>>::allocateClass();
        chem::SparseLargeSquareMatrix_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SparseLargeSquareMatrix_O::static_className() % (void*)(chem::SparseLargeSquareMatrix_O::static_allocator) );
    classchem__SparseLargeSquareMatrix_Oval->setCreator(chem::SparseLargeSquareMatrix_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SparseLargeSquareMatrix_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SparseLargeSquareMatrix_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StepReport_Oval]"));
    core::BuiltInClass_sp classchem__StepReport_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StepReport_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StepReport_Oval,_lisp,chem::StepReport_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StepReport_O>::id,chem::StepReport_O::static_classSymbol());
    chem::StepReport_O::___staticClass = classchem__StepReport_Oval;
#ifdef USE_MPS
    chem::StepReport_O::static_Kind = gctools::GCKind<chem::StepReport_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StepReport_Oval,chem::StepReport_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StepReport_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StepReport_O>>::allocateClass();
        chem::StepReport_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StepReport_O::static_className() % (void*)(chem::StepReport_O::static_allocator) );
    classchem__StepReport_Oval->setCreator(chem::StepReport_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StepReport_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StepReport_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__StereoisomerVirtualAtom_Oval]"));
    core::BuiltInClass_sp classchem__StereoisomerVirtualAtom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__StereoisomerVirtualAtom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__StereoisomerVirtualAtom_Oval,_lisp,chem::StereoisomerVirtualAtom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::StereoisomerVirtualAtom_O>::id,chem::StereoisomerVirtualAtom_O::static_classSymbol());
    chem::StereoisomerVirtualAtom_O::___staticClass = classchem__StereoisomerVirtualAtom_Oval;
#ifdef USE_MPS
    chem::StereoisomerVirtualAtom_O::static_Kind = gctools::GCKind<chem::StereoisomerVirtualAtom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__StereoisomerVirtualAtom_Oval,chem::StereoisomerVirtualAtom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::StereoisomerVirtualAtom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::StereoisomerVirtualAtom_O>>::allocateClass();
        chem::StereoisomerVirtualAtom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::StereoisomerVirtualAtom_O::static_className() % (void*)(chem::StereoisomerVirtualAtom_O::static_allocator) );
    classchem__StereoisomerVirtualAtom_Oval->setCreator(chem::StereoisomerVirtualAtom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::StereoisomerVirtualAtom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__StereoisomerVirtualAtom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Stereoisomer_Oval]"));
    core::BuiltInClass_sp classchem__Stereoisomer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Stereoisomer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Stereoisomer_Oval,_lisp,chem::Stereoisomer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Stereoisomer_O>::id,chem::Stereoisomer_O::static_classSymbol());
    chem::Stereoisomer_O::___staticClass = classchem__Stereoisomer_Oval;
#ifdef USE_MPS
    chem::Stereoisomer_O::static_Kind = gctools::GCKind<chem::Stereoisomer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Stereoisomer_Oval,chem::Stereoisomer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Stereoisomer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Stereoisomer_O>>::allocateClass();
        chem::Stereoisomer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Stereoisomer_O::static_className() % (void*)(chem::Stereoisomer_O::static_allocator) );
    classchem__Stereoisomer_Oval->setCreator(chem::Stereoisomer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Stereoisomer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Stereoisomer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SuperposableConformationCollection_Oval]"));
    core::BuiltInClass_sp classchem__SuperposableConformationCollection_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SuperposableConformationCollection_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SuperposableConformationCollection_Oval,_lisp,chem::SuperposableConformationCollection_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SuperposableConformationCollection_O>::id,chem::SuperposableConformationCollection_O::static_classSymbol());
    chem::SuperposableConformationCollection_O::___staticClass = classchem__SuperposableConformationCollection_Oval;
#ifdef USE_MPS
    chem::SuperposableConformationCollection_O::static_Kind = gctools::GCKind<chem::SuperposableConformationCollection_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SuperposableConformationCollection_Oval,chem::SuperposableConformationCollection_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SuperposableConformationCollection_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SuperposableConformationCollection_O>>::allocateClass();
        chem::SuperposableConformationCollection_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SuperposableConformationCollection_O::static_className() % (void*)(chem::SuperposableConformationCollection_O::static_allocator) );
    classchem__SuperposableConformationCollection_Oval->setCreator(chem::SuperposableConformationCollection_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SuperposableConformationCollection_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SuperposableConformationCollection_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SuperposeSelectedAtoms_Oval]"));
    core::BuiltInClass_sp classchem__SuperposeSelectedAtoms_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SuperposeSelectedAtoms_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SuperposeSelectedAtoms_Oval,_lisp,chem::SuperposeSelectedAtoms_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SuperposeSelectedAtoms_O>::id,chem::SuperposeSelectedAtoms_O::static_classSymbol());
    chem::SuperposeSelectedAtoms_O::___staticClass = classchem__SuperposeSelectedAtoms_Oval;
#ifdef USE_MPS
    chem::SuperposeSelectedAtoms_O::static_Kind = gctools::GCKind<chem::SuperposeSelectedAtoms_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SuperposeSelectedAtoms_Oval,chem::SuperposeSelectedAtoms_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SuperposeSelectedAtoms_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SuperposeSelectedAtoms_O>>::allocateClass();
        chem::SuperposeSelectedAtoms_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SuperposeSelectedAtoms_O::static_className() % (void*)(chem::SuperposeSelectedAtoms_O::static_allocator) );
    classchem__SuperposeSelectedAtoms_Oval->setCreator(chem::SuperposeSelectedAtoms_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SuperposeSelectedAtoms_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SuperposeSelectedAtoms_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ZMatrixAngleInternal_Oval]"));
    core::BuiltInClass_sp classchem__ZMatrixAngleInternal_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ZMatrixAngleInternal_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ZMatrixAngleInternal_Oval,_lisp,chem::ZMatrixAngleInternal_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ZMatrixAngleInternal_O>::id,chem::ZMatrixAngleInternal_O::static_classSymbol());
    chem::ZMatrixAngleInternal_O::___staticClass = classchem__ZMatrixAngleInternal_Oval;
#ifdef USE_MPS
    chem::ZMatrixAngleInternal_O::static_Kind = gctools::GCKind<chem::ZMatrixAngleInternal_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ZMatrixAngleInternal_Oval,chem::ZMatrixAngleInternal_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ZMatrixAngleInternal_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ZMatrixAngleInternal_O>>::allocateClass();
        chem::ZMatrixAngleInternal_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ZMatrixAngleInternal_O::static_className() % (void*)(chem::ZMatrixAngleInternal_O::static_allocator) );
    classchem__ZMatrixAngleInternal_Oval->setCreator(chem::ZMatrixAngleInternal_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ZMatrixAngleInternal_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ZMatrixAngleInternal_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ZMatrixBondInternal_Oval]"));
    core::BuiltInClass_sp classchem__ZMatrixBondInternal_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ZMatrixBondInternal_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ZMatrixBondInternal_Oval,_lisp,chem::ZMatrixBondInternal_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ZMatrixBondInternal_O>::id,chem::ZMatrixBondInternal_O::static_classSymbol());
    chem::ZMatrixBondInternal_O::___staticClass = classchem__ZMatrixBondInternal_Oval;
#ifdef USE_MPS
    chem::ZMatrixBondInternal_O::static_Kind = gctools::GCKind<chem::ZMatrixBondInternal_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ZMatrixBondInternal_Oval,chem::ZMatrixBondInternal_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ZMatrixBondInternal_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ZMatrixBondInternal_O>>::allocateClass();
        chem::ZMatrixBondInternal_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ZMatrixBondInternal_O::static_className() % (void*)(chem::ZMatrixBondInternal_O::static_allocator) );
    classchem__ZMatrixBondInternal_Oval->setCreator(chem::ZMatrixBondInternal_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ZMatrixBondInternal_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ZMatrixBondInternal_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ZMatrixDihedralInternal_Oval]"));
    core::BuiltInClass_sp classchem__ZMatrixDihedralInternal_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ZMatrixDihedralInternal_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ZMatrixDihedralInternal_Oval,_lisp,chem::ZMatrixDihedralInternal_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ZMatrixDihedralInternal_O>::id,chem::ZMatrixDihedralInternal_O::static_classSymbol());
    chem::ZMatrixDihedralInternal_O::___staticClass = classchem__ZMatrixDihedralInternal_Oval;
#ifdef USE_MPS
    chem::ZMatrixDihedralInternal_O::static_Kind = gctools::GCKind<chem::ZMatrixDihedralInternal_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ZMatrixDihedralInternal_Oval,chem::ZMatrixDihedralInternal_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ZMatrixDihedralInternal_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ZMatrixDihedralInternal_O>>::allocateClass();
        chem::ZMatrixDihedralInternal_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ZMatrixDihedralInternal_O::static_className() % (void*)(chem::ZMatrixDihedralInternal_O::static_allocator) );
    classchem__ZMatrixDihedralInternal_Oval->setCreator(chem::ZMatrixDihedralInternal_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ZMatrixDihedralInternal_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ZMatrixDihedralInternal_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AfterMatchBondTest_Oval]"));
    core::BuiltInClass_sp classchem__AfterMatchBondTest_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AfterMatchBondTest_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AfterMatchBondTest_Oval,_lisp,chem::AfterMatchBondTest_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AfterMatchBondTest_O>::id,chem::AfterMatchBondTest_O::static_classSymbol());
    chem::AfterMatchBondTest_O::___staticClass = classchem__AfterMatchBondTest_Oval;
#ifdef USE_MPS
    chem::AfterMatchBondTest_O::static_Kind = gctools::GCKind<chem::AfterMatchBondTest_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AfterMatchBondTest_Oval,chem::AfterMatchBondTest_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AfterMatchBondTest_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AfterMatchBondTest_O>>::allocateClass();
        chem::AfterMatchBondTest_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AfterMatchBondTest_O::static_className() % (void*)(chem::AfterMatchBondTest_O::static_allocator) );
    classchem__AfterMatchBondTest_Oval->setCreator(chem::AfterMatchBondTest_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AfterMatchBondTest_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AfterMatchBondTest_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AntechamberBondTest_Oval]"));
    core::BuiltInClass_sp classchem__AntechamberBondTest_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AntechamberBondTest_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AntechamberBondTest_Oval,_lisp,chem::AntechamberBondTest_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AntechamberBondTest_O>::id,chem::AntechamberBondTest_O::static_classSymbol());
    chem::AntechamberBondTest_O::___staticClass = classchem__AntechamberBondTest_Oval;
#ifdef USE_MPS
    chem::AntechamberBondTest_O::static_Kind = gctools::GCKind<chem::AntechamberBondTest_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AntechamberBondTest_Oval,chem::AntechamberBondTest_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AntechamberBondTest_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AntechamberBondTest_O>>::allocateClass();
        chem::AntechamberBondTest_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AntechamberBondTest_O::static_className() % (void*)(chem::AntechamberBondTest_O::static_allocator) );
    classchem__AntechamberBondTest_Oval->setCreator(chem::AntechamberBondTest_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AntechamberBondTest_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AntechamberBondTest_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomOrBondMatchNode_Oval]"));
    core::BuiltInClass_sp classchem__AtomOrBondMatchNode_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomOrBondMatchNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomOrBondMatchNode_Oval,_lisp,chem::AtomOrBondMatchNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomOrBondMatchNode_O>::id,chem::AtomOrBondMatchNode_O::static_classSymbol());
    chem::AtomOrBondMatchNode_O::___staticClass = classchem__AtomOrBondMatchNode_Oval;
#ifdef USE_MPS
    chem::AtomOrBondMatchNode_O::static_Kind = gctools::GCKind<chem::AtomOrBondMatchNode_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomOrBondMatchNode_Oval,chem::AtomOrBondMatchNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomOrBondMatchNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomOrBondMatchNode_O>>::allocateClass();
        chem::AtomOrBondMatchNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomOrBondMatchNode_O::static_className() % (void*)(chem::AtomOrBondMatchNode_O::static_allocator) );
    classchem__AtomOrBondMatchNode_Oval->setCreator(chem::AtomOrBondMatchNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomOrBondMatchNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomOrBondMatchNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__BondTest_Oval]"));
    core::BuiltInClass_sp classchem__BondTest_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__BondTest_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__BondTest_Oval,_lisp,chem::BondTest_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::BondTest_O>::id,chem::BondTest_O::static_classSymbol());
    chem::BondTest_O::___staticClass = classchem__BondTest_Oval;
#ifdef USE_MPS
    chem::BondTest_O::static_Kind = gctools::GCKind<chem::BondTest_O>::Kind;
#endif
    core::af_setf_findClass(classchem__BondTest_Oval,chem::BondTest_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::BondTest_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::BondTest_O>>::allocateClass();
        chem::BondTest_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::BondTest_O::static_className() % (void*)(chem::BondTest_O::static_allocator) );
    classchem__BondTest_Oval->setCreator(chem::BondTest_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::BondTest_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__BondTest_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Branch_Oval]"));
    core::BuiltInClass_sp classchem__Branch_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Branch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Branch_Oval,_lisp,chem::Branch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Branch_O>::id,chem::Branch_O::static_classSymbol());
    chem::Branch_O::___staticClass = classchem__Branch_Oval;
#ifdef USE_MPS
    chem::Branch_O::static_Kind = gctools::GCKind<chem::Branch_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Branch_Oval,chem::Branch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Branch_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Branch_O>>::allocateClass();
        chem::Branch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Branch_O::static_className() % (void*)(chem::Branch_O::static_allocator) );
    classchem__Branch_Oval->setCreator(chem::Branch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Branch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Branch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Chain_Oval]"));
    core::BuiltInClass_sp classchem__Chain_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Chain_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Chain_Oval,_lisp,chem::Chain_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Chain_O>::id,chem::Chain_O::static_classSymbol());
    chem::Chain_O::___staticClass = classchem__Chain_Oval;
#ifdef USE_MPS
    chem::Chain_O::static_Kind = gctools::GCKind<chem::Chain_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Chain_Oval,chem::Chain_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Chain_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Chain_O>>::allocateClass();
        chem::Chain_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Chain_O::static_className() % (void*)(chem::Chain_O::static_allocator) );
    classchem__Chain_Oval->setCreator(chem::Chain_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Chain_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Chain_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EntityNameSetWithCap_Oval]"));
    core::BuiltInClass_sp classchem__EntityNameSetWithCap_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EntityNameSetWithCap_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EntityNameSetWithCap_Oval,_lisp,chem::EntityNameSetWithCap_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EntityNameSetWithCap_O>::id,chem::EntityNameSetWithCap_O::static_classSymbol());
    chem::EntityNameSetWithCap_O::___staticClass = classchem__EntityNameSetWithCap_Oval;
#ifdef USE_MPS
    chem::EntityNameSetWithCap_O::static_Kind = gctools::GCKind<chem::EntityNameSetWithCap_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EntityNameSetWithCap_Oval,chem::EntityNameSetWithCap_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EntityNameSetWithCap_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EntityNameSetWithCap_O>>::allocateClass();
        chem::EntityNameSetWithCap_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EntityNameSetWithCap_O::static_className() % (void*)(chem::EntityNameSetWithCap_O::static_allocator) );
    classchem__EntityNameSetWithCap_Oval->setCreator(chem::EntityNameSetWithCap_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EntityNameSetWithCap_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EntityNameSetWithCap_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__EntityNameSet_Oval]"));
    core::BuiltInClass_sp classchem__EntityNameSet_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__EntityNameSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__EntityNameSet_Oval,_lisp,chem::EntityNameSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::EntityNameSet_O>::id,chem::EntityNameSet_O::static_classSymbol());
    chem::EntityNameSet_O::___staticClass = classchem__EntityNameSet_Oval;
#ifdef USE_MPS
    chem::EntityNameSet_O::static_Kind = gctools::GCKind<chem::EntityNameSet_O>::Kind;
#endif
    core::af_setf_findClass(classchem__EntityNameSet_Oval,chem::EntityNameSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::EntityNameSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::EntityNameSet_O>>::allocateClass();
        chem::EntityNameSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::EntityNameSet_O::static_className() % (void*)(chem::EntityNameSet_O::static_allocator) );
    classchem__EntityNameSet_Oval->setCreator(chem::EntityNameSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::EntityNameSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__EntityNameSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__InPlug_Oval]"));
    core::BuiltInClass_sp classchem__InPlug_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__InPlug_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__InPlug_Oval,_lisp,chem::InPlug_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::InPlug_O>::id,chem::InPlug_O::static_classSymbol());
    chem::InPlug_O::___staticClass = classchem__InPlug_Oval;
#ifdef USE_MPS
    chem::InPlug_O::static_Kind = gctools::GCKind<chem::InPlug_O>::Kind;
#endif
    core::af_setf_findClass(classchem__InPlug_Oval,chem::InPlug_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::InPlug_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::InPlug_O>>::allocateClass();
        chem::InPlug_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::InPlug_O::static_className() % (void*)(chem::InPlug_O::static_allocator) );
    classchem__InPlug_Oval->setCreator(chem::InPlug_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::InPlug_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__InPlug_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__IterateAtoms_Oval]"));
    core::BuiltInClass_sp classchem__IterateAtoms_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__IterateAtoms_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__IterateAtoms_Oval,_lisp,chem::IterateAtoms_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::IterateAtoms_O>::id,chem::IterateAtoms_O::static_classSymbol());
    chem::IterateAtoms_O::___staticClass = classchem__IterateAtoms_Oval;
#ifdef USE_MPS
    chem::IterateAtoms_O::static_Kind = gctools::GCKind<chem::IterateAtoms_O>::Kind;
#endif
    core::af_setf_findClass(classchem__IterateAtoms_Oval,chem::IterateAtoms_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::IterateAtoms_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::IterateAtoms_O>>::allocateClass();
        chem::IterateAtoms_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::IterateAtoms_O::static_className() % (void*)(chem::IterateAtoms_O::static_allocator) );
    classchem__IterateAtoms_Oval->setCreator(chem::IterateAtoms_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::IterateAtoms_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__IterateAtoms_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__IterateBonds_Oval]"));
    core::BuiltInClass_sp classchem__IterateBonds_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__IterateBonds_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__IterateBonds_Oval,_lisp,chem::IterateBonds_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::IterateBonds_O>::id,chem::IterateBonds_O::static_classSymbol());
    chem::IterateBonds_O::___staticClass = classchem__IterateBonds_Oval;
#ifdef USE_MPS
    chem::IterateBonds_O::static_Kind = gctools::GCKind<chem::IterateBonds_O>::Kind;
#endif
    core::af_setf_findClass(classchem__IterateBonds_Oval,chem::IterateBonds_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::IterateBonds_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::IterateBonds_O>>::allocateClass();
        chem::IterateBonds_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::IterateBonds_O::static_className() % (void*)(chem::IterateBonds_O::static_allocator) );
    classchem__IterateBonds_Oval->setCreator(chem::IterateBonds_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::IterateBonds_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__IterateBonds_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__IterateResidues_Oval]"));
    core::BuiltInClass_sp classchem__IterateResidues_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__IterateResidues_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__IterateResidues_Oval,_lisp,chem::IterateResidues_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::IterateResidues_O>::id,chem::IterateResidues_O::static_classSymbol());
    chem::IterateResidues_O::___staticClass = classchem__IterateResidues_Oval;
#ifdef USE_MPS
    chem::IterateResidues_O::static_Kind = gctools::GCKind<chem::IterateResidues_O>::Kind;
#endif
    core::af_setf_findClass(classchem__IterateResidues_Oval,chem::IterateResidues_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::IterateResidues_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::IterateResidues_O>>::allocateClass();
        chem::IterateResidues_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::IterateResidues_O::static_className() % (void*)(chem::IterateResidues_O::static_allocator) );
    classchem__IterateResidues_Oval->setCreator(chem::IterateResidues_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::IterateResidues_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__IterateResidues_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MultiMonomer_Oval]"));
    core::BuiltInClass_sp classchem__MultiMonomer_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MultiMonomer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MultiMonomer_Oval,_lisp,chem::MultiMonomer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MultiMonomer_O>::id,chem::MultiMonomer_O::static_classSymbol());
    chem::MultiMonomer_O::___staticClass = classchem__MultiMonomer_Oval;
#ifdef USE_MPS
    chem::MultiMonomer_O::static_Kind = gctools::GCKind<chem::MultiMonomer_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MultiMonomer_Oval,chem::MultiMonomer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MultiMonomer_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MultiMonomer_O>>::allocateClass();
        chem::MultiMonomer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MultiMonomer_O::static_className() % (void*)(chem::MultiMonomer_O::static_allocator) );
    classchem__MultiMonomer_Oval->setCreator(chem::MultiMonomer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MultiMonomer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MultiMonomer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__OutPlug_Oval]"));
    core::BuiltInClass_sp classchem__OutPlug_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__OutPlug_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__OutPlug_Oval,_lisp,chem::OutPlug_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::OutPlug_O>::id,chem::OutPlug_O::static_classSymbol());
    chem::OutPlug_O::___staticClass = classchem__OutPlug_Oval;
#ifdef USE_MPS
    chem::OutPlug_O::static_Kind = gctools::GCKind<chem::OutPlug_O>::Kind;
#endif
    core::af_setf_findClass(classchem__OutPlug_Oval,chem::OutPlug_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::OutPlug_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::OutPlug_O>>::allocateClass();
        chem::OutPlug_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::OutPlug_O::static_className() % (void*)(chem::OutPlug_O::static_allocator) );
    classchem__OutPlug_Oval->setCreator(chem::OutPlug_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::OutPlug_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__OutPlug_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__VirtualAtom_Oval]"));
    core::BuiltInClass_sp classchem__VirtualAtom_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__VirtualAtom_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__VirtualAtom_Oval,_lisp,chem::VirtualAtom_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::VirtualAtom_O>::id,chem::VirtualAtom_O::static_classSymbol());
    chem::VirtualAtom_O::___staticClass = classchem__VirtualAtom_Oval;
#ifdef USE_MPS
    chem::VirtualAtom_O::static_Kind = gctools::GCKind<chem::VirtualAtom_O>::Kind;
#endif
    core::af_setf_findClass(classchem__VirtualAtom_Oval,chem::VirtualAtom_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::VirtualAtom_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::VirtualAtom_O>>::allocateClass();
        chem::VirtualAtom_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::VirtualAtom_O::static_className() % (void*)(chem::VirtualAtom_O::static_allocator) );
    classchem__VirtualAtom_Oval->setCreator(chem::VirtualAtom_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::VirtualAtom_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__VirtualAtom_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AntechamberFocusAtomMatch_Oval]"));
    core::BuiltInClass_sp classchem__AntechamberFocusAtomMatch_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AntechamberFocusAtomMatch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AntechamberFocusAtomMatch_Oval,_lisp,chem::AntechamberFocusAtomMatch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AntechamberFocusAtomMatch_O>::id,chem::AntechamberFocusAtomMatch_O::static_classSymbol());
    chem::AntechamberFocusAtomMatch_O::___staticClass = classchem__AntechamberFocusAtomMatch_Oval;
#ifdef USE_MPS
    chem::AntechamberFocusAtomMatch_O::static_Kind = gctools::GCKind<chem::AntechamberFocusAtomMatch_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AntechamberFocusAtomMatch_Oval,chem::AntechamberFocusAtomMatch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AntechamberFocusAtomMatch_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AntechamberFocusAtomMatch_O>>::allocateClass();
        chem::AntechamberFocusAtomMatch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AntechamberFocusAtomMatch_O::static_className() % (void*)(chem::AntechamberFocusAtomMatch_O::static_allocator) );
    classchem__AntechamberFocusAtomMatch_Oval->setCreator(chem::AntechamberFocusAtomMatch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AntechamberFocusAtomMatch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AntechamberFocusAtomMatch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AtomTest_Oval]"));
    core::BuiltInClass_sp classchem__AtomTest_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AtomTest_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AtomTest_Oval,_lisp,chem::AtomTest_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AtomTest_O>::id,chem::AtomTest_O::static_classSymbol());
    chem::AtomTest_O::___staticClass = classchem__AtomTest_Oval;
#ifdef USE_MPS
    chem::AtomTest_O::static_Kind = gctools::GCKind<chem::AtomTest_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AtomTest_Oval,chem::AtomTest_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AtomTest_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AtomTest_O>>::allocateClass();
        chem::AtomTest_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AtomTest_O::static_className() % (void*)(chem::AtomTest_O::static_allocator) );
    classchem__AtomTest_Oval->setCreator(chem::AtomTest_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AtomTest_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AtomTest_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Logical_Oval]"));
    core::BuiltInClass_sp classchem__Logical_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Logical_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Logical_Oval,_lisp,chem::Logical_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Logical_O>::id,chem::Logical_O::static_classSymbol());
    chem::Logical_O::___staticClass = classchem__Logical_Oval;
#ifdef USE_MPS
    chem::Logical_O::static_Kind = gctools::GCKind<chem::Logical_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Logical_Oval,chem::Logical_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Logical_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Logical_O>>::allocateClass();
        chem::Logical_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Logical_O::static_className() % (void*)(chem::Logical_O::static_allocator) );
    classchem__Logical_Oval->setCreator(chem::Logical_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Logical_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Logical_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Mate_Oval]"));
    core::BuiltInClass_sp classchem__Mate_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Mate_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Mate_Oval,_lisp,chem::Mate_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Mate_O>::id,chem::Mate_O::static_classSymbol());
    chem::Mate_O::___staticClass = classchem__Mate_Oval;
#ifdef USE_MPS
    chem::Mate_O::static_Kind = gctools::GCKind<chem::Mate_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Mate_Oval,chem::Mate_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Mate_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Mate_O>>::allocateClass();
        chem::Mate_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Mate_O::static_className() % (void*)(chem::Mate_O::static_allocator) );
    classchem__Mate_Oval->setCreator(chem::Mate_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Mate_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Mate_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__MonomerPack_Oval]"));
    core::BuiltInClass_sp classchem__MonomerPack_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__MonomerPack_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__MonomerPack_Oval,_lisp,chem::MonomerPack_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::MonomerPack_O>::id,chem::MonomerPack_O::static_classSymbol());
    chem::MonomerPack_O::___staticClass = classchem__MonomerPack_Oval;
#ifdef USE_MPS
    chem::MonomerPack_O::static_Kind = gctools::GCKind<chem::MonomerPack_O>::Kind;
#endif
    core::af_setf_findClass(classchem__MonomerPack_Oval,chem::MonomerPack_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::MonomerPack_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::MonomerPack_O>>::allocateClass();
        chem::MonomerPack_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::MonomerPack_O::static_className() % (void*)(chem::MonomerPack_O::static_allocator) );
    classchem__MonomerPack_Oval->setCreator(chem::MonomerPack_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::MonomerPack_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__MonomerPack_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RepresentedEntityNameSet_Oval]"));
    core::BuiltInClass_sp classchem__RepresentedEntityNameSet_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RepresentedEntityNameSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RepresentedEntityNameSet_Oval,_lisp,chem::RepresentedEntityNameSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RepresentedEntityNameSet_O>::id,chem::RepresentedEntityNameSet_O::static_classSymbol());
    chem::RepresentedEntityNameSet_O::___staticClass = classchem__RepresentedEntityNameSet_Oval;
#ifdef USE_MPS
    chem::RepresentedEntityNameSet_O::static_Kind = gctools::GCKind<chem::RepresentedEntityNameSet_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RepresentedEntityNameSet_Oval,chem::RepresentedEntityNameSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RepresentedEntityNameSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RepresentedEntityNameSet_O>>::allocateClass();
        chem::RepresentedEntityNameSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RepresentedEntityNameSet_O::static_className() % (void*)(chem::RepresentedEntityNameSet_O::static_allocator) );
    classchem__RepresentedEntityNameSet_Oval->setCreator(chem::RepresentedEntityNameSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RepresentedEntityNameSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RepresentedEntityNameSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__ResidueTest_Oval]"));
    core::BuiltInClass_sp classchem__ResidueTest_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__ResidueTest_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__ResidueTest_Oval,_lisp,chem::ResidueTest_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::ResidueTest_O>::id,chem::ResidueTest_O::static_classSymbol());
    chem::ResidueTest_O::___staticClass = classchem__ResidueTest_Oval;
#ifdef USE_MPS
    chem::ResidueTest_O::static_Kind = gctools::GCKind<chem::ResidueTest_O>::Kind;
#endif
    core::af_setf_findClass(classchem__ResidueTest_Oval,chem::ResidueTest_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::ResidueTest_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::ResidueTest_O>>::allocateClass();
        chem::ResidueTest_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::ResidueTest_O::static_className() % (void*)(chem::ResidueTest_O::static_allocator) );
    classchem__ResidueTest_Oval->setCreator(chem::ResidueTest_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::ResidueTest_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__ResidueTest_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RingClosingMate_Oval]"));
    core::BuiltInClass_sp classchem__RingClosingMate_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RingClosingMate_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RingClosingMate_Oval,_lisp,chem::RingClosingMate_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RingClosingMate_O>::id,chem::RingClosingMate_O::static_classSymbol());
    chem::RingClosingMate_O::___staticClass = classchem__RingClosingMate_Oval;
#ifdef USE_MPS
    chem::RingClosingMate_O::static_Kind = gctools::GCKind<chem::RingClosingMate_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RingClosingMate_Oval,chem::RingClosingMate_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RingClosingMate_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RingClosingMate_O>>::allocateClass();
        chem::RingClosingMate_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RingClosingMate_O::static_className() % (void*)(chem::RingClosingMate_O::static_allocator) );
    classchem__RingClosingMate_Oval->setCreator(chem::RingClosingMate_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RingClosingMate_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RingClosingMate_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RingClosingPlug_Oval]"));
    core::BuiltInClass_sp classchem__RingClosingPlug_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RingClosingPlug_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RingClosingPlug_Oval,_lisp,chem::RingClosingPlug_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RingClosingPlug_O>::id,chem::RingClosingPlug_O::static_classSymbol());
    chem::RingClosingPlug_O::___staticClass = classchem__RingClosingPlug_Oval;
#ifdef USE_MPS
    chem::RingClosingPlug_O::static_Kind = gctools::GCKind<chem::RingClosingPlug_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RingClosingPlug_Oval,chem::RingClosingPlug_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RingClosingPlug_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RingClosingPlug_O>>::allocateClass();
        chem::RingClosingPlug_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RingClosingPlug_O::static_className() % (void*)(chem::RingClosingPlug_O::static_allocator) );
    classchem__RingClosingPlug_Oval->setCreator(chem::RingClosingPlug_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RingClosingPlug_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RingClosingPlug_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__RingTest_Oval]"));
    core::BuiltInClass_sp classchem__RingTest_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__RingTest_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__RingTest_Oval,_lisp,chem::RingTest_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::RingTest_O>::id,chem::RingTest_O::static_classSymbol());
    chem::RingTest_O::___staticClass = classchem__RingTest_Oval;
#ifdef USE_MPS
    chem::RingTest_O::static_Kind = gctools::GCKind<chem::RingTest_O>::Kind;
#endif
    core::af_setf_findClass(classchem__RingTest_Oval,chem::RingTest_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::RingTest_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::RingTest_O>>::allocateClass();
        chem::RingTest_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::RingTest_O::static_className() % (void*)(chem::RingTest_O::static_allocator) );
    classchem__RingTest_Oval->setCreator(chem::RingTest_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::RingTest_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__RingTest_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__Root_Oval]"));
    core::BuiltInClass_sp classchem__Root_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__Root_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__Root_Oval,_lisp,chem::Root_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::Root_O>::id,chem::Root_O::static_classSymbol());
    chem::Root_O::___staticClass = classchem__Root_Oval;
#ifdef USE_MPS
    chem::Root_O::static_Kind = gctools::GCKind<chem::Root_O>::Kind;
#endif
    core::af_setf_findClass(classchem__Root_Oval,chem::Root_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::Root_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::Root_O>>::allocateClass();
        chem::Root_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::Root_O::static_className() % (void*)(chem::Root_O::static_allocator) );
    classchem__Root_Oval->setCreator(chem::Root_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::Root_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__Root_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__TagSet_Oval]"));
    core::BuiltInClass_sp classchem__TagSet_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__TagSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__TagSet_Oval,_lisp,chem::TagSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::TagSet_O>::id,chem::TagSet_O::static_classSymbol());
    chem::TagSet_O::___staticClass = classchem__TagSet_Oval;
#ifdef USE_MPS
    chem::TagSet_O::static_Kind = gctools::GCKind<chem::TagSet_O>::Kind;
#endif
    core::af_setf_findClass(classchem__TagSet_Oval,chem::TagSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::TagSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::TagSet_O>>::allocateClass();
        chem::TagSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::TagSet_O::static_className() % (void*)(chem::TagSet_O::static_allocator) );
    classchem__TagSet_Oval->setCreator(chem::TagSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::TagSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__TagSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__AntechamberRoot_Oval]"));
    core::BuiltInClass_sp classchem__AntechamberRoot_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__AntechamberRoot_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__AntechamberRoot_Oval,_lisp,chem::AntechamberRoot_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::AntechamberRoot_O>::id,chem::AntechamberRoot_O::static_classSymbol());
    chem::AntechamberRoot_O::___staticClass = classchem__AntechamberRoot_Oval;
#ifdef USE_MPS
    chem::AntechamberRoot_O::static_Kind = gctools::GCKind<chem::AntechamberRoot_O>::Kind;
#endif
    core::af_setf_findClass(classchem__AntechamberRoot_Oval,chem::AntechamberRoot_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::AntechamberRoot_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::AntechamberRoot_O>>::allocateClass();
        chem::AntechamberRoot_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::AntechamberRoot_O::static_className() % (void*)(chem::AntechamberRoot_O::static_allocator) );
    classchem__AntechamberRoot_Oval->setCreator(chem::AntechamberRoot_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::AntechamberRoot_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__AntechamberRoot_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

    LOG(BF("Creating class[classchem__SmartsRoot_Oval]"));
    core::BuiltInClass_sp classchem__SmartsRoot_Oval = core::BuiltInClass_O::createUncollectable();
    classchem__SmartsRoot_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classchem__SmartsRoot_Oval,_lisp,chem::SmartsRoot_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<chem::SmartsRoot_O>::id,chem::SmartsRoot_O::static_classSymbol());
    chem::SmartsRoot_O::___staticClass = classchem__SmartsRoot_Oval;
#ifdef USE_MPS
    chem::SmartsRoot_O::static_Kind = gctools::GCKind<chem::SmartsRoot_O>::Kind;
#endif
    core::af_setf_findClass(classchem__SmartsRoot_Oval,chem::SmartsRoot_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        core::LispObjectCreator<chem::SmartsRoot_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<chem::SmartsRoot_O>>::allocateClass();
        chem::SmartsRoot_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% chem::SmartsRoot_O::static_className() % (void*)(chem::SmartsRoot_O::static_allocator) );
    classchem__SmartsRoot_Oval->setCreator(chem::SmartsRoot_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % chem::SmartsRoot_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classchem__SmartsRoot_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif // CREATE_CLASS
#undef CREATE_CLASS
#ifdef DUMP_INFO_CLASS // {
// Depends on nothing

    LOG(BF("---    dump_info   --- className: chem::AGEdge_O @ %X") % classchem__AGEdge_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AGEdge_O::static_className() % chem::AGEdge_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AGVertex_O @ %X") % classchem__AGVertex_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AGVertex_O::static_className() % chem::AGVertex_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AbstractLargeSquareMatrix_O @ %X") % classchem__AbstractLargeSquareMatrix_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AbstractLargeSquareMatrix_O::static_className() % chem::AbstractLargeSquareMatrix_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Alias_O @ %X") % classchem__Alias_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Alias_O::static_className() % chem::Alias_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Angle_O @ %X") % classchem__Angle_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Angle_O::static_className() % chem::Angle_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomIdToAtomMap_O @ %X") % classchem__AtomIdToAtomMap_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomIdToAtomMap_O::static_className() % chem::AtomIdToAtomMap_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomId_O @ %X") % classchem__AtomId_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomId_O::static_className() % chem::AtomId_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomIndexer_O @ %X") % classchem__AtomIndexer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomIndexer_O::static_className() % chem::AtomIndexer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomReferenceBase_O @ %X") % classchem__AtomReferenceBase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomReferenceBase_O::static_className() % chem::AtomReferenceBase_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomTable_O @ %X") % classchem__AtomTable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomTable_O::static_className() % chem::AtomTable_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::BondList_O @ %X") % classchem__BondList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::BondList_O::static_className() % chem::BondList_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Bond_O @ %X") % classchem__Bond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Bond_O::static_className() % chem::Bond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CDBond_O @ %X") % classchem__CDBond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CDBond_O::static_className() % chem::CDBond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CDFragment_O @ %X") % classchem__CDFragment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CDFragment_O::static_className() % chem::CDFragment_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CDNode_O @ %X") % classchem__CDNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CDNode_O::static_className() % chem::CDNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CDText_O @ %X") % classchem__CDText_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CDText_O::static_className() % chem::CDText_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CalculatePosition_O @ %X") % classchem__CalculatePosition_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CalculatePosition_O::static_className() % chem::CalculatePosition_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CandoDatabase_O @ %X") % classchem__CandoDatabase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CandoDatabase_O::static_className() % chem::CandoDatabase_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ChemDraw_O @ %X") % classchem__ChemDraw_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ChemDraw_O::static_className() % chem::ChemDraw_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ChemInfoMatch_O @ %X") % classchem__ChemInfoMatch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ChemInfoMatch_O::static_className() % chem::ChemInfoMatch_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ChemInfoNode_O @ %X") % classchem__ChemInfoNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ChemInfoNode_O::static_className() % chem::ChemInfoNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ChemInfo_O @ %X") % classchem__ChemInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ChemInfo_O::static_className() % chem::ChemInfo_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CipPrioritizer_O @ %X") % classchem__CipPrioritizer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CipPrioritizer_O::static_className() % chem::CipPrioritizer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Command_O @ %X") % classchem__Command_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Command_O::static_className() % chem::Command_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ComplexRestraint_O @ %X") % classchem__ComplexRestraint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ComplexRestraint_O::static_className() % chem::ComplexRestraint_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConformationCollectionEntry_O @ %X") % classchem__ConformationCollectionEntry_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConformationCollectionEntry_O::static_className() % chem::ConformationCollectionEntry_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConformationCollection_O @ %X") % classchem__ConformationCollection_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConformationCollection_O::static_className() % chem::ConformationCollection_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConformationExplorerEntryStage_O @ %X") % classchem__ConformationExplorerEntryStage_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConformationExplorerEntryStage_O::static_className() % chem::ConformationExplorerEntryStage_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConformationExplorerEntry_O @ %X") % classchem__ConformationExplorerEntry_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConformationExplorerEntry_O::static_className() % chem::ConformationExplorerEntry_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConformationExplorerMatch_O @ %X") % classchem__ConformationExplorerMatch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConformationExplorerMatch_O::static_className() % chem::ConformationExplorerMatch_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConformationExplorer_O @ %X") % classchem__ConformationExplorer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConformationExplorer_O::static_className() % chem::ConformationExplorer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConstitutionAtom_O @ %X") % classchem__ConstitutionAtom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConstitutionAtom_O::static_className() % chem::ConstitutionAtom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConstitutionAtoms_O @ %X") % classchem__ConstitutionAtoms_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConstitutionAtoms_O::static_className() % chem::ConstitutionAtoms_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConstitutionBond_O @ %X") % classchem__ConstitutionBond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConstitutionBond_O::static_className() % chem::ConstitutionBond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CoordinateSystem_O @ %X") % classchem__CoordinateSystem_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CoordinateSystem_O::static_className() % chem::CoordinateSystem_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Coupling_O @ %X") % classchem__Coupling_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Coupling_O::static_className() % chem::Coupling_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyComponent_O @ %X") % classchem__EnergyComponent_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyComponent_O::static_className() % chem::EnergyComponent_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyFunction_O @ %X") % classchem__EnergyFunction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyFunction_O::static_className() % chem::EnergyFunction_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Entity_O @ %X") % classchem__Entity_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Entity_O::static_className() % chem::Entity_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFBaseDb_O @ %X") % classchem__FFBaseDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFBaseDb_O::static_className() % chem::FFBaseDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFNonbondCrossTermTable_O @ %X") % classchem__FFNonbondCrossTermTable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFNonbondCrossTermTable_O::static_className() % chem::FFNonbondCrossTermTable_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFParameter_O @ %X") % classchem__FFParameter_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFParameter_O::static_className() % chem::FFParameter_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ForceField_O @ %X") % classchem__ForceField_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ForceField_O::static_className() % chem::ForceField_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FrameRecognizer_O @ %X") % classchem__FrameRecognizer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FrameRecognizer_O::static_className() % chem::FrameRecognizer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ImproperTorsion_O @ %X") % classchem__ImproperTorsion_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ImproperTorsion_O::static_className() % chem::ImproperTorsion_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::InfoDb_O @ %X") % classchem__InfoDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::InfoDb_O::static_className() % chem::InfoDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MacroModelFile_O @ %X") % classchem__MacroModelFile_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MacroModelFile_O::static_className() % chem::MacroModelFile_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MapOfMonomerNamesToAtomIndexers_O @ %X") % classchem__MapOfMonomerNamesToAtomIndexers_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MapOfMonomerNamesToAtomIndexers_O::static_className() % chem::MapOfMonomerNamesToAtomIndexers_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Matter_O @ %X") % classchem__Matter_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Matter_O::static_className() % chem::Matter_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MinimizerLog_O @ %X") % classchem__MinimizerLog_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MinimizerLog_O::static_className() % chem::MinimizerLog_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Minimizer_O @ %X") % classchem__Minimizer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Minimizer_O::static_className() % chem::Minimizer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MonomerContext_O @ %X") % classchem__MonomerContext_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MonomerContext_O::static_className() % chem::MonomerContext_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MonomerCoordinates_O @ %X") % classchem__MonomerCoordinates_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MonomerCoordinates_O::static_className() % chem::MonomerCoordinates_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::NVector_O @ %X") % classchem__NVector_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::NVector_O::static_className() % chem::NVector_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::NumericalFunction_O @ %X") % classchem__NumericalFunction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::NumericalFunction_O::static_className() % chem::NumericalFunction_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::OligomerPart_Base_O @ %X") % classchem__OligomerPart_Base_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::OligomerPart_Base_O::static_className() % chem::OligomerPart_Base_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Oligomer_O @ %X") % classchem__Oligomer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Oligomer_O::static_className() % chem::Oligomer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::OneMonomer_O @ %X") % classchem__OneMonomer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::OneMonomer_O::static_className() % chem::OneMonomer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::OneTypeRule_O @ %X") % classchem__OneTypeRule_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::OneTypeRule_O::static_className() % chem::OneTypeRule_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::PathMessage_O @ %X") % classchem__PathMessage_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::PathMessage_O::static_className() % chem::PathMessage_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::PdbMonomerConnectivity_O @ %X") % classchem__PdbMonomerConnectivity_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::PdbMonomerConnectivity_O::static_className() % chem::PdbMonomerConnectivity_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::PdbMonomerDatabase_O @ %X") % classchem__PdbMonomerDatabase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::PdbMonomerDatabase_O::static_className() % chem::PdbMonomerDatabase_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::PdbReader_O @ %X") % classchem__PdbReader_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::PdbReader_O::static_className() % chem::PdbReader_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::PdbWriter_O @ %X") % classchem__PdbWriter_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::PdbWriter_O::static_className() % chem::PdbWriter_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Plug_O @ %X") % classchem__Plug_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Plug_O::static_className() % chem::Plug_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ProperTorsion_O @ %X") % classchem__ProperTorsion_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ProperTorsion_O::static_className() % chem::ProperTorsion_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ReadAmberParameters_O @ %X") % classchem__ReadAmberParameters_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ReadAmberParameters_O::static_className() % chem::ReadAmberParameters_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ReportBase_O @ %X") % classchem__ReportBase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ReportBase_O::static_className() % chem::ReportBase_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ResidueList_O @ %X") % classchem__ResidueList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ResidueList_O::static_className() % chem::ResidueList_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintList_O @ %X") % classchem__RestraintList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintList_O::static_className() % chem::RestraintList_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Restraint_O @ %X") % classchem__Restraint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Restraint_O::static_className() % chem::Restraint_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RingFinder_O @ %X") % classchem__RingFinder_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RingFinder_O::static_className() % chem::RingFinder_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Rotamer_O @ %X") % classchem__Rotamer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Rotamer_O::static_className() % chem::Rotamer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SearchStatistics_O @ %X") % classchem__SearchStatistics_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SearchStatistics_O::static_className() % chem::SearchStatistics_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SpanningLoop_O @ %X") % classchem__SpanningLoop_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SpanningLoop_O::static_className() % chem::SpanningLoop_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SpecificContextSet_O @ %X") % classchem__SpecificContextSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SpecificContextSet_O::static_className() % chem::SpecificContextSet_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SpecificContext_O @ %X") % classchem__SpecificContext_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SpecificContext_O::static_className() % chem::SpecificContext_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StatusTracker_O @ %X") % classchem__StatusTracker_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StatusTracker_O::static_className() % chem::StatusTracker_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StereoConfiguration_O @ %X") % classchem__StereoConfiguration_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StereoConfiguration_O::static_className() % chem::StereoConfiguration_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StereoInformation_O @ %X") % classchem__StereoInformation_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StereoInformation_O::static_className() % chem::StereoInformation_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StereoisomerAtom_O @ %X") % classchem__StereoisomerAtom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StereoisomerAtom_O::static_className() % chem::StereoisomerAtom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StereoisomerAtoms_O @ %X") % classchem__StereoisomerAtoms_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StereoisomerAtoms_O::static_className() % chem::StereoisomerAtoms_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StructureComparer_O @ %X") % classchem__StructureComparer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StructureComparer_O::static_className() % chem::StructureComparer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Structure_Old_ListEntry_O @ %X") % classchem__Structure_Old_ListEntry_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Structure_Old_ListEntry_O::static_className() % chem::Structure_Old_ListEntry_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Structure_Old_List_O @ %X") % classchem__Structure_Old_List_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Structure_Old_List_O::static_className() % chem::Structure_Old_List_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SuperposeEngine_O @ %X") % classchem__SuperposeEngine_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SuperposeEngine_O::static_className() % chem::SuperposeEngine_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Topology_O @ %X") % classchem__Topology_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Topology_O::static_className() % chem::Topology_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::TrajectoryFrame_O @ %X") % classchem__TrajectoryFrame_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::TrajectoryFrame_O::static_className() % chem::TrajectoryFrame_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Trajectory_O @ %X") % classchem__Trajectory_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Trajectory_O::static_className() % chem::Trajectory_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::TwisterDriver_O @ %X") % classchem__TwisterDriver_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::TwisterDriver_O::static_className() % chem::TwisterDriver_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Twister_O @ %X") % classchem__Twister_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Twister_O::static_className() % chem::Twister_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::TypeAssignmentRules_O @ %X") % classchem__TypeAssignmentRules_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::TypeAssignmentRules_O::static_className() % chem::TypeAssignmentRules_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::VirtualSphere_O @ %X") % classchem__VirtualSphere_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::VirtualSphere_O::static_className() % chem::VirtualSphere_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::WildElementDict_O @ %X") % classchem__WildElementDict_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::WildElementDict_O::static_className() % chem::WildElementDict_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ZMatrixEntry_O @ %X") % classchem__ZMatrixEntry_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ZMatrixEntry_O::static_className() % chem::ZMatrixEntry_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ZMatrixInternal_O @ %X") % classchem__ZMatrixInternal_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ZMatrixInternal_O::static_className() % chem::ZMatrixInternal_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ZMatrix_O @ %X") % classchem__ZMatrix_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ZMatrix_O::static_className() % chem::ZMatrix_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Aggregate_O @ %X") % classchem__Aggregate_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Aggregate_O::static_className() % chem::Aggregate_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomReference_O @ %X") % classchem__AtomReference_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomReference_O::static_className() % chem::AtomReference_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Atom_O @ %X") % classchem__Atom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Atom_O::static_className() % chem::Atom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::BondListMatchNode_O @ %X") % classchem__BondListMatchNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::BondListMatchNode_O::static_className() % chem::BondListMatchNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::BondMatchNode_O @ %X") % classchem__BondMatchNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::BondMatchNode_O::static_className() % chem::BondMatchNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CalculatePositionAlongBond_O @ %X") % classchem__CalculatePositionAlongBond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CalculatePositionAlongBond_O::static_className() % chem::CalculatePositionAlongBond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CalculatePositionRelativeToOrigin_O @ %X") % classchem__CalculatePositionRelativeToOrigin_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CalculatePositionRelativeToOrigin_O::static_className() % chem::CalculatePositionRelativeToOrigin_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::CalculatePositionUsingInternals_O @ %X") % classchem__CalculatePositionUsingInternals_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::CalculatePositionUsingInternals_O::static_className() % chem::CalculatePositionUsingInternals_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ConstitutionVirtualAtom_O @ %X") % classchem__ConstitutionVirtualAtom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ConstitutionVirtualAtom_O::static_className() % chem::ConstitutionVirtualAtom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Constitution_O @ %X") % classchem__Constitution_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Constitution_O::static_className() % chem::Constitution_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::DirectionalCoupling_O @ %X") % classchem__DirectionalCoupling_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::DirectionalCoupling_O::static_className() % chem::DirectionalCoupling_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyAnchorRestraint_O @ %X") % classchem__EnergyAnchorRestraint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyAnchorRestraint_O::static_className() % chem::EnergyAnchorRestraint_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyAngle_O @ %X") % classchem__EnergyAngle_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyAngle_O::static_className() % chem::EnergyAngle_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyChiralRestraint_O @ %X") % classchem__EnergyChiralRestraint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyChiralRestraint_O::static_className() % chem::EnergyChiralRestraint_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyDihedral_O @ %X") % classchem__EnergyDihedral_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyDihedral_O::static_className() % chem::EnergyDihedral_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyFixedNonbondRestraint_O @ %X") % classchem__EnergyFixedNonbondRestraint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyFixedNonbondRestraint_O::static_className() % chem::EnergyFixedNonbondRestraint_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyImproperRestraint_O @ %X") % classchem__EnergyImproperRestraint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyImproperRestraint_O::static_className() % chem::EnergyImproperRestraint_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyNonbond_O @ %X") % classchem__EnergyNonbond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyNonbond_O::static_className() % chem::EnergyNonbond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EnergyStretch_O @ %X") % classchem__EnergyStretch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EnergyStretch_O::static_className() % chem::EnergyStretch_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EntityNameSetBase_O @ %X") % classchem__EntityNameSetBase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EntityNameSetBase_O::static_className() % chem::EntityNameSetBase_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFAngleDb_O @ %X") % classchem__FFAngleDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFAngleDb_O::static_className() % chem::FFAngleDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFAngle_O @ %X") % classchem__FFAngle_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFAngle_O::static_className() % chem::FFAngle_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFItorDb_O @ %X") % classchem__FFItorDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFItorDb_O::static_className() % chem::FFItorDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFItor_O @ %X") % classchem__FFItor_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFItor_O::static_className() % chem::FFItor_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFNonbondDb_O @ %X") % classchem__FFNonbondDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFNonbondDb_O::static_className() % chem::FFNonbondDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFNonbond_O @ %X") % classchem__FFNonbond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFNonbond_O::static_className() % chem::FFNonbond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFPtorDb_O @ %X") % classchem__FFPtorDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFPtorDb_O::static_className() % chem::FFPtorDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFPtor_O @ %X") % classchem__FFPtor_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFPtor_O::static_className() % chem::FFPtor_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFStretchDb_O @ %X") % classchem__FFStretchDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFStretchDb_O::static_className() % chem::FFStretchDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFStretch_O @ %X") % classchem__FFStretch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFStretch_O::static_className() % chem::FFStretch_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFTypesDb_O @ %X") % classchem__FFTypesDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFTypesDb_O::static_className() % chem::FFTypesDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFVdwDb_O @ %X") % classchem__FFVdwDb_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFVdwDb_O::static_className() % chem::FFVdwDb_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FFVdw_O @ %X") % classchem__FFVdw_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FFVdw_O::static_className() % chem::FFVdw_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ForceMatchReport_O @ %X") % classchem__ForceMatchReport_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ForceMatchReport_O::static_className() % chem::ForceMatchReport_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::FullLargeSquareMatrix_O @ %X") % classchem__FullLargeSquareMatrix_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::FullLargeSquareMatrix_O::static_className() % chem::FullLargeSquareMatrix_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::IterateMatter_O @ %X") % classchem__IterateMatter_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::IterateMatter_O::static_className() % chem::IterateMatter_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::IterateRestraints_O @ %X") % classchem__IterateRestraints_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::IterateRestraints_O::static_className() % chem::IterateRestraints_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::JumpPlug_O @ %X") % classchem__JumpPlug_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::JumpPlug_O::static_className() % chem::JumpPlug_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MessageReport_O @ %X") % classchem__MessageReport_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MessageReport_O::static_className() % chem::MessageReport_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Molecule_O @ %X") % classchem__Molecule_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Molecule_O::static_className() % chem::Molecule_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Monomer_O @ %X") % classchem__Monomer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Monomer_O::static_className() % chem::Monomer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::OligomerPart_Link_O @ %X") % classchem__OligomerPart_Link_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::OligomerPart_Link_O::static_className() % chem::OligomerPart_Link_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::OligomerPart_Monomer_O @ %X") % classchem__OligomerPart_Monomer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::OligomerPart_Monomer_O::static_className() % chem::OligomerPart_Monomer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::PlugWithMates_O @ %X") % classchem__PlugWithMates_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::PlugWithMates_O::static_className() % chem::PlugWithMates_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Residue_O @ %X") % classchem__Residue_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Residue_O::static_className() % chem::Residue_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestrainedExoCyclicAtom_O @ %X") % classchem__RestrainedExoCyclicAtom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestrainedExoCyclicAtom_O::static_className() % chem::RestrainedExoCyclicAtom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestrainedPiBond_O @ %X") % classchem__RestrainedPiBond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestrainedPiBond_O::static_className() % chem::RestrainedPiBond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintAnchor_O @ %X") % classchem__RestraintAnchor_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintAnchor_O::static_className() % chem::RestraintAnchor_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintAngle_O @ %X") % classchem__RestraintAngle_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintAngle_O::static_className() % chem::RestraintAngle_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintChiral_O @ %X") % classchem__RestraintChiral_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintChiral_O::static_className() % chem::RestraintChiral_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintDihedral_O @ %X") % classchem__RestraintDihedral_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintDihedral_O::static_className() % chem::RestraintDihedral_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintDistance_O @ %X") % classchem__RestraintDistance_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintDistance_O::static_className() % chem::RestraintDistance_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RestraintFixedNonbond_O @ %X") % classchem__RestraintFixedNonbond_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RestraintFixedNonbond_O::static_className() % chem::RestraintFixedNonbond_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RingCoupling_O @ %X") % classchem__RingCoupling_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RingCoupling_O::static_className() % chem::RingCoupling_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RootMatchNode_O @ %X") % classchem__RootMatchNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RootMatchNode_O::static_className() % chem::RootMatchNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SparseLargeSquareMatrix_O @ %X") % classchem__SparseLargeSquareMatrix_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SparseLargeSquareMatrix_O::static_className() % chem::SparseLargeSquareMatrix_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StepReport_O @ %X") % classchem__StepReport_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StepReport_O::static_className() % chem::StepReport_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::StereoisomerVirtualAtom_O @ %X") % classchem__StereoisomerVirtualAtom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::StereoisomerVirtualAtom_O::static_className() % chem::StereoisomerVirtualAtom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Stereoisomer_O @ %X") % classchem__Stereoisomer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Stereoisomer_O::static_className() % chem::Stereoisomer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SuperposableConformationCollection_O @ %X") % classchem__SuperposableConformationCollection_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SuperposableConformationCollection_O::static_className() % chem::SuperposableConformationCollection_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SuperposeSelectedAtoms_O @ %X") % classchem__SuperposeSelectedAtoms_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SuperposeSelectedAtoms_O::static_className() % chem::SuperposeSelectedAtoms_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ZMatrixAngleInternal_O @ %X") % classchem__ZMatrixAngleInternal_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ZMatrixAngleInternal_O::static_className() % chem::ZMatrixAngleInternal_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ZMatrixBondInternal_O @ %X") % classchem__ZMatrixBondInternal_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ZMatrixBondInternal_O::static_className() % chem::ZMatrixBondInternal_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ZMatrixDihedralInternal_O @ %X") % classchem__ZMatrixDihedralInternal_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ZMatrixDihedralInternal_O::static_className() % chem::ZMatrixDihedralInternal_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AfterMatchBondTest_O @ %X") % classchem__AfterMatchBondTest_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AfterMatchBondTest_O::static_className() % chem::AfterMatchBondTest_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AntechamberBondTest_O @ %X") % classchem__AntechamberBondTest_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AntechamberBondTest_O::static_className() % chem::AntechamberBondTest_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomOrBondMatchNode_O @ %X") % classchem__AtomOrBondMatchNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomOrBondMatchNode_O::static_className() % chem::AtomOrBondMatchNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::BondTest_O @ %X") % classchem__BondTest_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::BondTest_O::static_className() % chem::BondTest_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Branch_O @ %X") % classchem__Branch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Branch_O::static_className() % chem::Branch_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Chain_O @ %X") % classchem__Chain_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Chain_O::static_className() % chem::Chain_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EntityNameSetWithCap_O @ %X") % classchem__EntityNameSetWithCap_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EntityNameSetWithCap_O::static_className() % chem::EntityNameSetWithCap_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::EntityNameSet_O @ %X") % classchem__EntityNameSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::EntityNameSet_O::static_className() % chem::EntityNameSet_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::InPlug_O @ %X") % classchem__InPlug_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::InPlug_O::static_className() % chem::InPlug_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::IterateAtoms_O @ %X") % classchem__IterateAtoms_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::IterateAtoms_O::static_className() % chem::IterateAtoms_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::IterateBonds_O @ %X") % classchem__IterateBonds_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::IterateBonds_O::static_className() % chem::IterateBonds_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::IterateResidues_O @ %X") % classchem__IterateResidues_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::IterateResidues_O::static_className() % chem::IterateResidues_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MultiMonomer_O @ %X") % classchem__MultiMonomer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MultiMonomer_O::static_className() % chem::MultiMonomer_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::OutPlug_O @ %X") % classchem__OutPlug_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::OutPlug_O::static_className() % chem::OutPlug_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::VirtualAtom_O @ %X") % classchem__VirtualAtom_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::VirtualAtom_O::static_className() % chem::VirtualAtom_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AntechamberFocusAtomMatch_O @ %X") % classchem__AntechamberFocusAtomMatch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AntechamberFocusAtomMatch_O::static_className() % chem::AntechamberFocusAtomMatch_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AtomTest_O @ %X") % classchem__AtomTest_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AtomTest_O::static_className() % chem::AtomTest_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Logical_O @ %X") % classchem__Logical_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Logical_O::static_className() % chem::Logical_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Mate_O @ %X") % classchem__Mate_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Mate_O::static_className() % chem::Mate_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::MonomerPack_O @ %X") % classchem__MonomerPack_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::MonomerPack_O::static_className() % chem::MonomerPack_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RepresentedEntityNameSet_O @ %X") % classchem__RepresentedEntityNameSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RepresentedEntityNameSet_O::static_className() % chem::RepresentedEntityNameSet_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::ResidueTest_O @ %X") % classchem__ResidueTest_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::ResidueTest_O::static_className() % chem::ResidueTest_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RingClosingMate_O @ %X") % classchem__RingClosingMate_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RingClosingMate_O::static_className() % chem::RingClosingMate_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RingClosingPlug_O @ %X") % classchem__RingClosingPlug_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RingClosingPlug_O::static_className() % chem::RingClosingPlug_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::RingTest_O @ %X") % classchem__RingTest_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::RingTest_O::static_className() % chem::RingTest_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::Root_O @ %X") % classchem__Root_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::Root_O::static_className() % chem::Root_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::TagSet_O @ %X") % classchem__TagSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::TagSet_O::static_className() % chem::TagSet_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::AntechamberRoot_O @ %X") % classchem__AntechamberRoot_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::AntechamberRoot_O::static_className() % chem::AntechamberRoot_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: chem::SmartsRoot_O @ %X") % classchem__SmartsRoot_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % chem::SmartsRoot_O::static_className() % chem::SmartsRoot_O::static_classSymbol() );
#endif // } DUMP_INFO_CLASS
#undef DUMP_INFO_CLASS
#if defined(DEFINE_BASE_CLASSES) || defined(ALL_STAGES) // {
// Depends on nothing
classchem__AGEdge_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AGVertex_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AbstractLargeSquareMatrix_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Alias_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Angle_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AtomIdToAtomMap_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AtomId_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AtomIndexer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AtomReferenceBase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__AtomTable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__BondList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Bond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CDBond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CDFragment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CDNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CDText_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CalculatePosition_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CandoDatabase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ChemDraw_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ChemInfoMatch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ChemInfoNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ChemInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CipPrioritizer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Command_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ComplexRestraint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConformationCollectionEntry_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConformationCollection_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConformationExplorerEntryStage_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConformationExplorerEntry_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConformationExplorerMatch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConformationExplorer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConstitutionAtom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConstitutionAtoms_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ConstitutionBond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__CoordinateSystem_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Coupling_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__EnergyComponent_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__EnergyFunction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Entity_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__FFBaseDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__FFNonbondCrossTermTable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__FFParameter_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ForceField_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__FrameRecognizer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ImproperTorsion_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__InfoDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__MacroModelFile_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__MapOfMonomerNamesToAtomIndexers_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Matter_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__MinimizerLog_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Minimizer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__MonomerContext_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__MonomerCoordinates_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__NVector_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__NumericalFunction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__OligomerPart_Base_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Oligomer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__OneMonomer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__OneTypeRule_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__PathMessage_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__PdbMonomerConnectivity_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__PdbMonomerDatabase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__PdbReader_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__PdbWriter_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Plug_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ProperTorsion_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ReadAmberParameters_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ReportBase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ResidueList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__RestraintList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
classchem__Restraint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__RingFinder_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Rotamer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__SearchStatistics_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__SpanningLoop_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__SpecificContextSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__SpecificContext_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__StatusTracker_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__StereoConfiguration_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__StereoInformation_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__StereoisomerAtom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__StereoisomerAtoms_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__StructureComparer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Structure_Old_ListEntry_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Structure_Old_List_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__SuperposeEngine_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Topology_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__TrajectoryFrame_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Trajectory_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__TwisterDriver_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Twister_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__TypeAssignmentRules_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__VirtualSphere_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__WildElementDict_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ZMatrixEntry_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ZMatrixInternal_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__ZMatrix_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classchem__Aggregate_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Matter_O::static_classSymbol());
classchem__AtomReference_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomReferenceBase_O::static_classSymbol());
classchem__Atom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Matter_O::static_classSymbol());
classchem__BondListMatchNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ChemInfoNode_O::static_classSymbol());
classchem__BondMatchNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ChemInfoNode_O::static_classSymbol());
classchem__CalculatePositionAlongBond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::CalculatePosition_O::static_classSymbol());
classchem__CalculatePositionRelativeToOrigin_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::CalculatePosition_O::static_classSymbol());
classchem__CalculatePositionUsingInternals_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::CalculatePosition_O::static_classSymbol());
classchem__ConstitutionVirtualAtom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ConstitutionAtom_O::static_classSymbol());
classchem__Constitution_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Entity_O::static_classSymbol());
classchem__DirectionalCoupling_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Coupling_O::static_classSymbol());
classchem__EnergyAnchorRestraint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyAngle_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyChiralRestraint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyDihedral_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyFixedNonbondRestraint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyImproperRestraint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyNonbond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EnergyStretch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());
classchem__EntityNameSetBase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Entity_O::static_classSymbol());
classchem__FFAngleDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFAngle_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());
classchem__FFItorDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFItor_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());
classchem__FFNonbondDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFNonbond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());
classchem__FFPtorDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFPtor_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());
classchem__FFStretchDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFStretch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());
classchem__FFTypesDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFVdwDb_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());
classchem__FFVdw_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());
classchem__ForceMatchReport_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ReportBase_O::static_classSymbol());
classchem__FullLargeSquareMatrix_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AbstractLargeSquareMatrix_O::static_classSymbol());
classchem__IterateMatter_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
classchem__IterateRestraints_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
classchem__JumpPlug_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Plug_O::static_classSymbol());
classchem__MessageReport_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ReportBase_O::static_classSymbol());
classchem__Molecule_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Matter_O::static_classSymbol());
classchem__Monomer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Entity_O::static_classSymbol());
classchem__OligomerPart_Link_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::OligomerPart_Base_O::static_classSymbol());
classchem__OligomerPart_Monomer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::OligomerPart_Base_O::static_classSymbol());
classchem__PlugWithMates_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Plug_O::static_classSymbol());
classchem__Residue_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Matter_O::static_classSymbol());
classchem__RestrainedExoCyclicAtom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ComplexRestraint_O::static_classSymbol());
classchem__RestrainedPiBond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ComplexRestraint_O::static_classSymbol());
classchem__RestraintAnchor_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());
classchem__RestraintAngle_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());
classchem__RestraintChiral_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());
classchem__RestraintDihedral_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());
classchem__RestraintDistance_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());
classchem__RestraintFixedNonbond_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());
classchem__RingCoupling_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Coupling_O::static_classSymbol());
classchem__RootMatchNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ChemInfoNode_O::static_classSymbol());
classchem__SparseLargeSquareMatrix_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AbstractLargeSquareMatrix_O::static_classSymbol());
classchem__StepReport_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ReportBase_O::static_classSymbol());
classchem__StereoisomerVirtualAtom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::StereoisomerAtom_O::static_classSymbol());
classchem__Stereoisomer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Entity_O::static_classSymbol());
classchem__SuperposableConformationCollection_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ConformationCollection_O::static_classSymbol());
classchem__SuperposeSelectedAtoms_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::SuperposeEngine_O::static_classSymbol());
classchem__ZMatrixAngleInternal_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ZMatrixInternal_O::static_classSymbol());
classchem__ZMatrixBondInternal_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ZMatrixInternal_O::static_classSymbol());
classchem__ZMatrixDihedralInternal_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::ZMatrixInternal_O::static_classSymbol());
classchem__AfterMatchBondTest_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::RootMatchNode_O::static_classSymbol());
classchem__AntechamberBondTest_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::BondMatchNode_O::static_classSymbol());
classchem__AtomOrBondMatchNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::BondMatchNode_O::static_classSymbol());
classchem__BondTest_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::BondMatchNode_O::static_classSymbol());
classchem__Branch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::BondListMatchNode_O::static_classSymbol());
classchem__Chain_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::BondListMatchNode_O::static_classSymbol());
classchem__EntityNameSetWithCap_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EntityNameSetBase_O::static_classSymbol());
classchem__EntityNameSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EntityNameSetBase_O::static_classSymbol());
classchem__InPlug_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::PlugWithMates_O::static_classSymbol());
classchem__IterateAtoms_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::IterateMatter_O::static_classSymbol());
classchem__IterateBonds_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::IterateMatter_O::static_classSymbol());
classchem__IterateResidues_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::IterateMatter_O::static_classSymbol());
classchem__MultiMonomer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Monomer_O::static_classSymbol());
classchem__OutPlug_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::PlugWithMates_O::static_classSymbol());
classchem__VirtualAtom_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Atom_O::static_classSymbol());
classchem__AntechamberFocusAtomMatch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__AtomTest_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__Logical_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__Mate_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EntityNameSetWithCap_O::static_classSymbol());
classchem__MonomerPack_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EntityNameSet_O::static_classSymbol());
classchem__RepresentedEntityNameSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EntityNameSet_O::static_classSymbol());
classchem__ResidueTest_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__RingClosingMate_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::EntityNameSetWithCap_O::static_classSymbol());
classchem__RingClosingPlug_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::OutPlug_O::static_classSymbol());
classchem__RingTest_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__Root_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__TagSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());
classchem__AntechamberRoot_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Root_O::static_classSymbol());
classchem__SmartsRoot_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(chem::Root_O::static_classSymbol());
#endif // } DEFINE_BASE_CLASSES
#undef DEFINE_BASE_CLASSES
#if defined(DEFINE_CLASS_NAMES) || defined(ALL_STAGES) // {
// Depends on nothing

    classchem__AGEdge_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AGEdge_O::static_classSymbol());

    classchem__AGVertex_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AGVertex_O::static_classSymbol());

    classchem__AbstractLargeSquareMatrix_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AbstractLargeSquareMatrix_O::static_classSymbol());

    classchem__Alias_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Alias_O::static_classSymbol());

    classchem__Angle_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Angle_O::static_classSymbol());

    classchem__AtomIdToAtomMap_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomIdToAtomMap_O::static_classSymbol());

    classchem__AtomId_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomId_O::static_classSymbol());

    classchem__AtomIndexer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomIndexer_O::static_classSymbol());

    classchem__AtomReferenceBase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomReferenceBase_O::static_classSymbol());

    classchem__AtomTable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomTable_O::static_classSymbol());

    classchem__BondList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::BondList_O::static_classSymbol());

    classchem__Bond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Bond_O::static_classSymbol());

    classchem__CDBond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CDBond_O::static_classSymbol());

    classchem__CDFragment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CDFragment_O::static_classSymbol());

    classchem__CDNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CDNode_O::static_classSymbol());

    classchem__CDText_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CDText_O::static_classSymbol());

    classchem__CalculatePosition_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CalculatePosition_O::static_classSymbol());

    classchem__CandoDatabase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CandoDatabase_O::static_classSymbol());

    classchem__ChemDraw_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ChemDraw_O::static_classSymbol());

    classchem__ChemInfoMatch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ChemInfoMatch_O::static_classSymbol());

    classchem__ChemInfoNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ChemInfoNode_O::static_classSymbol());

    classchem__ChemInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ChemInfo_O::static_classSymbol());

    classchem__CipPrioritizer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CipPrioritizer_O::static_classSymbol());

    classchem__Command_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Command_O::static_classSymbol());

    classchem__ComplexRestraint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ComplexRestraint_O::static_classSymbol());

    classchem__ConformationCollectionEntry_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConformationCollectionEntry_O::static_classSymbol());

    classchem__ConformationCollection_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConformationCollection_O::static_classSymbol());

    classchem__ConformationExplorerEntryStage_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConformationExplorerEntryStage_O::static_classSymbol());

    classchem__ConformationExplorerEntry_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConformationExplorerEntry_O::static_classSymbol());

    classchem__ConformationExplorerMatch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConformationExplorerMatch_O::static_classSymbol());

    classchem__ConformationExplorer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConformationExplorer_O::static_classSymbol());

    classchem__ConstitutionAtom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConstitutionAtom_O::static_classSymbol());

    classchem__ConstitutionAtoms_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConstitutionAtoms_O::static_classSymbol());

    classchem__ConstitutionBond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConstitutionBond_O::static_classSymbol());

    classchem__CoordinateSystem_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CoordinateSystem_O::static_classSymbol());

    classchem__Coupling_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Coupling_O::static_classSymbol());

    classchem__EnergyComponent_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyComponent_O::static_classSymbol());

    classchem__EnergyFunction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyFunction_O::static_classSymbol());

    classchem__Entity_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Entity_O::static_classSymbol());

    classchem__FFBaseDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFBaseDb_O::static_classSymbol());

    classchem__FFNonbondCrossTermTable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFNonbondCrossTermTable_O::static_classSymbol());

    classchem__FFParameter_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFParameter_O::static_classSymbol());

    classchem__ForceField_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ForceField_O::static_classSymbol());

    classchem__FrameRecognizer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FrameRecognizer_O::static_classSymbol());

    classchem__ImproperTorsion_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ImproperTorsion_O::static_classSymbol());

    classchem__InfoDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::InfoDb_O::static_classSymbol());

    classchem__MacroModelFile_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MacroModelFile_O::static_classSymbol());

    classchem__MapOfMonomerNamesToAtomIndexers_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MapOfMonomerNamesToAtomIndexers_O::static_classSymbol());

    classchem__Matter_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Matter_O::static_classSymbol());

    classchem__MinimizerLog_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MinimizerLog_O::static_classSymbol());

    classchem__Minimizer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Minimizer_O::static_classSymbol());

    classchem__MonomerContext_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MonomerContext_O::static_classSymbol());

    classchem__MonomerCoordinates_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MonomerCoordinates_O::static_classSymbol());

    classchem__NVector_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::NVector_O::static_classSymbol());

    classchem__NumericalFunction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::NumericalFunction_O::static_classSymbol());

    classchem__OligomerPart_Base_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::OligomerPart_Base_O::static_classSymbol());

    classchem__Oligomer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Oligomer_O::static_classSymbol());

    classchem__OneMonomer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::OneMonomer_O::static_classSymbol());

    classchem__OneTypeRule_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::OneTypeRule_O::static_classSymbol());

    classchem__PathMessage_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::PathMessage_O::static_classSymbol());

    classchem__PdbMonomerConnectivity_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::PdbMonomerConnectivity_O::static_classSymbol());

    classchem__PdbMonomerDatabase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::PdbMonomerDatabase_O::static_classSymbol());

    classchem__PdbReader_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::PdbReader_O::static_classSymbol());

    classchem__PdbWriter_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::PdbWriter_O::static_classSymbol());

    classchem__Plug_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Plug_O::static_classSymbol());

    classchem__ProperTorsion_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ProperTorsion_O::static_classSymbol());

    classchem__ReadAmberParameters_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ReadAmberParameters_O::static_classSymbol());

    classchem__ReportBase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ReportBase_O::static_classSymbol());

    classchem__ResidueList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ResidueList_O::static_classSymbol());

    classchem__RestraintList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintList_O::static_classSymbol());

    classchem__Restraint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Restraint_O::static_classSymbol());

    classchem__RingFinder_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RingFinder_O::static_classSymbol());

    classchem__Rotamer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Rotamer_O::static_classSymbol());

    classchem__SearchStatistics_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SearchStatistics_O::static_classSymbol());

    classchem__SpanningLoop_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SpanningLoop_O::static_classSymbol());

    classchem__SpecificContextSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SpecificContextSet_O::static_classSymbol());

    classchem__SpecificContext_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SpecificContext_O::static_classSymbol());

    classchem__StatusTracker_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StatusTracker_O::static_classSymbol());

    classchem__StereoConfiguration_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StereoConfiguration_O::static_classSymbol());

    classchem__StereoInformation_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StereoInformation_O::static_classSymbol());

    classchem__StereoisomerAtom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StereoisomerAtom_O::static_classSymbol());

    classchem__StereoisomerAtoms_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StereoisomerAtoms_O::static_classSymbol());

    classchem__StructureComparer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StructureComparer_O::static_classSymbol());

    classchem__Structure_Old_ListEntry_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Structure_Old_ListEntry_O::static_classSymbol());

    classchem__Structure_Old_List_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Structure_Old_List_O::static_classSymbol());

    classchem__SuperposeEngine_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SuperposeEngine_O::static_classSymbol());

    classchem__Topology_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Topology_O::static_classSymbol());

    classchem__TrajectoryFrame_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::TrajectoryFrame_O::static_classSymbol());

    classchem__Trajectory_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Trajectory_O::static_classSymbol());

    classchem__TwisterDriver_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::TwisterDriver_O::static_classSymbol());

    classchem__Twister_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Twister_O::static_classSymbol());

    classchem__TypeAssignmentRules_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::TypeAssignmentRules_O::static_classSymbol());

    classchem__VirtualSphere_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::VirtualSphere_O::static_classSymbol());

    classchem__WildElementDict_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::WildElementDict_O::static_classSymbol());

    classchem__ZMatrixEntry_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ZMatrixEntry_O::static_classSymbol());

    classchem__ZMatrixInternal_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ZMatrixInternal_O::static_classSymbol());

    classchem__ZMatrix_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ZMatrix_O::static_classSymbol());

    classchem__Aggregate_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Aggregate_O::static_classSymbol());

    classchem__AtomReference_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomReference_O::static_classSymbol());

    classchem__Atom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Atom_O::static_classSymbol());

    classchem__BondListMatchNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::BondListMatchNode_O::static_classSymbol());

    classchem__BondMatchNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::BondMatchNode_O::static_classSymbol());

    classchem__CalculatePositionAlongBond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CalculatePositionAlongBond_O::static_classSymbol());

    classchem__CalculatePositionRelativeToOrigin_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CalculatePositionRelativeToOrigin_O::static_classSymbol());

    classchem__CalculatePositionUsingInternals_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::CalculatePositionUsingInternals_O::static_classSymbol());

    classchem__ConstitutionVirtualAtom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ConstitutionVirtualAtom_O::static_classSymbol());

    classchem__Constitution_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Constitution_O::static_classSymbol());

    classchem__DirectionalCoupling_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::DirectionalCoupling_O::static_classSymbol());

    classchem__EnergyAnchorRestraint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyAnchorRestraint_O::static_classSymbol());

    classchem__EnergyAngle_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyAngle_O::static_classSymbol());

    classchem__EnergyChiralRestraint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyChiralRestraint_O::static_classSymbol());

    classchem__EnergyDihedral_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyDihedral_O::static_classSymbol());

    classchem__EnergyFixedNonbondRestraint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyFixedNonbondRestraint_O::static_classSymbol());

    classchem__EnergyImproperRestraint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyImproperRestraint_O::static_classSymbol());

    classchem__EnergyNonbond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyNonbond_O::static_classSymbol());

    classchem__EnergyStretch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EnergyStretch_O::static_classSymbol());

    classchem__EntityNameSetBase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EntityNameSetBase_O::static_classSymbol());

    classchem__FFAngleDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFAngleDb_O::static_classSymbol());

    classchem__FFAngle_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFAngle_O::static_classSymbol());

    classchem__FFItorDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFItorDb_O::static_classSymbol());

    classchem__FFItor_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFItor_O::static_classSymbol());

    classchem__FFNonbondDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFNonbondDb_O::static_classSymbol());

    classchem__FFNonbond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFNonbond_O::static_classSymbol());

    classchem__FFPtorDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFPtorDb_O::static_classSymbol());

    classchem__FFPtor_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFPtor_O::static_classSymbol());

    classchem__FFStretchDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFStretchDb_O::static_classSymbol());

    classchem__FFStretch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFStretch_O::static_classSymbol());

    classchem__FFTypesDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFTypesDb_O::static_classSymbol());

    classchem__FFVdwDb_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFVdwDb_O::static_classSymbol());

    classchem__FFVdw_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FFVdw_O::static_classSymbol());

    classchem__ForceMatchReport_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ForceMatchReport_O::static_classSymbol());

    classchem__FullLargeSquareMatrix_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::FullLargeSquareMatrix_O::static_classSymbol());

    classchem__IterateMatter_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::IterateMatter_O::static_classSymbol());

    classchem__IterateRestraints_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::IterateRestraints_O::static_classSymbol());

    classchem__JumpPlug_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::JumpPlug_O::static_classSymbol());

    classchem__MessageReport_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MessageReport_O::static_classSymbol());

    classchem__Molecule_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Molecule_O::static_classSymbol());

    classchem__Monomer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Monomer_O::static_classSymbol());

    classchem__OligomerPart_Link_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::OligomerPart_Link_O::static_classSymbol());

    classchem__OligomerPart_Monomer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::OligomerPart_Monomer_O::static_classSymbol());

    classchem__PlugWithMates_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::PlugWithMates_O::static_classSymbol());

    classchem__Residue_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Residue_O::static_classSymbol());

    classchem__RestrainedExoCyclicAtom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestrainedExoCyclicAtom_O::static_classSymbol());

    classchem__RestrainedPiBond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestrainedPiBond_O::static_classSymbol());

    classchem__RestraintAnchor_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintAnchor_O::static_classSymbol());

    classchem__RestraintAngle_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintAngle_O::static_classSymbol());

    classchem__RestraintChiral_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintChiral_O::static_classSymbol());

    classchem__RestraintDihedral_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintDihedral_O::static_classSymbol());

    classchem__RestraintDistance_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintDistance_O::static_classSymbol());

    classchem__RestraintFixedNonbond_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RestraintFixedNonbond_O::static_classSymbol());

    classchem__RingCoupling_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RingCoupling_O::static_classSymbol());

    classchem__RootMatchNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RootMatchNode_O::static_classSymbol());

    classchem__SparseLargeSquareMatrix_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SparseLargeSquareMatrix_O::static_classSymbol());

    classchem__StepReport_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StepReport_O::static_classSymbol());

    classchem__StereoisomerVirtualAtom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::StereoisomerVirtualAtom_O::static_classSymbol());

    classchem__Stereoisomer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Stereoisomer_O::static_classSymbol());

    classchem__SuperposableConformationCollection_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SuperposableConformationCollection_O::static_classSymbol());

    classchem__SuperposeSelectedAtoms_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SuperposeSelectedAtoms_O::static_classSymbol());

    classchem__ZMatrixAngleInternal_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ZMatrixAngleInternal_O::static_classSymbol());

    classchem__ZMatrixBondInternal_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ZMatrixBondInternal_O::static_classSymbol());

    classchem__ZMatrixDihedralInternal_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ZMatrixDihedralInternal_O::static_classSymbol());

    classchem__AfterMatchBondTest_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AfterMatchBondTest_O::static_classSymbol());

    classchem__AntechamberBondTest_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AntechamberBondTest_O::static_classSymbol());

    classchem__AtomOrBondMatchNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomOrBondMatchNode_O::static_classSymbol());

    classchem__BondTest_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::BondTest_O::static_classSymbol());

    classchem__Branch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Branch_O::static_classSymbol());

    classchem__Chain_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Chain_O::static_classSymbol());

    classchem__EntityNameSetWithCap_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EntityNameSetWithCap_O::static_classSymbol());

    classchem__EntityNameSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::EntityNameSet_O::static_classSymbol());

    classchem__InPlug_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::InPlug_O::static_classSymbol());

    classchem__IterateAtoms_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::IterateAtoms_O::static_classSymbol());

    classchem__IterateBonds_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::IterateBonds_O::static_classSymbol());

    classchem__IterateResidues_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::IterateResidues_O::static_classSymbol());

    classchem__MultiMonomer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MultiMonomer_O::static_classSymbol());

    classchem__OutPlug_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::OutPlug_O::static_classSymbol());

    classchem__VirtualAtom_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::VirtualAtom_O::static_classSymbol());

    classchem__AntechamberFocusAtomMatch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AntechamberFocusAtomMatch_O::static_classSymbol());

    classchem__AtomTest_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AtomTest_O::static_classSymbol());

    classchem__Logical_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Logical_O::static_classSymbol());

    classchem__Mate_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Mate_O::static_classSymbol());

    classchem__MonomerPack_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::MonomerPack_O::static_classSymbol());

    classchem__RepresentedEntityNameSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RepresentedEntityNameSet_O::static_classSymbol());

    classchem__ResidueTest_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::ResidueTest_O::static_classSymbol());

    classchem__RingClosingMate_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RingClosingMate_O::static_classSymbol());

    classchem__RingClosingPlug_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RingClosingPlug_O::static_classSymbol());

    classchem__RingTest_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::RingTest_O::static_classSymbol());

    classchem__Root_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::Root_O::static_classSymbol());

    classchem__TagSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::TagSet_O::static_classSymbol());

    classchem__AntechamberRoot_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::AntechamberRoot_O::static_classSymbol());

    classchem__SmartsRoot_Oval->__setupStage3NameAndCalculateClassPrecedenceList(chem::SmartsRoot_O::static_classSymbol());
#endif // } DEFINE_CLASS_NAMES
#undef DEFINE_CLASS_NAMES
#if defined(EXPOSE_TO_CANDO) || defined(ALL_STAGES)
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AGEdge_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AGEdge_O");
    chem::Register_chem__AGEdge_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AGVertex_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AGVertex_O");
    chem::Register_chem__AGVertex_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AbstractLargeSquareMatrix_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AbstractLargeSquareMatrix_O");
    chem::Register_chem__AbstractLargeSquareMatrix_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Alias_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Alias_O");
    chem::Register_chem__Alias_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Angle_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Angle_O");
    chem::Register_chem__Angle_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomIdToAtomMap_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomIdToAtomMap_O");
    chem::Register_chem__AtomIdToAtomMap_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomId_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomId_O");
    chem::Register_chem__AtomId_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomIndexer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomIndexer_O");
    chem::Register_chem__AtomIndexer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomReferenceBase_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomReferenceBase_O");
    chem::Register_chem__AtomReferenceBase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomTable_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomTable_O");
    chem::Register_chem__AtomTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__BondList_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__BondList_O");
    chem::Register_chem__BondList_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Bond_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Bond_O");
    chem::Register_chem__Bond_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CDBond_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CDBond_O");
    chem::Register_chem__CDBond_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CDFragment_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CDFragment_O");
    chem::Register_chem__CDFragment_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CDNode_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CDNode_O");
    chem::Register_chem__CDNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CDText_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CDText_O");
    chem::Register_chem__CDText_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CalculatePosition_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CalculatePosition_O");
    chem::Register_chem__CalculatePosition_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CandoDatabase_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CandoDatabase_O");
    chem::Register_chem__CandoDatabase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ChemDraw_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ChemDraw_O");
    chem::Register_chem__ChemDraw_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ChemInfoMatch_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ChemInfoMatch_O");
    chem::Register_chem__ChemInfoMatch_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ChemInfoNode_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ChemInfoNode_O");
    chem::Register_chem__ChemInfoNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ChemInfo_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ChemInfo_O");
    chem::Register_chem__ChemInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CipPrioritizer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CipPrioritizer_O");
    chem::Register_chem__CipPrioritizer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Command_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Command_O");
    chem::Register_chem__Command_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ComplexRestraint_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ComplexRestraint_O");
    chem::Register_chem__ComplexRestraint_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConformationCollectionEntry_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConformationCollectionEntry_O");
    chem::Register_chem__ConformationCollectionEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConformationCollection_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConformationCollection_O");
    chem::Register_chem__ConformationCollection_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConformationExplorerEntryStage_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConformationExplorerEntryStage_O");
    chem::Register_chem__ConformationExplorerEntryStage_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConformationExplorerEntry_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConformationExplorerEntry_O");
    chem::Register_chem__ConformationExplorerEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConformationExplorerMatch_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConformationExplorerMatch_O");
    chem::Register_chem__ConformationExplorerMatch_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConformationExplorer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConformationExplorer_O");
    chem::Register_chem__ConformationExplorer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConstitutionAtom_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConstitutionAtom_O");
    chem::Register_chem__ConstitutionAtom_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConstitutionAtoms_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConstitutionAtoms_O");
    chem::Register_chem__ConstitutionAtoms_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConstitutionBond_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConstitutionBond_O");
    chem::Register_chem__ConstitutionBond_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CoordinateSystem_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CoordinateSystem_O");
    chem::Register_chem__CoordinateSystem_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Coupling_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Coupling_O");
    chem::Register_chem__Coupling_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyComponent_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyComponent_O");
    chem::Register_chem__EnergyComponent_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyFunction_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyFunction_O");
    chem::Register_chem__EnergyFunction_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Entity_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Entity_O");
    chem::Register_chem__Entity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFBaseDb_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFBaseDb_O");
    chem::Register_chem__FFBaseDb_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFNonbondCrossTermTable_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFNonbondCrossTermTable_O");
    chem::Register_chem__FFNonbondCrossTermTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFParameter_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFParameter_O");
    chem::Register_chem__FFParameter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ForceField_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ForceField_O");
    chem::Register_chem__ForceField_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FrameRecognizer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FrameRecognizer_O");
    chem::Register_chem__FrameRecognizer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ImproperTorsion_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ImproperTorsion_O");
    chem::Register_chem__ImproperTorsion_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__InfoDb_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__InfoDb_O");
    chem::Register_chem__InfoDb_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MacroModelFile_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MacroModelFile_O");
    chem::Register_chem__MacroModelFile_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MapOfMonomerNamesToAtomIndexers_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MapOfMonomerNamesToAtomIndexers_O");
    chem::Register_chem__MapOfMonomerNamesToAtomIndexers_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Matter_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Matter_O");
    chem::Register_chem__Matter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MinimizerLog_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MinimizerLog_O");
    chem::Register_chem__MinimizerLog_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Minimizer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Minimizer_O");
    chem::Register_chem__Minimizer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MonomerContext_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MonomerContext_O");
    chem::Register_chem__MonomerContext_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MonomerCoordinates_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MonomerCoordinates_O");
    chem::Register_chem__MonomerCoordinates_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__NVector_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__NVector_O");
    chem::Register_chem__NVector_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__NumericalFunction_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__NumericalFunction_O");
    chem::Register_chem__NumericalFunction_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__OligomerPart_Base_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__OligomerPart_Base_O");
    chem::Register_chem__OligomerPart_Base_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Oligomer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Oligomer_O");
    chem::Register_chem__Oligomer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__OneMonomer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__OneMonomer_O");
    chem::Register_chem__OneMonomer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__OneTypeRule_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__OneTypeRule_O");
    chem::Register_chem__OneTypeRule_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__PathMessage_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__PathMessage_O");
    chem::Register_chem__PathMessage_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__PdbMonomerConnectivity_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__PdbMonomerConnectivity_O");
    chem::Register_chem__PdbMonomerConnectivity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__PdbMonomerDatabase_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__PdbMonomerDatabase_O");
    chem::Register_chem__PdbMonomerDatabase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__PdbReader_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__PdbReader_O");
    chem::Register_chem__PdbReader_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__PdbWriter_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__PdbWriter_O");
    chem::Register_chem__PdbWriter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Plug_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Plug_O");
    chem::Register_chem__Plug_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ProperTorsion_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ProperTorsion_O");
    chem::Register_chem__ProperTorsion_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ReadAmberParameters_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ReadAmberParameters_O");
    chem::Register_chem__ReadAmberParameters_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ReportBase_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ReportBase_O");
    chem::Register_chem__ReportBase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ResidueList_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ResidueList_O");
    chem::Register_chem__ResidueList_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintList_O(core::Lisp_sp); // base(s): set(['core::Vector_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintList_O");
    chem::Register_chem__RestraintList_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Restraint_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Restraint_O");
    chem::Register_chem__Restraint_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RingFinder_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RingFinder_O");
    chem::Register_chem__RingFinder_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Rotamer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Rotamer_O");
    chem::Register_chem__Rotamer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SearchStatistics_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SearchStatistics_O");
    chem::Register_chem__SearchStatistics_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SpanningLoop_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SpanningLoop_O");
    chem::Register_chem__SpanningLoop_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SpecificContextSet_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SpecificContextSet_O");
    chem::Register_chem__SpecificContextSet_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SpecificContext_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SpecificContext_O");
    chem::Register_chem__SpecificContext_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StatusTracker_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StatusTracker_O");
    chem::Register_chem__StatusTracker_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StereoConfiguration_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StereoConfiguration_O");
    chem::Register_chem__StereoConfiguration_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StereoInformation_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StereoInformation_O");
    chem::Register_chem__StereoInformation_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StereoisomerAtom_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StereoisomerAtom_O");
    chem::Register_chem__StereoisomerAtom_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StereoisomerAtoms_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StereoisomerAtoms_O");
    chem::Register_chem__StereoisomerAtoms_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StructureComparer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StructureComparer_O");
    chem::Register_chem__StructureComparer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Structure_Old_ListEntry_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Structure_Old_ListEntry_O");
    chem::Register_chem__Structure_Old_ListEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Structure_Old_List_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Structure_Old_List_O");
    chem::Register_chem__Structure_Old_List_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SuperposeEngine_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SuperposeEngine_O");
    chem::Register_chem__SuperposeEngine_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Topology_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Topology_O");
    chem::Register_chem__Topology_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__TrajectoryFrame_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__TrajectoryFrame_O");
    chem::Register_chem__TrajectoryFrame_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Trajectory_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Trajectory_O");
    chem::Register_chem__Trajectory_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__TwisterDriver_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__TwisterDriver_O");
    chem::Register_chem__TwisterDriver_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Twister_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Twister_O");
    chem::Register_chem__Twister_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__TypeAssignmentRules_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__TypeAssignmentRules_O");
    chem::Register_chem__TypeAssignmentRules_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__VirtualSphere_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__VirtualSphere_O");
    chem::Register_chem__VirtualSphere_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__WildElementDict_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__WildElementDict_O");
    chem::Register_chem__WildElementDict_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ZMatrixEntry_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ZMatrixEntry_O");
    chem::Register_chem__ZMatrixEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ZMatrixInternal_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ZMatrixInternal_O");
    chem::Register_chem__ZMatrixInternal_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ZMatrix_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ZMatrix_O");
    chem::Register_chem__ZMatrix_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Aggregate_O(core::Lisp_sp); // base(s): set(['chem::Matter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Aggregate_O");
    chem::Register_chem__Aggregate_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomReference_O(core::Lisp_sp); // base(s): set(['chem::AtomReferenceBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomReference_O");
    chem::Register_chem__AtomReference_O(_lisp); // base(s): set(['chem::AtomReferenceBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Atom_O(core::Lisp_sp); // base(s): set(['chem::Matter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Atom_O");
    chem::Register_chem__Atom_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__BondListMatchNode_O(core::Lisp_sp); // base(s): set(['chem::ChemInfoNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__BondListMatchNode_O");
    chem::Register_chem__BondListMatchNode_O(_lisp); // base(s): set(['chem::ChemInfoNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__BondMatchNode_O(core::Lisp_sp); // base(s): set(['chem::ChemInfoNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__BondMatchNode_O");
    chem::Register_chem__BondMatchNode_O(_lisp); // base(s): set(['chem::ChemInfoNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CalculatePositionAlongBond_O(core::Lisp_sp); // base(s): set(['chem::CalculatePosition_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CalculatePositionAlongBond_O");
    chem::Register_chem__CalculatePositionAlongBond_O(_lisp); // base(s): set(['chem::CalculatePosition_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CalculatePositionRelativeToOrigin_O(core::Lisp_sp); // base(s): set(['chem::CalculatePosition_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CalculatePositionRelativeToOrigin_O");
    chem::Register_chem__CalculatePositionRelativeToOrigin_O(_lisp); // base(s): set(['chem::CalculatePosition_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__CalculatePositionUsingInternals_O(core::Lisp_sp); // base(s): set(['chem::CalculatePosition_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__CalculatePositionUsingInternals_O");
    chem::Register_chem__CalculatePositionUsingInternals_O(_lisp); // base(s): set(['chem::CalculatePosition_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ConstitutionVirtualAtom_O(core::Lisp_sp); // base(s): set(['chem::ConstitutionAtom_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ConstitutionVirtualAtom_O");
    chem::Register_chem__ConstitutionVirtualAtom_O(_lisp); // base(s): set(['chem::ConstitutionAtom_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Constitution_O(core::Lisp_sp); // base(s): set(['chem::Entity_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Constitution_O");
    chem::Register_chem__Constitution_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__DirectionalCoupling_O(core::Lisp_sp); // base(s): set(['chem::Coupling_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__DirectionalCoupling_O");
    chem::Register_chem__DirectionalCoupling_O(_lisp); // base(s): set(['chem::Coupling_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyAnchorRestraint_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyAnchorRestraint_O");
    chem::Register_chem__EnergyAnchorRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyAngle_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyAngle_O");
    chem::Register_chem__EnergyAngle_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyChiralRestraint_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyChiralRestraint_O");
    chem::Register_chem__EnergyChiralRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyDihedral_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyDihedral_O");
    chem::Register_chem__EnergyDihedral_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyFixedNonbondRestraint_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyFixedNonbondRestraint_O");
    chem::Register_chem__EnergyFixedNonbondRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyImproperRestraint_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyImproperRestraint_O");
    chem::Register_chem__EnergyImproperRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyNonbond_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyNonbond_O");
    chem::Register_chem__EnergyNonbond_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EnergyStretch_O(core::Lisp_sp); // base(s): set(['chem::EnergyComponent_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EnergyStretch_O");
    chem::Register_chem__EnergyStretch_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EntityNameSetBase_O(core::Lisp_sp); // base(s): set(['chem::Entity_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EntityNameSetBase_O");
    chem::Register_chem__EntityNameSetBase_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFAngleDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFAngleDb_O");
    chem::Register_chem__FFAngleDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFAngle_O(core::Lisp_sp); // base(s): set(['chem::FFParameter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFAngle_O");
    chem::Register_chem__FFAngle_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFItorDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFItorDb_O");
    chem::Register_chem__FFItorDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFItor_O(core::Lisp_sp); // base(s): set(['chem::FFParameter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFItor_O");
    chem::Register_chem__FFItor_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFNonbondDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFNonbondDb_O");
    chem::Register_chem__FFNonbondDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFNonbond_O(core::Lisp_sp); // base(s): set(['chem::FFParameter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFNonbond_O");
    chem::Register_chem__FFNonbond_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFPtorDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFPtorDb_O");
    chem::Register_chem__FFPtorDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFPtor_O(core::Lisp_sp); // base(s): set(['chem::FFParameter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFPtor_O");
    chem::Register_chem__FFPtor_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFStretchDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFStretchDb_O");
    chem::Register_chem__FFStretchDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFStretch_O(core::Lisp_sp); // base(s): set(['chem::FFParameter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFStretch_O");
    chem::Register_chem__FFStretch_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFTypesDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFTypesDb_O");
    chem::Register_chem__FFTypesDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFVdwDb_O(core::Lisp_sp); // base(s): set(['chem::FFBaseDb_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFVdwDb_O");
    chem::Register_chem__FFVdwDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FFVdw_O(core::Lisp_sp); // base(s): set(['chem::FFParameter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FFVdw_O");
    chem::Register_chem__FFVdw_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ForceMatchReport_O(core::Lisp_sp); // base(s): set(['chem::ReportBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ForceMatchReport_O");
    chem::Register_chem__ForceMatchReport_O(_lisp); // base(s): set(['chem::ReportBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__FullLargeSquareMatrix_O(core::Lisp_sp); // base(s): set(['chem::AbstractLargeSquareMatrix_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__FullLargeSquareMatrix_O");
    chem::Register_chem__FullLargeSquareMatrix_O(_lisp); // base(s): set(['chem::AbstractLargeSquareMatrix_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__IterateMatter_O(core::Lisp_sp); // base(s): set(['core::Iterator_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__IterateMatter_O");
    chem::Register_chem__IterateMatter_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__IterateRestraints_O(core::Lisp_sp); // base(s): set(['core::Iterator_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__IterateRestraints_O");
    chem::Register_chem__IterateRestraints_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__JumpPlug_O(core::Lisp_sp); // base(s): set(['chem::Plug_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__JumpPlug_O");
    chem::Register_chem__JumpPlug_O(_lisp); // base(s): set(['chem::Plug_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MessageReport_O(core::Lisp_sp); // base(s): set(['chem::ReportBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MessageReport_O");
    chem::Register_chem__MessageReport_O(_lisp); // base(s): set(['chem::ReportBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Molecule_O(core::Lisp_sp); // base(s): set(['chem::Matter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Molecule_O");
    chem::Register_chem__Molecule_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Monomer_O(core::Lisp_sp); // base(s): set(['chem::Entity_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Monomer_O");
    chem::Register_chem__Monomer_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__OligomerPart_Link_O(core::Lisp_sp); // base(s): set(['chem::OligomerPart_Base_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__OligomerPart_Link_O");
    chem::Register_chem__OligomerPart_Link_O(_lisp); // base(s): set(['chem::OligomerPart_Base_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__OligomerPart_Monomer_O(core::Lisp_sp); // base(s): set(['chem::OligomerPart_Base_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__OligomerPart_Monomer_O");
    chem::Register_chem__OligomerPart_Monomer_O(_lisp); // base(s): set(['chem::OligomerPart_Base_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__PlugWithMates_O(core::Lisp_sp); // base(s): set(['chem::Plug_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__PlugWithMates_O");
    chem::Register_chem__PlugWithMates_O(_lisp); // base(s): set(['chem::Plug_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Residue_O(core::Lisp_sp); // base(s): set(['chem::Matter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Residue_O");
    chem::Register_chem__Residue_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestrainedExoCyclicAtom_O(core::Lisp_sp); // base(s): set(['chem::ComplexRestraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestrainedExoCyclicAtom_O");
    chem::Register_chem__RestrainedExoCyclicAtom_O(_lisp); // base(s): set(['chem::ComplexRestraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestrainedPiBond_O(core::Lisp_sp); // base(s): set(['chem::ComplexRestraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestrainedPiBond_O");
    chem::Register_chem__RestrainedPiBond_O(_lisp); // base(s): set(['chem::ComplexRestraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintAnchor_O(core::Lisp_sp); // base(s): set(['chem::Restraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintAnchor_O");
    chem::Register_chem__RestraintAnchor_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintAngle_O(core::Lisp_sp); // base(s): set(['chem::Restraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintAngle_O");
    chem::Register_chem__RestraintAngle_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintChiral_O(core::Lisp_sp); // base(s): set(['chem::Restraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintChiral_O");
    chem::Register_chem__RestraintChiral_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintDihedral_O(core::Lisp_sp); // base(s): set(['chem::Restraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintDihedral_O");
    chem::Register_chem__RestraintDihedral_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintDistance_O(core::Lisp_sp); // base(s): set(['chem::Restraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintDistance_O");
    chem::Register_chem__RestraintDistance_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RestraintFixedNonbond_O(core::Lisp_sp); // base(s): set(['chem::Restraint_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RestraintFixedNonbond_O");
    chem::Register_chem__RestraintFixedNonbond_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RingCoupling_O(core::Lisp_sp); // base(s): set(['chem::Coupling_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RingCoupling_O");
    chem::Register_chem__RingCoupling_O(_lisp); // base(s): set(['chem::Coupling_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RootMatchNode_O(core::Lisp_sp); // base(s): set(['chem::ChemInfoNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RootMatchNode_O");
    chem::Register_chem__RootMatchNode_O(_lisp); // base(s): set(['chem::ChemInfoNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SparseLargeSquareMatrix_O(core::Lisp_sp); // base(s): set(['chem::AbstractLargeSquareMatrix_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SparseLargeSquareMatrix_O");
    chem::Register_chem__SparseLargeSquareMatrix_O(_lisp); // base(s): set(['chem::AbstractLargeSquareMatrix_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StepReport_O(core::Lisp_sp); // base(s): set(['chem::ReportBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StepReport_O");
    chem::Register_chem__StepReport_O(_lisp); // base(s): set(['chem::ReportBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__StereoisomerVirtualAtom_O(core::Lisp_sp); // base(s): set(['chem::StereoisomerAtom_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__StereoisomerVirtualAtom_O");
    chem::Register_chem__StereoisomerVirtualAtom_O(_lisp); // base(s): set(['chem::StereoisomerAtom_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Stereoisomer_O(core::Lisp_sp); // base(s): set(['chem::Entity_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Stereoisomer_O");
    chem::Register_chem__Stereoisomer_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SuperposableConformationCollection_O(core::Lisp_sp); // base(s): set(['chem::ConformationCollection_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SuperposableConformationCollection_O");
    chem::Register_chem__SuperposableConformationCollection_O(_lisp); // base(s): set(['chem::ConformationCollection_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SuperposeSelectedAtoms_O(core::Lisp_sp); // base(s): set(['chem::SuperposeEngine_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SuperposeSelectedAtoms_O");
    chem::Register_chem__SuperposeSelectedAtoms_O(_lisp); // base(s): set(['chem::SuperposeEngine_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ZMatrixAngleInternal_O(core::Lisp_sp); // base(s): set(['chem::ZMatrixInternal_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ZMatrixAngleInternal_O");
    chem::Register_chem__ZMatrixAngleInternal_O(_lisp); // base(s): set(['chem::ZMatrixInternal_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ZMatrixBondInternal_O(core::Lisp_sp); // base(s): set(['chem::ZMatrixInternal_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ZMatrixBondInternal_O");
    chem::Register_chem__ZMatrixBondInternal_O(_lisp); // base(s): set(['chem::ZMatrixInternal_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ZMatrixDihedralInternal_O(core::Lisp_sp); // base(s): set(['chem::ZMatrixInternal_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ZMatrixDihedralInternal_O");
    chem::Register_chem__ZMatrixDihedralInternal_O(_lisp); // base(s): set(['chem::ZMatrixInternal_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AfterMatchBondTest_O(core::Lisp_sp); // base(s): set(['chem::RootMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AfterMatchBondTest_O");
    chem::Register_chem__AfterMatchBondTest_O(_lisp); // base(s): set(['chem::RootMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AntechamberBondTest_O(core::Lisp_sp); // base(s): set(['chem::BondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AntechamberBondTest_O");
    chem::Register_chem__AntechamberBondTest_O(_lisp); // base(s): set(['chem::BondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomOrBondMatchNode_O(core::Lisp_sp); // base(s): set(['chem::BondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomOrBondMatchNode_O");
    chem::Register_chem__AtomOrBondMatchNode_O(_lisp); // base(s): set(['chem::BondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__BondTest_O(core::Lisp_sp); // base(s): set(['chem::BondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__BondTest_O");
    chem::Register_chem__BondTest_O(_lisp); // base(s): set(['chem::BondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Branch_O(core::Lisp_sp); // base(s): set(['chem::BondListMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Branch_O");
    chem::Register_chem__Branch_O(_lisp); // base(s): set(['chem::BondListMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Chain_O(core::Lisp_sp); // base(s): set(['chem::BondListMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Chain_O");
    chem::Register_chem__Chain_O(_lisp); // base(s): set(['chem::BondListMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EntityNameSetWithCap_O(core::Lisp_sp); // base(s): set(['chem::EntityNameSetBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EntityNameSetWithCap_O");
    chem::Register_chem__EntityNameSetWithCap_O(_lisp); // base(s): set(['chem::EntityNameSetBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__EntityNameSet_O(core::Lisp_sp); // base(s): set(['chem::EntityNameSetBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__EntityNameSet_O");
    chem::Register_chem__EntityNameSet_O(_lisp); // base(s): set(['chem::EntityNameSetBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__InPlug_O(core::Lisp_sp); // base(s): set(['chem::PlugWithMates_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__InPlug_O");
    chem::Register_chem__InPlug_O(_lisp); // base(s): set(['chem::PlugWithMates_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__IterateAtoms_O(core::Lisp_sp); // base(s): set(['chem::IterateMatter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__IterateAtoms_O");
    chem::Register_chem__IterateAtoms_O(_lisp); // base(s): set(['chem::IterateMatter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__IterateBonds_O(core::Lisp_sp); // base(s): set(['chem::IterateMatter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__IterateBonds_O");
    chem::Register_chem__IterateBonds_O(_lisp); // base(s): set(['chem::IterateMatter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__IterateResidues_O(core::Lisp_sp); // base(s): set(['chem::IterateMatter_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__IterateResidues_O");
    chem::Register_chem__IterateResidues_O(_lisp); // base(s): set(['chem::IterateMatter_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MultiMonomer_O(core::Lisp_sp); // base(s): set(['chem::Monomer_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MultiMonomer_O");
    chem::Register_chem__MultiMonomer_O(_lisp); // base(s): set(['chem::Monomer_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__OutPlug_O(core::Lisp_sp); // base(s): set(['chem::PlugWithMates_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__OutPlug_O");
    chem::Register_chem__OutPlug_O(_lisp); // base(s): set(['chem::PlugWithMates_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__VirtualAtom_O(core::Lisp_sp); // base(s): set(['chem::Atom_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__VirtualAtom_O");
    chem::Register_chem__VirtualAtom_O(_lisp); // base(s): set(['chem::Atom_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AntechamberFocusAtomMatch_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AntechamberFocusAtomMatch_O");
    chem::Register_chem__AntechamberFocusAtomMatch_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AtomTest_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AtomTest_O");
    chem::Register_chem__AtomTest_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Logical_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Logical_O");
    chem::Register_chem__Logical_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Mate_O(core::Lisp_sp); // base(s): set(['chem::EntityNameSetWithCap_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Mate_O");
    chem::Register_chem__Mate_O(_lisp); // base(s): set(['chem::EntityNameSetWithCap_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__MonomerPack_O(core::Lisp_sp); // base(s): set(['chem::EntityNameSet_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__MonomerPack_O");
    chem::Register_chem__MonomerPack_O(_lisp); // base(s): set(['chem::EntityNameSet_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RepresentedEntityNameSet_O(core::Lisp_sp); // base(s): set(['chem::EntityNameSet_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RepresentedEntityNameSet_O");
    chem::Register_chem__RepresentedEntityNameSet_O(_lisp); // base(s): set(['chem::EntityNameSet_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__ResidueTest_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__ResidueTest_O");
    chem::Register_chem__ResidueTest_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RingClosingMate_O(core::Lisp_sp); // base(s): set(['chem::EntityNameSetWithCap_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RingClosingMate_O");
    chem::Register_chem__RingClosingMate_O(_lisp); // base(s): set(['chem::EntityNameSetWithCap_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RingClosingPlug_O(core::Lisp_sp); // base(s): set(['chem::OutPlug_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RingClosingPlug_O");
    chem::Register_chem__RingClosingPlug_O(_lisp); // base(s): set(['chem::OutPlug_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__RingTest_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__RingTest_O");
    chem::Register_chem__RingTest_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__Root_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__Root_O");
    chem::Register_chem__Root_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__TagSet_O(core::Lisp_sp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__TagSet_O");
    chem::Register_chem__TagSet_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__AntechamberRoot_O(core::Lisp_sp); // base(s): set(['chem::Root_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__AntechamberRoot_O");
    chem::Register_chem__AntechamberRoot_O(_lisp); // base(s): set(['chem::Root_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
#ifdef EXTERN_REGISTER
extern void Register_chem__SmartsRoot_O(core::Lisp_sp); // base(s): set(['chem::Root_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_chem__SmartsRoot_O");
    chem::Register_chem__SmartsRoot_O(_lisp); // base(s): set(['chem::Root_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ChemPkg
#endif // EXPOSE_TO_CANDO
#undef EXPOSE_TO_CANDO
#ifdef EXPOSE_TO_PYTHON
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AGEdge_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AGEdge_O");
	Call_exposePython_chem__AGEdge_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AGVertex_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AGVertex_O");
	Call_exposePython_chem__AGVertex_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AbstractLargeSquareMatrix_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AbstractLargeSquareMatrix_O");
	Call_exposePython_chem__AbstractLargeSquareMatrix_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Alias_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Alias_O");
	Call_exposePython_chem__Alias_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Angle_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Angle_O");
	Call_exposePython_chem__Angle_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomIdToAtomMap_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomIdToAtomMap_O");
	Call_exposePython_chem__AtomIdToAtomMap_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomId_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomId_O");
	Call_exposePython_chem__AtomId_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomIndexer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomIndexer_O");
	Call_exposePython_chem__AtomIndexer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomReferenceBase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomReferenceBase_O");
	Call_exposePython_chem__AtomReferenceBase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomTable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomTable_O");
	Call_exposePython_chem__AtomTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__BondList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__BondList_O");
	Call_exposePython_chem__BondList_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Bond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Bond_O");
	Call_exposePython_chem__Bond_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CDBond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CDBond_O");
	Call_exposePython_chem__CDBond_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CDFragment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CDFragment_O");
	Call_exposePython_chem__CDFragment_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CDNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CDNode_O");
	Call_exposePython_chem__CDNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CDText_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CDText_O");
	Call_exposePython_chem__CDText_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CalculatePosition_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CalculatePosition_O");
	Call_exposePython_chem__CalculatePosition_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CandoDatabase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CandoDatabase_O");
	Call_exposePython_chem__CandoDatabase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ChemDraw_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ChemDraw_O");
	Call_exposePython_chem__ChemDraw_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ChemInfoMatch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ChemInfoMatch_O");
	Call_exposePython_chem__ChemInfoMatch_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ChemInfoNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ChemInfoNode_O");
	Call_exposePython_chem__ChemInfoNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ChemInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ChemInfo_O");
	Call_exposePython_chem__ChemInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CipPrioritizer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CipPrioritizer_O");
	Call_exposePython_chem__CipPrioritizer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Command_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Command_O");
	Call_exposePython_chem__Command_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ComplexRestraint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ComplexRestraint_O");
	Call_exposePython_chem__ComplexRestraint_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConformationCollectionEntry_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConformationCollectionEntry_O");
	Call_exposePython_chem__ConformationCollectionEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConformationCollection_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConformationCollection_O");
	Call_exposePython_chem__ConformationCollection_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConformationExplorerEntryStage_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConformationExplorerEntryStage_O");
	Call_exposePython_chem__ConformationExplorerEntryStage_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConformationExplorerEntry_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConformationExplorerEntry_O");
	Call_exposePython_chem__ConformationExplorerEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConformationExplorerMatch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConformationExplorerMatch_O");
	Call_exposePython_chem__ConformationExplorerMatch_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConformationExplorer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConformationExplorer_O");
	Call_exposePython_chem__ConformationExplorer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConstitutionAtom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConstitutionAtom_O");
	Call_exposePython_chem__ConstitutionAtom_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConstitutionAtoms_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConstitutionAtoms_O");
	Call_exposePython_chem__ConstitutionAtoms_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConstitutionBond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConstitutionBond_O");
	Call_exposePython_chem__ConstitutionBond_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CoordinateSystem_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CoordinateSystem_O");
	Call_exposePython_chem__CoordinateSystem_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Coupling_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Coupling_O");
	Call_exposePython_chem__Coupling_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyComponent_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyComponent_O");
	Call_exposePython_chem__EnergyComponent_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyFunction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyFunction_O");
	Call_exposePython_chem__EnergyFunction_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Entity_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Entity_O");
	Call_exposePython_chem__Entity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFBaseDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFBaseDb_O");
	Call_exposePython_chem__FFBaseDb_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFNonbondCrossTermTable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFNonbondCrossTermTable_O");
	Call_exposePython_chem__FFNonbondCrossTermTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFParameter_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFParameter_O");
	Call_exposePython_chem__FFParameter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ForceField_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ForceField_O");
	Call_exposePython_chem__ForceField_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FrameRecognizer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FrameRecognizer_O");
	Call_exposePython_chem__FrameRecognizer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ImproperTorsion_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ImproperTorsion_O");
	Call_exposePython_chem__ImproperTorsion_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__InfoDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__InfoDb_O");
	Call_exposePython_chem__InfoDb_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MacroModelFile_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MacroModelFile_O");
	Call_exposePython_chem__MacroModelFile_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MapOfMonomerNamesToAtomIndexers_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MapOfMonomerNamesToAtomIndexers_O");
	Call_exposePython_chem__MapOfMonomerNamesToAtomIndexers_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Matter_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Matter_O");
	Call_exposePython_chem__Matter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MinimizerLog_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MinimizerLog_O");
	Call_exposePython_chem__MinimizerLog_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Minimizer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Minimizer_O");
	Call_exposePython_chem__Minimizer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MonomerContext_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MonomerContext_O");
	Call_exposePython_chem__MonomerContext_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MonomerCoordinates_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MonomerCoordinates_O");
	Call_exposePython_chem__MonomerCoordinates_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__NVector_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__NVector_O");
	Call_exposePython_chem__NVector_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__NumericalFunction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__NumericalFunction_O");
	Call_exposePython_chem__NumericalFunction_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__OligomerPart_Base_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__OligomerPart_Base_O");
	Call_exposePython_chem__OligomerPart_Base_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Oligomer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Oligomer_O");
	Call_exposePython_chem__Oligomer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__OneMonomer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__OneMonomer_O");
	Call_exposePython_chem__OneMonomer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__OneTypeRule_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__OneTypeRule_O");
	Call_exposePython_chem__OneTypeRule_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__PathMessage_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__PathMessage_O");
	Call_exposePython_chem__PathMessage_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__PdbMonomerConnectivity_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__PdbMonomerConnectivity_O");
	Call_exposePython_chem__PdbMonomerConnectivity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__PdbMonomerDatabase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__PdbMonomerDatabase_O");
	Call_exposePython_chem__PdbMonomerDatabase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__PdbReader_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__PdbReader_O");
	Call_exposePython_chem__PdbReader_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__PdbWriter_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__PdbWriter_O");
	Call_exposePython_chem__PdbWriter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Plug_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Plug_O");
	Call_exposePython_chem__Plug_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ProperTorsion_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ProperTorsion_O");
	Call_exposePython_chem__ProperTorsion_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ReadAmberParameters_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ReadAmberParameters_O");
	Call_exposePython_chem__ReadAmberParameters_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ReportBase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ReportBase_O");
	Call_exposePython_chem__ReportBase_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ResidueList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ResidueList_O");
	Call_exposePython_chem__ResidueList_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintList_O");
	Call_exposePython_chem__RestraintList_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Restraint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Restraint_O");
	Call_exposePython_chem__Restraint_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RingFinder_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RingFinder_O");
	Call_exposePython_chem__RingFinder_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Rotamer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Rotamer_O");
	Call_exposePython_chem__Rotamer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SearchStatistics_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SearchStatistics_O");
	Call_exposePython_chem__SearchStatistics_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SpanningLoop_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SpanningLoop_O");
	Call_exposePython_chem__SpanningLoop_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SpecificContextSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SpecificContextSet_O");
	Call_exposePython_chem__SpecificContextSet_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SpecificContext_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SpecificContext_O");
	Call_exposePython_chem__SpecificContext_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StatusTracker_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StatusTracker_O");
	Call_exposePython_chem__StatusTracker_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StereoConfiguration_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StereoConfiguration_O");
	Call_exposePython_chem__StereoConfiguration_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StereoInformation_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StereoInformation_O");
	Call_exposePython_chem__StereoInformation_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StereoisomerAtom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StereoisomerAtom_O");
	Call_exposePython_chem__StereoisomerAtom_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StereoisomerAtoms_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StereoisomerAtoms_O");
	Call_exposePython_chem__StereoisomerAtoms_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StructureComparer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StructureComparer_O");
	Call_exposePython_chem__StructureComparer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Structure_Old_ListEntry_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Structure_Old_ListEntry_O");
	Call_exposePython_chem__Structure_Old_ListEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Structure_Old_List_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Structure_Old_List_O");
	Call_exposePython_chem__Structure_Old_List_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SuperposeEngine_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SuperposeEngine_O");
	Call_exposePython_chem__SuperposeEngine_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Topology_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Topology_O");
	Call_exposePython_chem__Topology_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__TrajectoryFrame_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__TrajectoryFrame_O");
	Call_exposePython_chem__TrajectoryFrame_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Trajectory_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Trajectory_O");
	Call_exposePython_chem__Trajectory_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__TwisterDriver_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__TwisterDriver_O");
	Call_exposePython_chem__TwisterDriver_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Twister_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Twister_O");
	Call_exposePython_chem__Twister_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__TypeAssignmentRules_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__TypeAssignmentRules_O");
	Call_exposePython_chem__TypeAssignmentRules_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__VirtualSphere_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__VirtualSphere_O");
	Call_exposePython_chem__VirtualSphere_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__WildElementDict_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__WildElementDict_O");
	Call_exposePython_chem__WildElementDict_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ZMatrixEntry_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ZMatrixEntry_O");
	Call_exposePython_chem__ZMatrixEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ZMatrixInternal_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ZMatrixInternal_O");
	Call_exposePython_chem__ZMatrixInternal_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ZMatrix_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ZMatrix_O");
	Call_exposePython_chem__ZMatrix_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Aggregate_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Aggregate_O");
	Call_exposePython_chem__Aggregate_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomReference_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomReference_O");
	Call_exposePython_chem__AtomReference_O(_lisp); // base(s): set(['chem::AtomReferenceBase_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Atom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Atom_O");
	Call_exposePython_chem__Atom_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__BondListMatchNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__BondListMatchNode_O");
	Call_exposePython_chem__BondListMatchNode_O(_lisp); // base(s): set(['chem::ChemInfoNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__BondMatchNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__BondMatchNode_O");
	Call_exposePython_chem__BondMatchNode_O(_lisp); // base(s): set(['chem::ChemInfoNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CalculatePositionAlongBond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CalculatePositionAlongBond_O");
	Call_exposePython_chem__CalculatePositionAlongBond_O(_lisp); // base(s): set(['chem::CalculatePosition_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CalculatePositionRelativeToOrigin_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CalculatePositionRelativeToOrigin_O");
	Call_exposePython_chem__CalculatePositionRelativeToOrigin_O(_lisp); // base(s): set(['chem::CalculatePosition_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__CalculatePositionUsingInternals_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__CalculatePositionUsingInternals_O");
	Call_exposePython_chem__CalculatePositionUsingInternals_O(_lisp); // base(s): set(['chem::CalculatePosition_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ConstitutionVirtualAtom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ConstitutionVirtualAtom_O");
	Call_exposePython_chem__ConstitutionVirtualAtom_O(_lisp); // base(s): set(['chem::ConstitutionAtom_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Constitution_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Constitution_O");
	Call_exposePython_chem__Constitution_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__DirectionalCoupling_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__DirectionalCoupling_O");
	Call_exposePython_chem__DirectionalCoupling_O(_lisp); // base(s): set(['chem::Coupling_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyAnchorRestraint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyAnchorRestraint_O");
	Call_exposePython_chem__EnergyAnchorRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyAngle_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyAngle_O");
	Call_exposePython_chem__EnergyAngle_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyChiralRestraint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyChiralRestraint_O");
	Call_exposePython_chem__EnergyChiralRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyDihedral_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyDihedral_O");
	Call_exposePython_chem__EnergyDihedral_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyFixedNonbondRestraint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyFixedNonbondRestraint_O");
	Call_exposePython_chem__EnergyFixedNonbondRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyImproperRestraint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyImproperRestraint_O");
	Call_exposePython_chem__EnergyImproperRestraint_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyNonbond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyNonbond_O");
	Call_exposePython_chem__EnergyNonbond_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EnergyStretch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EnergyStretch_O");
	Call_exposePython_chem__EnergyStretch_O(_lisp); // base(s): set(['chem::EnergyComponent_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EntityNameSetBase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EntityNameSetBase_O");
	Call_exposePython_chem__EntityNameSetBase_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFAngleDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFAngleDb_O");
	Call_exposePython_chem__FFAngleDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFAngle_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFAngle_O");
	Call_exposePython_chem__FFAngle_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFItorDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFItorDb_O");
	Call_exposePython_chem__FFItorDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFItor_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFItor_O");
	Call_exposePython_chem__FFItor_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFNonbondDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFNonbondDb_O");
	Call_exposePython_chem__FFNonbondDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFNonbond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFNonbond_O");
	Call_exposePython_chem__FFNonbond_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFPtorDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFPtorDb_O");
	Call_exposePython_chem__FFPtorDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFPtor_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFPtor_O");
	Call_exposePython_chem__FFPtor_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFStretchDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFStretchDb_O");
	Call_exposePython_chem__FFStretchDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFStretch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFStretch_O");
	Call_exposePython_chem__FFStretch_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFTypesDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFTypesDb_O");
	Call_exposePython_chem__FFTypesDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFVdwDb_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFVdwDb_O");
	Call_exposePython_chem__FFVdwDb_O(_lisp); // base(s): set(['chem::FFBaseDb_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FFVdw_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FFVdw_O");
	Call_exposePython_chem__FFVdw_O(_lisp); // base(s): set(['chem::FFParameter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ForceMatchReport_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ForceMatchReport_O");
	Call_exposePython_chem__ForceMatchReport_O(_lisp); // base(s): set(['chem::ReportBase_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__FullLargeSquareMatrix_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__FullLargeSquareMatrix_O");
	Call_exposePython_chem__FullLargeSquareMatrix_O(_lisp); // base(s): set(['chem::AbstractLargeSquareMatrix_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__IterateMatter_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__IterateMatter_O");
	Call_exposePython_chem__IterateMatter_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__IterateRestraints_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__IterateRestraints_O");
	Call_exposePython_chem__IterateRestraints_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__JumpPlug_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__JumpPlug_O");
	Call_exposePython_chem__JumpPlug_O(_lisp); // base(s): set(['chem::Plug_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MessageReport_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MessageReport_O");
	Call_exposePython_chem__MessageReport_O(_lisp); // base(s): set(['chem::ReportBase_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Molecule_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Molecule_O");
	Call_exposePython_chem__Molecule_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Monomer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Monomer_O");
	Call_exposePython_chem__Monomer_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__OligomerPart_Link_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__OligomerPart_Link_O");
	Call_exposePython_chem__OligomerPart_Link_O(_lisp); // base(s): set(['chem::OligomerPart_Base_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__OligomerPart_Monomer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__OligomerPart_Monomer_O");
	Call_exposePython_chem__OligomerPart_Monomer_O(_lisp); // base(s): set(['chem::OligomerPart_Base_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__PlugWithMates_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__PlugWithMates_O");
	Call_exposePython_chem__PlugWithMates_O(_lisp); // base(s): set(['chem::Plug_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Residue_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Residue_O");
	Call_exposePython_chem__Residue_O(_lisp); // base(s): set(['chem::Matter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestrainedExoCyclicAtom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestrainedExoCyclicAtom_O");
	Call_exposePython_chem__RestrainedExoCyclicAtom_O(_lisp); // base(s): set(['chem::ComplexRestraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestrainedPiBond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestrainedPiBond_O");
	Call_exposePython_chem__RestrainedPiBond_O(_lisp); // base(s): set(['chem::ComplexRestraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintAnchor_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintAnchor_O");
	Call_exposePython_chem__RestraintAnchor_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintAngle_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintAngle_O");
	Call_exposePython_chem__RestraintAngle_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintChiral_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintChiral_O");
	Call_exposePython_chem__RestraintChiral_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintDihedral_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintDihedral_O");
	Call_exposePython_chem__RestraintDihedral_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintDistance_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintDistance_O");
	Call_exposePython_chem__RestraintDistance_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RestraintFixedNonbond_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RestraintFixedNonbond_O");
	Call_exposePython_chem__RestraintFixedNonbond_O(_lisp); // base(s): set(['chem::Restraint_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RingCoupling_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RingCoupling_O");
	Call_exposePython_chem__RingCoupling_O(_lisp); // base(s): set(['chem::Coupling_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RootMatchNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RootMatchNode_O");
	Call_exposePython_chem__RootMatchNode_O(_lisp); // base(s): set(['chem::ChemInfoNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SparseLargeSquareMatrix_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SparseLargeSquareMatrix_O");
	Call_exposePython_chem__SparseLargeSquareMatrix_O(_lisp); // base(s): set(['chem::AbstractLargeSquareMatrix_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StepReport_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StepReport_O");
	Call_exposePython_chem__StepReport_O(_lisp); // base(s): set(['chem::ReportBase_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__StereoisomerVirtualAtom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__StereoisomerVirtualAtom_O");
	Call_exposePython_chem__StereoisomerVirtualAtom_O(_lisp); // base(s): set(['chem::StereoisomerAtom_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Stereoisomer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Stereoisomer_O");
	Call_exposePython_chem__Stereoisomer_O(_lisp); // base(s): set(['chem::Entity_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SuperposableConformationCollection_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SuperposableConformationCollection_O");
	Call_exposePython_chem__SuperposableConformationCollection_O(_lisp); // base(s): set(['chem::ConformationCollection_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SuperposeSelectedAtoms_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SuperposeSelectedAtoms_O");
	Call_exposePython_chem__SuperposeSelectedAtoms_O(_lisp); // base(s): set(['chem::SuperposeEngine_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ZMatrixAngleInternal_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ZMatrixAngleInternal_O");
	Call_exposePython_chem__ZMatrixAngleInternal_O(_lisp); // base(s): set(['chem::ZMatrixInternal_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ZMatrixBondInternal_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ZMatrixBondInternal_O");
	Call_exposePython_chem__ZMatrixBondInternal_O(_lisp); // base(s): set(['chem::ZMatrixInternal_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ZMatrixDihedralInternal_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ZMatrixDihedralInternal_O");
	Call_exposePython_chem__ZMatrixDihedralInternal_O(_lisp); // base(s): set(['chem::ZMatrixInternal_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AfterMatchBondTest_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AfterMatchBondTest_O");
	Call_exposePython_chem__AfterMatchBondTest_O(_lisp); // base(s): set(['chem::RootMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AntechamberBondTest_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AntechamberBondTest_O");
	Call_exposePython_chem__AntechamberBondTest_O(_lisp); // base(s): set(['chem::BondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomOrBondMatchNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomOrBondMatchNode_O");
	Call_exposePython_chem__AtomOrBondMatchNode_O(_lisp); // base(s): set(['chem::BondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__BondTest_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__BondTest_O");
	Call_exposePython_chem__BondTest_O(_lisp); // base(s): set(['chem::BondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Branch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Branch_O");
	Call_exposePython_chem__Branch_O(_lisp); // base(s): set(['chem::BondListMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Chain_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Chain_O");
	Call_exposePython_chem__Chain_O(_lisp); // base(s): set(['chem::BondListMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EntityNameSetWithCap_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EntityNameSetWithCap_O");
	Call_exposePython_chem__EntityNameSetWithCap_O(_lisp); // base(s): set(['chem::EntityNameSetBase_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__EntityNameSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__EntityNameSet_O");
	Call_exposePython_chem__EntityNameSet_O(_lisp); // base(s): set(['chem::EntityNameSetBase_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__InPlug_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__InPlug_O");
	Call_exposePython_chem__InPlug_O(_lisp); // base(s): set(['chem::PlugWithMates_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__IterateAtoms_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__IterateAtoms_O");
	Call_exposePython_chem__IterateAtoms_O(_lisp); // base(s): set(['chem::IterateMatter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__IterateBonds_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__IterateBonds_O");
	Call_exposePython_chem__IterateBonds_O(_lisp); // base(s): set(['chem::IterateMatter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__IterateResidues_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__IterateResidues_O");
	Call_exposePython_chem__IterateResidues_O(_lisp); // base(s): set(['chem::IterateMatter_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MultiMonomer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MultiMonomer_O");
	Call_exposePython_chem__MultiMonomer_O(_lisp); // base(s): set(['chem::Monomer_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__OutPlug_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__OutPlug_O");
	Call_exposePython_chem__OutPlug_O(_lisp); // base(s): set(['chem::PlugWithMates_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__VirtualAtom_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__VirtualAtom_O");
	Call_exposePython_chem__VirtualAtom_O(_lisp); // base(s): set(['chem::Atom_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AntechamberFocusAtomMatch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AntechamberFocusAtomMatch_O");
	Call_exposePython_chem__AntechamberFocusAtomMatch_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AtomTest_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AtomTest_O");
	Call_exposePython_chem__AtomTest_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Logical_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Logical_O");
	Call_exposePython_chem__Logical_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Mate_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Mate_O");
	Call_exposePython_chem__Mate_O(_lisp); // base(s): set(['chem::EntityNameSetWithCap_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__MonomerPack_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__MonomerPack_O");
	Call_exposePython_chem__MonomerPack_O(_lisp); // base(s): set(['chem::EntityNameSet_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RepresentedEntityNameSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RepresentedEntityNameSet_O");
	Call_exposePython_chem__RepresentedEntityNameSet_O(_lisp); // base(s): set(['chem::EntityNameSet_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__ResidueTest_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__ResidueTest_O");
	Call_exposePython_chem__ResidueTest_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RingClosingMate_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RingClosingMate_O");
	Call_exposePython_chem__RingClosingMate_O(_lisp); // base(s): set(['chem::EntityNameSetWithCap_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RingClosingPlug_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RingClosingPlug_O");
	Call_exposePython_chem__RingClosingPlug_O(_lisp); // base(s): set(['chem::OutPlug_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__RingTest_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__RingTest_O");
	Call_exposePython_chem__RingTest_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__Root_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__Root_O");
	Call_exposePython_chem__Root_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__TagSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__TagSet_O");
	Call_exposePython_chem__TagSet_O(_lisp); // base(s): set(['chem::AtomOrBondMatchNode_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__AntechamberRoot_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__AntechamberRoot_O");
	Call_exposePython_chem__AntechamberRoot_O(_lisp); // base(s): set(['chem::Root_O'])
}
#endif // ifdef Use_ChemPkg
#ifdef Use_ChemPkg
extern void Call_exposePython_chem__SmartsRoot_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: chem__SmartsRoot_O");
	Call_exposePython_chem__SmartsRoot_O(_lisp); // base(s): set(['chem::Root_O'])
}
#endif // ifdef Use_ChemPkg
#endif // EXPOSE_TO_PYTHON
#undef EXPOSE_TO_PYTHON
#if defined(EXPAND_CLASS_MACROS)
_CLASS_MACRO(chem::AGEdge_O)
_CLASS_MACRO(chem::AGVertex_O)
_CLASS_MACRO(chem::AbstractLargeSquareMatrix_O)
_CLASS_MACRO(chem::Alias_O)
_CLASS_MACRO(chem::Angle_O)
_CLASS_MACRO(chem::AtomIdToAtomMap_O)
_CLASS_MACRO(chem::AtomId_O)
_CLASS_MACRO(chem::AtomIndexer_O)
_CLASS_MACRO(chem::AtomReferenceBase_O)
_CLASS_MACRO(chem::AtomTable_O)
_CLASS_MACRO(chem::BondList_O)
_CLASS_MACRO(chem::Bond_O)
_CLASS_MACRO(chem::CDBond_O)
_CLASS_MACRO(chem::CDFragment_O)
_CLASS_MACRO(chem::CDNode_O)
_CLASS_MACRO(chem::CDText_O)
_CLASS_MACRO(chem::CalculatePosition_O)
_CLASS_MACRO(chem::CandoDatabase_O)
_CLASS_MACRO(chem::ChemDraw_O)
_CLASS_MACRO(chem::ChemInfoMatch_O)
_CLASS_MACRO(chem::ChemInfoNode_O)
_CLASS_MACRO(chem::ChemInfo_O)
_CLASS_MACRO(chem::CipPrioritizer_O)
_CLASS_MACRO(chem::Command_O)
_CLASS_MACRO(chem::ComplexRestraint_O)
_CLASS_MACRO(chem::ConformationCollectionEntry_O)
_CLASS_MACRO(chem::ConformationCollection_O)
_CLASS_MACRO(chem::ConformationExplorerEntryStage_O)
_CLASS_MACRO(chem::ConformationExplorerEntry_O)
_CLASS_MACRO(chem::ConformationExplorerMatch_O)
_CLASS_MACRO(chem::ConformationExplorer_O)
_CLASS_MACRO(chem::ConstitutionAtom_O)
_CLASS_MACRO(chem::ConstitutionAtoms_O)
_CLASS_MACRO(chem::ConstitutionBond_O)
_CLASS_MACRO(chem::CoordinateSystem_O)
_CLASS_MACRO(chem::Coupling_O)
_CLASS_MACRO(chem::EnergyComponent_O)
_CLASS_MACRO(chem::EnergyFunction_O)
_CLASS_MACRO(chem::Entity_O)
_CLASS_MACRO(chem::FFBaseDb_O)
_CLASS_MACRO(chem::FFNonbondCrossTermTable_O)
_CLASS_MACRO(chem::FFParameter_O)
_CLASS_MACRO(chem::ForceField_O)
_CLASS_MACRO(chem::FrameRecognizer_O)
_CLASS_MACRO(chem::ImproperTorsion_O)
_CLASS_MACRO(chem::InfoDb_O)
_CLASS_MACRO(chem::MacroModelFile_O)
_CLASS_MACRO(chem::MapOfMonomerNamesToAtomIndexers_O)
_CLASS_MACRO(chem::Matter_O)
_CLASS_MACRO(chem::MinimizerLog_O)
_CLASS_MACRO(chem::Minimizer_O)
_CLASS_MACRO(chem::MonomerContext_O)
_CLASS_MACRO(chem::MonomerCoordinates_O)
_CLASS_MACRO(chem::NVector_O)
_CLASS_MACRO(chem::NumericalFunction_O)
_CLASS_MACRO(chem::OligomerPart_Base_O)
_CLASS_MACRO(chem::Oligomer_O)
_CLASS_MACRO(chem::OneMonomer_O)
_CLASS_MACRO(chem::OneTypeRule_O)
_CLASS_MACRO(chem::PathMessage_O)
_CLASS_MACRO(chem::PdbMonomerConnectivity_O)
_CLASS_MACRO(chem::PdbMonomerDatabase_O)
_CLASS_MACRO(chem::PdbReader_O)
_CLASS_MACRO(chem::PdbWriter_O)
_CLASS_MACRO(chem::Plug_O)
_CLASS_MACRO(chem::ProperTorsion_O)
_CLASS_MACRO(chem::ReadAmberParameters_O)
_CLASS_MACRO(chem::ReportBase_O)
_CLASS_MACRO(chem::ResidueList_O)
_CLASS_MACRO(chem::RestraintList_O)
_CLASS_MACRO(chem::Restraint_O)
_CLASS_MACRO(chem::RingFinder_O)
_CLASS_MACRO(chem::Rotamer_O)
_CLASS_MACRO(chem::SearchStatistics_O)
_CLASS_MACRO(chem::SpanningLoop_O)
_CLASS_MACRO(chem::SpecificContextSet_O)
_CLASS_MACRO(chem::SpecificContext_O)
_CLASS_MACRO(chem::StatusTracker_O)
_CLASS_MACRO(chem::StereoConfiguration_O)
_CLASS_MACRO(chem::StereoInformation_O)
_CLASS_MACRO(chem::StereoisomerAtom_O)
_CLASS_MACRO(chem::StereoisomerAtoms_O)
_CLASS_MACRO(chem::StructureComparer_O)
_CLASS_MACRO(chem::Structure_Old_ListEntry_O)
_CLASS_MACRO(chem::Structure_Old_List_O)
_CLASS_MACRO(chem::SuperposeEngine_O)
_CLASS_MACRO(chem::Topology_O)
_CLASS_MACRO(chem::TrajectoryFrame_O)
_CLASS_MACRO(chem::Trajectory_O)
_CLASS_MACRO(chem::TwisterDriver_O)
_CLASS_MACRO(chem::Twister_O)
_CLASS_MACRO(chem::TypeAssignmentRules_O)
_CLASS_MACRO(chem::VirtualSphere_O)
_CLASS_MACRO(chem::WildElementDict_O)
_CLASS_MACRO(chem::ZMatrixEntry_O)
_CLASS_MACRO(chem::ZMatrixInternal_O)
_CLASS_MACRO(chem::ZMatrix_O)
_CLASS_MACRO(chem::Aggregate_O)
_CLASS_MACRO(chem::AtomReference_O)
_CLASS_MACRO(chem::Atom_O)
_CLASS_MACRO(chem::BondListMatchNode_O)
_CLASS_MACRO(chem::BondMatchNode_O)
_CLASS_MACRO(chem::CalculatePositionAlongBond_O)
_CLASS_MACRO(chem::CalculatePositionRelativeToOrigin_O)
_CLASS_MACRO(chem::CalculatePositionUsingInternals_O)
_CLASS_MACRO(chem::ConstitutionVirtualAtom_O)
_CLASS_MACRO(chem::Constitution_O)
_CLASS_MACRO(chem::DirectionalCoupling_O)
_CLASS_MACRO(chem::EnergyAnchorRestraint_O)
_CLASS_MACRO(chem::EnergyAngle_O)
_CLASS_MACRO(chem::EnergyChiralRestraint_O)
_CLASS_MACRO(chem::EnergyDihedral_O)
_CLASS_MACRO(chem::EnergyFixedNonbondRestraint_O)
_CLASS_MACRO(chem::EnergyImproperRestraint_O)
_CLASS_MACRO(chem::EnergyNonbond_O)
_CLASS_MACRO(chem::EnergyStretch_O)
_CLASS_MACRO(chem::EntityNameSetBase_O)
_CLASS_MACRO(chem::FFAngleDb_O)
_CLASS_MACRO(chem::FFAngle_O)
_CLASS_MACRO(chem::FFItorDb_O)
_CLASS_MACRO(chem::FFItor_O)
_CLASS_MACRO(chem::FFNonbondDb_O)
_CLASS_MACRO(chem::FFNonbond_O)
_CLASS_MACRO(chem::FFPtorDb_O)
_CLASS_MACRO(chem::FFPtor_O)
_CLASS_MACRO(chem::FFStretchDb_O)
_CLASS_MACRO(chem::FFStretch_O)
_CLASS_MACRO(chem::FFTypesDb_O)
_CLASS_MACRO(chem::FFVdwDb_O)
_CLASS_MACRO(chem::FFVdw_O)
_CLASS_MACRO(chem::ForceMatchReport_O)
_CLASS_MACRO(chem::FullLargeSquareMatrix_O)
_CLASS_MACRO(chem::IterateMatter_O)
_CLASS_MACRO(chem::IterateRestraints_O)
_CLASS_MACRO(chem::JumpPlug_O)
_CLASS_MACRO(chem::MessageReport_O)
_CLASS_MACRO(chem::Molecule_O)
_CLASS_MACRO(chem::Monomer_O)
_CLASS_MACRO(chem::OligomerPart_Link_O)
_CLASS_MACRO(chem::OligomerPart_Monomer_O)
_CLASS_MACRO(chem::PlugWithMates_O)
_CLASS_MACRO(chem::Residue_O)
_CLASS_MACRO(chem::RestrainedExoCyclicAtom_O)
_CLASS_MACRO(chem::RestrainedPiBond_O)
_CLASS_MACRO(chem::RestraintAnchor_O)
_CLASS_MACRO(chem::RestraintAngle_O)
_CLASS_MACRO(chem::RestraintChiral_O)
_CLASS_MACRO(chem::RestraintDihedral_O)
_CLASS_MACRO(chem::RestraintDistance_O)
_CLASS_MACRO(chem::RestraintFixedNonbond_O)
_CLASS_MACRO(chem::RingCoupling_O)
_CLASS_MACRO(chem::RootMatchNode_O)
_CLASS_MACRO(chem::SparseLargeSquareMatrix_O)
_CLASS_MACRO(chem::StepReport_O)
_CLASS_MACRO(chem::StereoisomerVirtualAtom_O)
_CLASS_MACRO(chem::Stereoisomer_O)
_CLASS_MACRO(chem::SuperposableConformationCollection_O)
_CLASS_MACRO(chem::SuperposeSelectedAtoms_O)
_CLASS_MACRO(chem::ZMatrixAngleInternal_O)
_CLASS_MACRO(chem::ZMatrixBondInternal_O)
_CLASS_MACRO(chem::ZMatrixDihedralInternal_O)
_CLASS_MACRO(chem::AfterMatchBondTest_O)
_CLASS_MACRO(chem::AntechamberBondTest_O)
_CLASS_MACRO(chem::AtomOrBondMatchNode_O)
_CLASS_MACRO(chem::BondTest_O)
_CLASS_MACRO(chem::Branch_O)
_CLASS_MACRO(chem::Chain_O)
_CLASS_MACRO(chem::EntityNameSetWithCap_O)
_CLASS_MACRO(chem::EntityNameSet_O)
_CLASS_MACRO(chem::InPlug_O)
_CLASS_MACRO(chem::IterateAtoms_O)
_CLASS_MACRO(chem::IterateBonds_O)
_CLASS_MACRO(chem::IterateResidues_O)
_CLASS_MACRO(chem::MultiMonomer_O)
_CLASS_MACRO(chem::OutPlug_O)
_CLASS_MACRO(chem::VirtualAtom_O)
_CLASS_MACRO(chem::AntechamberFocusAtomMatch_O)
_CLASS_MACRO(chem::AtomTest_O)
_CLASS_MACRO(chem::Logical_O)
_CLASS_MACRO(chem::Mate_O)
_CLASS_MACRO(chem::MonomerPack_O)
_CLASS_MACRO(chem::RepresentedEntityNameSet_O)
_CLASS_MACRO(chem::ResidueTest_O)
_CLASS_MACRO(chem::RingClosingMate_O)
_CLASS_MACRO(chem::RingClosingPlug_O)
_CLASS_MACRO(chem::RingTest_O)
_CLASS_MACRO(chem::Root_O)
_CLASS_MACRO(chem::TagSet_O)
_CLASS_MACRO(chem::AntechamberRoot_O)
_CLASS_MACRO(chem::SmartsRoot_O)
#endif // EXPAND_CLASS_MACROS
#undef ALL_STAGES
// ---------------- after class initializers
#ifdef EXPOSE_TO_PYTHON
#endif // EXPOSE_TO_PYTHON
