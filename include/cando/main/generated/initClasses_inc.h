// start
// define cpp macros: SET_SYMBOL, CREATE_CLASS, SET_CLASS, DEFINE_BASE_CLASSES, DEFINE_CLASS_NAMES, EXPOSE_TO_CANDO 
// define cpp macro: ALL_STAGES to get the effect of defining all of the macros above
// define cpp macro: EXPOSE_PYTHON to expose python
 // class ActivationFrame_O : public Environment_O // RuntimeVisibleEnvironment_O
 // class ValueFrame_O : public ActivationFrame_O {
 // class FunctionFrame_O : public ActivationFrame_O {
 // class TagbodyFrame_O : public ActivationFrame_O {
 // class Array_O : public T_O {
 // class ArrayObjects_O : public Array_O {
 // class Bignum_O : public Integer_O {
 // class Binder_O : public T_O {
 // class BitVector_O : public Vector_O {
 // class SimpleBitVector_O : public BitVector_O {
 // class BuiltInClass_O : public Class_O {
 // class Character_dummy_O : public T_O {
 // class CxxClass_O : public Class_O {
 // class CandoException_O : public T_O {
 // class Cons_O : public T_O {
 // class CxxObject_O : public T_O // StandardObject_O
 // class Environment_O : public T_O {
 // class LexicalEnvironment_O : public Environment_O {
 // class RuntimeVisibleEnvironment_O : public LexicalEnvironment_O {
 // class ValueEnvironment_O : public RuntimeVisibleEnvironment_O {
 // class FunctionValueEnvironment_O : public RuntimeVisibleEnvironment_O {
 // class CompileTimeEnvironment_O : public LexicalEnvironment_O {
 // class UnwindProtectEnvironment_O : public CompileTimeEnvironment_O {
 // class BlockEnvironment_O : public CompileTimeEnvironment_O {
 // class CatchEnvironment_O : public CompileTimeEnvironment_O {
 // class FunctionContainerEnvironment_O : public CompileTimeEnvironment_O {
 // class TagbodyEnvironment_O : public RuntimeVisibleEnvironment_O {
 // class MacroletEnvironment_O : public CompileTimeEnvironment_O {
 // class SymbolMacroletEnvironment_O : public CompileTimeEnvironment_O {
 // class StackValueEnvironment_O : public CompileTimeEnvironment_O {
 // class GlueEnvironment_O : public Environment_O {
 // class Function_O : public T_O {
 // class CompiledFunction_O : public Function_O {
 // class ExternalObject_O : public T_O // StandardObject_O
 // class ForeignData_O : public ExternalObject_O // StandardObject_O
 // class Path_O : public T_O {
 // class DirectoryIterator_O : public Iterator_O {
 // class RecursiveDirectoryIterator_O : public Iterator_O {
 // class DirectoryEntry_O : public T_O {
 // class FileStatus_O : public T_O {
 // class ForwardReferencedClass_O : public Class_O {
 // class FuncallableStandardClass_O : public StdClass_O {
 // class HashTable_O : public T_O {
 // class HashTableEq_O : public HashTable_O {
 // class HashTableEql_O : public HashTable_O {
 // class HashTableEqual_O : public HashTable_O {
 // class HashTableEqualp_O : public HashTable_O {
 // class Instance_O : public Function_O {
 // class IntArray_O : public T_O {
 // class Iterator_O : public T_O {
 // class LambdaListHandler_O : public T_O {
 // class VaList_dummy_O : public T_O {
 // class Stream_O : public T_O {
 // class AnsiStream_O : public Stream_O {
 // class FileStream_O : public AnsiStream_O {
 // class IOFileStream_O : public FileStream_O {
 // class IOStreamStream_O : public FileStream_O {
 // class StringStream_O : public AnsiStream_O {
 // class StringOutputStream_O : public StringStream_O {
 // class StringInputStream_O : public StringStream_O {
 // class SynonymStream_O : public AnsiStream_O {
 // class TwoWayStream_O : public AnsiStream_O {
 // class BroadcastStream_O : public AnsiStream_O {
 // class ConcatenatedStream_O : public AnsiStream_O {
 // class EchoStream_O : public AnsiStream_O {
 // class String_O : public Vector_O {
 // class Vector_O : public Array_O {
 // class LoadTimeValues_O : public T_O {
 // class Class_O : public Specializer_O {
 // class Metaobject_O : public StandardObject_O {
 // class VectorObjects_O : public Vector_O {
 // class VectorObjectsWithFillPtr_O : public VectorObjects_O {
 // class Null_O : public Symbol_O {
 // class Number_O : public T_O {
 // class Real_O : public Number_O {
 // class Rational_O : public Real_O {
 // class Integer_O : public Rational_O {
 // class Fixnum_dummy_O : public Integer_O {
 // class Float_O : public Real_O {
 // class ShortFloat_O : public Float_O {
 // class SingleFloat_dummy_O : public Float_O {
 // class DoubleFloat_O : public Float_O {
 // class LongFloat_O : public Float_O {
 // class Complex_O : public Number_O {
 // class Ratio_O : public Rational_O {
 // class T_O : public _RootDummyClass {
 // class Pathname_O : public T_O {
 // class LogicalPathname_O : public Pathname_O {
 // class SmallMultimap_O : public T_O {
 // class VectorObjectsWithFillPtr_O : public VectorObjects_O {
 // class Pointer_O : public T_O {
 // class PosixTime_O : public T_O {
 // class PosixTimeDuration_O : public T_O {
 // class InvocationHistoryFrameIterator_O : public T_O {
 // class Reader_O : public T_O {
 // class ReadTable_O : public T_O {
 // class Record_O : public T_O {
 // class RegexMatch_O : public core::T_O {
 // class Regex_O : public core::T_O {
 // class SNode_O : public T_O {
 // class LeafSNode_O : public SNode_O {
 // class BranchSNode_O : public SNode_O {
 // class Archive_O : public T_O {
 // class LoadArchive_O : public Archive_O {
 // class SaveArchive_O : public Archive_O {
 // class SexpLoadArchive_O : public LoadArchive_O {
 // class SexpSaveArchive_O : public SaveArchive_O {
 // class SingleDispatchEffectiveMethodFunction_O : public Function_O {
 // class SingleDispatchGenericFunction_O : public Function_O {
 // class SingleDispatchMethod_O : public T_O {
 // class SmallMap_O : public T_O {
 // class SourceFileInfo_O : public T_O {
 // class SourcePosInfo_O : public T_O {
 // class SourceManager_O : public T_O {
 // class SpecialForm_O : public Function_O {
 // class Specializer_O : public Metaobject_O {
 // class StandardClass_O : public StdClass_O {
 // class StandardObject_O : public T_O {
 // class StdClass_O : public Class_O {
 // class Str_O : public String_O {
 // class StrWithFillPtr_O : public Str_O {
 // class StructureClass_O : public Class_O {
 // class StructureObject_O : public T_O {
 // class LightUserData_O : public core::T_O // StandardObject_O
 // class UserData_O : public core::LightUserData_O // StandardObject_O
 // class SymbolToEnumConverter_O : public T_O {
 // class Symbol_O : public T_O {
 // class VectorObjects_O : public Vector_O {
 // class WeakHashTable_O : public T_O {
 // class WeakKeyHashTable_O : public WeakHashTable_O {
 // class WeakKeyMapping_O : public T_O {
 // class WeakPointer_O : public T_O {
 // class WrappedPointer_O : public core::T_O {
 // class ArrayDisplaced_O : public Array_O {
 // class MultiStringBuffer_O : public T_O {
 // class RandomState_O : public T_O {
 // class Package_O : public T_O {
 // class ClassRegistry_O : public core::T_O {
 // class ClassRep_O : public core::BuiltInClass_O {
 // class Pointer_O : public core::T_O {
 //class AstVisitor_O : public core::T_O, public clang::RecursiveASTVisitor<AstVisitor_O> {
 // class DebugInfo_O : public core::T_O {
 //class DIDescriptor_O : public DebugInfo_O, public llvm::DIDescriptor {
 //class DIScope_O : public DebugInfo_O, public llvm::DIScope {
 //class DIArray_O : public DebugInfo_O, public llvm::DIArray {
 //class DITypeArray_O : public DebugInfo_O, public llvm::DITypeArray {
 //class DIFile_O : public DebugInfo_O, public llvm::DIFile {
 //class DISubprogram_O : public DebugInfo_O, public llvm::DISubprogram {
 //class DIType_O : public DebugInfo_O, public llvm::DIType {
 //class DIDerivedType_O : public DebugInfo_O, public llvm::DIDerivedType {
 //class DIBasicType_O : public DebugInfo_O, public llvm::DIBasicType {
 //class DICompositeType_O : public DebugInfo_O, public llvm::DICompositeType {
 //class DISubroutineType_O : public DebugInfo_O, public llvm::DISubroutineType {
 //class DILexicalBlock_O : public DebugInfo_O, public llvm::DILexicalBlock {
 //class DICompileUnit_O : public DebugInfo_O, public llvm::DICompileUnit {
 // class DIBuilder_O : public core::ExternalObject_O {
 // class DebugLoc_O : public core::T_O {
 // class InsertPoint_O : public core::T_O {
 // class LLVMContext_O : public core::ExternalObject_O {
 // class Linker_O : public core::ExternalObject_O {
 // class Pass_O : public core::ExternalObject_O {
 // class AttributeSet_O : public core::T_O {
 // class Triple_O : public core::ExternalObject_O {
 // class TargetOptions_O : public core::ExternalObject_O {
 // class Target_O : public core::ExternalObject_O {
 // class MCSubtargetInfo_O : public core::ExternalObject_O {
 // class TargetSubtargetInfo_O : public MCSubtargetInfo_O {
 // class TargetMachine_O : public core::ExternalObject_O {
 // class LLVMTargetMachine_O : public TargetMachine_O {
 // class FunctionPass_O : public Pass_O {
 // class ModulePass_O : public Pass_O {
 // class ImmutablePass_O : public ModulePass_O {
 // class PassManagerBase_O : public core::ExternalObject_O {
 // class Value_O : public core::ExternalObject_O {
 // class Metadata_O : public core::ExternalObject_O {
 // class User_O : public Value_O {
 // class Attribute_O : public core::T_O {
 // class DataLayout_O : public core::ExternalObject_O {
 // class Constant_O : public User_O {
 // class ConstantArray_O : public Constant_O {
 // class BlockAddress_O : public Constant_O {
 // class ConstantDataSequential_O : public Constant_O {
 // class ConstantDataArray_O : public ConstantDataSequential_O {
 // class ConstantExpr_O : public Constant_O {
 // class GlobalValue_O : public Constant_O {
 // class GlobalVariable_O : public GlobalValue_O {
 // class ExecutionEngine_O : public core::ExternalObject_O {
 // class Module_O : public core::ExternalObject_O {
 // class DataLayoutPass_O : public ImmutablePass_O {
 // class TargetLibraryInfo_O : public ImmutablePass_O {
 // class FunctionPassManager_O : public PassManagerBase_O {
 // class PassManager_O : public PassManagerBase_O {
 // class EngineBuilder_O : public core::ExternalObject_O {
 // class PassManagerBuilder_O : public core::ExternalObject_O {
 // class APFloat_O : public core::ExternalObject_O {
 // class APInt_O : public core::ExternalObject_O {
 // class IRBuilderBase_O : public core::ExternalObject_O {
 // class IRBuilder_O : public IRBuilderBase_O {
 // class Instruction_O : public User_O {
 // class StoreInst_O : public Instruction_O {
 // class FenceInst_O : public Instruction_O {
 // class AtomicCmpXchgInst_O : public Instruction_O {
 // class AtomicRMWInst_O : public Instruction_O {
 // class PHINode_O : public Instruction_O {
 // class CallInst_O : public Instruction_O {
 // class LandingPadInst_O : public Instruction_O {
 // class UnaryInstruction_O : public Instruction_O {
 // class AllocaInst_O : public UnaryInstruction_O {
 // class VAArgInst_O : public UnaryInstruction_O {
 // class LoadInst_O : public UnaryInstruction_O {
 // class TerminatorInst_O : public Instruction_O {
 // class BranchInst_O : public TerminatorInst_O {
 // class SwitchInst_O : public TerminatorInst_O {
 // class IndirectBrInst_O : public TerminatorInst_O {
 // class InvokeInst_O : public TerminatorInst_O {
 // class ResumeInst_O : public TerminatorInst_O {
 // class UnreachableInst_O : public TerminatorInst_O {
 // class ReturnInst_O : public TerminatorInst_O {
 // class ConstantFP_O : public Constant_O {
 // class ConstantInt_O : public Constant_O {
 // class ConstantStruct_O : public Constant_O {
 // class UndefValue_O : public Constant_O {
 // class ConstantPointerNull_O : public Constant_O {
 // class MDNode_O : public Metadata_O {
 // class MDString_O : public Metadata_O {
 // class ValueAsMetadata_O : public Metadata_O {
 // class NamedMDNode_O : public core::ExternalObject_O {
 // class Function_O : public GlobalValue_O {
 // class BasicBlock_O : public Value_O {
 // class Argument_O : public Value_O {
 // class Type_O : public core::ExternalObject_O {
 // class FunctionType_O : public Type_O {
 // class IntegerType_O : public Type_O {
 // class CompositeType_O : public Type_O {
 // class StructType_O : public CompositeType_O {
 // class SequentialType_O : public CompositeType_O {
 // class PointerType_O : public SequentialType_O {
 // class ArrayType_O : public SequentialType_O {
 // class VectorType_O : public SequentialType_O {
 // class IterateCons_O : public core::Iterator_O
 // class ObjectSet_O : public T_O {
 // class QDomNode_O : public core::CxxObject_O
 // class StringList_O : public T_O {
 // class StringSet_O : public core::CxxObject_O {
 // class SymbolList_O : public core::CxxObject_O {
 // class SymbolSet_O : public T_O {
 // class BoundingBox_O : public core::CxxObject_O
 // class Color_O : public core::CxxObject_O
 // class CoordinateArray_O : public core::Vector_O
 // class CoordinateArrayWithHash_O : public CoordinateArray_O
 // class OMatrix_O : public core::CxxObject_O
 // class OVector2_O : public core::CxxObject_O
 // class OVector3_O : public core::CxxObject_O
 // class Dimension_O : public core::CxxObject_O
 // class NamedUnit_O : public Unit_O
 // class Quantity_O : public core::CxxObject_O
 // class Unit_O : public core::CxxObject_O
#ifdef HEADER_INCLUDES
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/objectSet.h"
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/stringList.h"
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/symbolSet.h"
#include "/home/meister/Development/clasp/src/core/include/object.h"
#include "/home/meister/Development/clasp/src/asttooling/include/astVisitor.h"
#include "/home/meister/Development/clasp/src/cffi/include/cffi.h"
#include "/home/meister/Development/clasp/src/clbind/include/class_registry.h"
#include "/home/meister/Development/clasp/src/core/include/serialize.h"
#include "/home/meister/Development/clasp/src/core/include/array.h"
#include "/home/meister/Development/clasp/src/core/include/binder.h"
#include "/home/meister/Development/clasp/src/core/include/conditions.h"
#include "/home/meister/Development/clasp/src/core/include/character.h"
#include "/home/meister/Development/clasp/src/core/include/cons.h"
#include "/home/meister/Development/clasp/src/core/include/cxxObject.h"
#include "/home/meister/Development/clasp/src/core/include/fileSystem.h"
#include "/home/meister/Development/clasp/src/core/include/environment.h"
#include "/home/meister/Development/clasp/src/core/include/externalObject.h"
#include "/home/meister/Development/clasp/src/core/include/executables.h"
#include "/home/meister/Development/clasp/src/core/include/hashTable.h"
#include "/home/meister/Development/clasp/src/core/include/intArray.h"
#include "/home/meister/Development/clasp/src/core/include/primitives.h"
#include "/home/meister/Development/clasp/src/core/include/iterator.h"
#include "/home/meister/Development/clasp/src/core/include/lambdaListHandler.h"
#include "/home/meister/Development/clasp/src/core/include/userData.h"
#include "/home/meister/Development/clasp/src/core/include/loadTimeValues.h"
#include "/home/meister/Development/clasp/src/core/include/multiStringBuffer.h"
#include "/home/meister/Development/clasp/src/core/include/numbers.h"
#include "/home/meister/Development/clasp/src/core/include/package.h"
#include "/home/meister/Development/clasp/src/core/include/pathname.h"
#include "/home/meister/Development/clasp/src/core/include/pointer.h"
#include "/home/meister/Development/clasp/src/core/include/posixTime.h"
#include "/home/meister/Development/clasp/src/core/include/random.h"
#include "/home/meister/Development/clasp/src/core/include/readtable.h"
#include "/home/meister/Development/clasp/src/core/include/reader.h"
#include "/home/meister/Development/clasp/src/core/include/record.h"
#include "/home/meister/Development/clasp/src/core/include/regex.h"
#include "/home/meister/Development/clasp/src/core/include/singleDispatchMethod.h"
#include "/home/meister/Development/clasp/src/core/include/smallMap.h"
#include "/home/meister/Development/clasp/src/core/include/smallMultimap.h"
#include "/home/meister/Development/clasp/src/core/include/sourceFileInfo.h"
#include "/home/meister/Development/clasp/src/core/include/standardObject.h"
#include "/home/meister/Development/clasp/src/core/include/lispStream.h"
#include "/home/meister/Development/clasp/src/core/include/structureObject.h"
#include "/home/meister/Development/clasp/src/core/include/symbolToEnumConverter.h"
#include "/home/meister/Development/clasp/src/core/include/symbol.h"
#include "/home/meister/Development/clasp/src/core/include/lispList.h"
#include "/home/meister/Development/clasp/src/core/include/weakHashTable.h"
#include "/home/meister/Development/clasp/src/core/include/weakKeyMapping.h"
#include "/home/meister/Development/clasp/src/core/include/weakPointer.h"
#include "/home/meister/Development/clasp/src/core/include/wrappedPointer.h"
#include "/home/meister/Development/clasp/src/llvmo/include/llvmoExpose.h"
#include "/home/meister/Development/clasp/src/llvmo/include/debugInfoExpose.h"
#include "/home/meister/Development/clasp/src/llvmo/include/debugLoc.h"
#include "/home/meister/Development/clasp/src/llvmo/include/insertPoint.h"
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/iterateCons.h"
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/quickDom.h"
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/stringSet.h"
#include "/home/meister/Development/clasp/projects/cando/src/adapt/include/symbolList.h"
#include "/home/meister/Development/clasp/src/core/include/activationFrame.h"
#include "/home/meister/Development/clasp/src/core/include/arrayDisplaced.h"
#include "/home/meister/Development/clasp/src/core/include/arrayObjects.h"
#include "/home/meister/Development/clasp/src/core/include/hashTableEq.h"
#include "/home/meister/Development/clasp/src/core/include/hashTableEql.h"
#include "/home/meister/Development/clasp/src/core/include/hashTableEqual.h"
#include "/home/meister/Development/clasp/src/core/include/hashTableEqualp.h"
#include "/home/meister/Development/clasp/src/core/include/instance.h"
#include "/home/meister/Development/clasp/src/core/include/metaobject.h"
#include "/home/meister/Development/clasp/src/core/include/null.h"
#include "/home/meister/Development/clasp/src/core/include/singleDispatchEffectiveMethodFunction.h"
#include "/home/meister/Development/clasp/src/core/include/singleDispatchGenericFunction.h"
#include "/home/meister/Development/clasp/src/core/include/specialForm.h"
#include "/home/meister/Development/clasp/src/core/include/lispVector.h"
#include "/home/meister/Development/clasp/projects/cando/src/geom/include/boundingBox.h"
#include "/home/meister/Development/clasp/projects/cando/src/geom/include/color.h"
#include "/home/meister/Development/clasp/projects/cando/src/geom/include/omatrix.h"
#include "/home/meister/Development/clasp/projects/cando/src/geom/include/ovector2.h"
#include "/home/meister/Development/clasp/projects/cando/src/geom/include/ovector3.h"
#include "/home/meister/Development/clasp/projects/cando/src/units/include/dimension.h"
#include "/home/meister/Development/clasp/projects/cando/src/units/include/quantity.h"
#include "/home/meister/Development/clasp/projects/cando/src/units/include/unit.h"
#include "/home/meister/Development/clasp/src/core/include/bitVector.h"
#include "/home/meister/Development/clasp/src/core/include/sexpLoadArchive.h"
#include "/home/meister/Development/clasp/src/core/include/sexpSaveArchive.h"
#include "/home/meister/Development/clasp/src/core/include/specializer.h"
#include "/home/meister/Development/clasp/src/core/include/lispString.h"
#include "/home/meister/Development/clasp/src/core/include/newVectorObjects.h"
#include "/home/meister/Development/clasp/projects/cando/src/geom/include/coordinateArray.h"
#include "/home/meister/Development/clasp/projects/cando/src/units/include/namedUnit.h"
#include "/home/meister/Development/clasp/src/core/include/metaClass.h"
#include "/home/meister/Development/clasp/src/core/include/str.h"
#include "/home/meister/Development/clasp/src/core/include/newVectorObjectsWithFillPtr.h"
#include "/home/meister/Development/clasp/src/core/include/bignum.h"
#include "/home/meister/Development/clasp/src/core/include/builtInClass.h"
#include "/home/meister/Development/clasp/src/core/include/cxxClass.h"
#include "/home/meister/Development/clasp/src/core/include/forwardReferencedClass.h"
#include "/home/meister/Development/clasp/src/core/include/stdClass.h"
#include "/home/meister/Development/clasp/src/core/include/strWithFillPtr.h"
#include "/home/meister/Development/clasp/src/core/include/structureClass.h"
#include "/home/meister/Development/clasp/src/clbind/include/class_rep.h"
#include "/home/meister/Development/clasp/src/core/include/funcallableStandardClass.h"
#include "/home/meister/Development/clasp/src/core/include/standardClass.h"
#endif // HEADER_INCLUDES
#undef HEADER_INCLUDES
#if defined(SET_SYMBOL) || defined(ALL_STAGES)
// requires LOOKUP_SYMBOL(pkg,symbolName) be defined
#ifdef Use_AdaptPkg
adapt::ObjectSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::ObjectSet_O::static_packageName(),adapt::ObjectSet_O::static_className()));
#endif
#ifdef Use_AdaptPkg
adapt::StringList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::StringList_O::static_packageName(),adapt::StringList_O::static_className()));
#endif
#ifdef Use_AdaptPkg
adapt::SymbolSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::SymbolSet_O::static_packageName(),adapt::SymbolSet_O::static_className()));
#endif
#ifdef Use_CorePkg
core::T_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::T_O::static_packageName(),core::T_O::static_className()));
#endif
#ifdef Use_AstToolingPkg
asttooling::AstVisitor_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(asttooling::AstVisitor_O::static_packageName(),asttooling::AstVisitor_O::static_className()));
#endif
#ifdef Use_CffiPkg
cffi::Pointer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(cffi::Pointer_O::static_packageName(),cffi::Pointer_O::static_className()));
#endif
#ifdef Use_ClbindPkg
clbind::ClassRegistry_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(clbind::ClassRegistry_O::static_packageName(),clbind::ClassRegistry_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Archive_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Archive_O::static_packageName(),core::Archive_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Array_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Array_O::static_packageName(),core::Array_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Binder_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Binder_O::static_packageName(),core::Binder_O::static_className()));
#endif
#ifdef Use_CorePkg
core::CandoException_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::CandoException_O::static_packageName(),core::CandoException_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Character_dummy_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Character_dummy_O::static_packageName(),core::Character_dummy_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Cons_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Cons_O::static_packageName(),core::Cons_O::static_className()));
#endif
#ifdef Use_CorePkg
core::CxxObject_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::CxxObject_O::static_packageName(),core::CxxObject_O::static_className()));
#endif
#ifdef Use_CorePkg
core::DirectoryEntry_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::DirectoryEntry_O::static_packageName(),core::DirectoryEntry_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Environment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Environment_O::static_packageName(),core::Environment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ExternalObject_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ExternalObject_O::static_packageName(),core::ExternalObject_O::static_className()));
#endif
#ifdef Use_CorePkg
core::FileStatus_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::FileStatus_O::static_packageName(),core::FileStatus_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Function_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Function_O::static_packageName(),core::Function_O::static_className()));
#endif
#ifdef Use_CorePkg
core::HashTable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::HashTable_O::static_packageName(),core::HashTable_O::static_className()));
#endif
#ifdef Use_CorePkg
core::IntArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::IntArray_O::static_packageName(),core::IntArray_O::static_className()));
#endif
#ifdef Use_CorePkg
core::InvocationHistoryFrameIterator_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::InvocationHistoryFrameIterator_O::static_packageName(),core::InvocationHistoryFrameIterator_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Iterator_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Iterator_O::static_packageName(),core::Iterator_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LambdaListHandler_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LambdaListHandler_O::static_packageName(),core::LambdaListHandler_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LightUserData_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LightUserData_O::static_packageName(),core::LightUserData_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LoadTimeValues_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LoadTimeValues_O::static_packageName(),core::LoadTimeValues_O::static_className()));
#endif
#ifdef Use_CorePkg
core::MultiStringBuffer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::MultiStringBuffer_O::static_packageName(),core::MultiStringBuffer_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Number_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Number_O::static_packageName(),core::Number_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Package_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Package_O::static_packageName(),core::Package_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Path_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Path_O::static_packageName(),core::Path_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Pathname_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Pathname_O::static_packageName(),core::Pathname_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Pointer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Pointer_O::static_packageName(),core::Pointer_O::static_className()));
#endif
#ifdef Use_CorePkg
core::PosixTimeDuration_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::PosixTimeDuration_O::static_packageName(),core::PosixTimeDuration_O::static_className()));
#endif
#ifdef Use_CorePkg
core::PosixTime_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::PosixTime_O::static_packageName(),core::PosixTime_O::static_className()));
#endif
#ifdef Use_CorePkg
core::RandomState_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::RandomState_O::static_packageName(),core::RandomState_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ReadTable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ReadTable_O::static_packageName(),core::ReadTable_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Reader_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Reader_O::static_packageName(),core::Reader_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Record_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Record_O::static_packageName(),core::Record_O::static_className()));
#endif
#ifdef Use_CorePkg
core::RegexMatch_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::RegexMatch_O::static_packageName(),core::RegexMatch_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Regex_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Regex_O::static_packageName(),core::Regex_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SNode_O::static_packageName(),core::SNode_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SingleDispatchMethod_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SingleDispatchMethod_O::static_packageName(),core::SingleDispatchMethod_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SmallMap_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SmallMap_O::static_packageName(),core::SmallMap_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SmallMultimap_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SmallMultimap_O::static_packageName(),core::SmallMultimap_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SourceFileInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SourceFileInfo_O::static_packageName(),core::SourceFileInfo_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SourceManager_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SourceManager_O::static_packageName(),core::SourceManager_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SourcePosInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SourcePosInfo_O::static_packageName(),core::SourcePosInfo_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StandardObject_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StandardObject_O::static_packageName(),core::StandardObject_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Stream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Stream_O::static_packageName(),core::Stream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StructureObject_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StructureObject_O::static_packageName(),core::StructureObject_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SymbolToEnumConverter_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SymbolToEnumConverter_O::static_packageName(),core::SymbolToEnumConverter_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Symbol_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Symbol_O::static_packageName(),core::Symbol_O::static_className()));
#endif
#ifdef Use_CorePkg
core::VaList_dummy_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::VaList_dummy_O::static_packageName(),core::VaList_dummy_O::static_className()));
#endif
#ifdef Use_CorePkg
core::WeakHashTable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::WeakHashTable_O::static_packageName(),core::WeakHashTable_O::static_className()));
#endif
#ifdef Use_CorePkg
core::WeakKeyMapping_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::WeakKeyMapping_O::static_packageName(),core::WeakKeyMapping_O::static_className()));
#endif
#ifdef Use_CorePkg
core::WeakPointer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::WeakPointer_O::static_packageName(),core::WeakPointer_O::static_className()));
#endif
#ifdef Use_CorePkg
core::WrappedPointer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::WrappedPointer_O::static_packageName(),core::WrappedPointer_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::AttributeSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::AttributeSet_O::static_packageName(),llvmo::AttributeSet_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Attribute_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Attribute_O::static_packageName(),llvmo::Attribute_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DebugInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DebugInfo_O::static_packageName(),llvmo::DebugInfo_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DebugLoc_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DebugLoc_O::static_packageName(),llvmo::DebugLoc_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::InsertPoint_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::InsertPoint_O::static_packageName(),llvmo::InsertPoint_O::static_className()));
#endif
#ifdef Use_AdaptPkg
adapt::IterateCons_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::IterateCons_O::static_packageName(),adapt::IterateCons_O::static_className()));
#endif
#ifdef Use_AdaptPkg
adapt::QDomNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::QDomNode_O::static_packageName(),adapt::QDomNode_O::static_className()));
#endif
#ifdef Use_AdaptPkg
adapt::StringSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::StringSet_O::static_packageName(),adapt::StringSet_O::static_className()));
#endif
#ifdef Use_AdaptPkg
adapt::SymbolList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::SymbolList_O::static_packageName(),adapt::SymbolList_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ActivationFrame_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ActivationFrame_O::static_packageName(),core::ActivationFrame_O::static_className()));
#endif
#ifdef Use_CorePkg
core::AnsiStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::AnsiStream_O::static_packageName(),core::AnsiStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ArrayDisplaced_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ArrayDisplaced_O::static_packageName(),core::ArrayDisplaced_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ArrayObjects_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ArrayObjects_O::static_packageName(),core::ArrayObjects_O::static_className()));
#endif
#ifdef Use_CorePkg
core::BranchSNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::BranchSNode_O::static_packageName(),core::BranchSNode_O::static_className()));
#endif
#ifdef Use_CorePkg
core::CompiledFunction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::CompiledFunction_O::static_packageName(),core::CompiledFunction_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Complex_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Complex_O::static_packageName(),core::Complex_O::static_className()));
#endif
#ifdef Use_CorePkg
core::DirectoryIterator_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::DirectoryIterator_O::static_packageName(),core::DirectoryIterator_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ForeignData_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ForeignData_O::static_packageName(),core::ForeignData_O::static_className()));
#endif
#ifdef Use_CorePkg
core::GlueEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::GlueEnvironment_O::static_packageName(),core::GlueEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::HashTableEq_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::HashTableEq_O::static_packageName(),core::HashTableEq_O::static_className()));
#endif
#ifdef Use_CorePkg
core::HashTableEql_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::HashTableEql_O::static_packageName(),core::HashTableEql_O::static_className()));
#endif
#ifdef Use_CorePkg
core::HashTableEqual_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::HashTableEqual_O::static_packageName(),core::HashTableEqual_O::static_className()));
#endif
#ifdef Use_CorePkg
core::HashTableEqualp_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::HashTableEqualp_O::static_packageName(),core::HashTableEqualp_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Instance_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Instance_O::static_packageName(),core::Instance_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LeafSNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LeafSNode_O::static_packageName(),core::LeafSNode_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LexicalEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LexicalEnvironment_O::static_packageName(),core::LexicalEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LoadArchive_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LoadArchive_O::static_packageName(),core::LoadArchive_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LogicalPathname_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LogicalPathname_O::static_packageName(),core::LogicalPathname_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Metaobject_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Metaobject_O::static_packageName(),core::Metaobject_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Null_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Null_O::static_packageName(),core::Null_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Real_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Real_O::static_packageName(),core::Real_O::static_className()));
#endif
#ifdef Use_CorePkg
core::RecursiveDirectoryIterator_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::RecursiveDirectoryIterator_O::static_packageName(),core::RecursiveDirectoryIterator_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SaveArchive_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SaveArchive_O::static_packageName(),core::SaveArchive_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SingleDispatchEffectiveMethodFunction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SingleDispatchEffectiveMethodFunction_O::static_packageName(),core::SingleDispatchEffectiveMethodFunction_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SingleDispatchGenericFunction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SingleDispatchGenericFunction_O::static_packageName(),core::SingleDispatchGenericFunction_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SpecialForm_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SpecialForm_O::static_packageName(),core::SpecialForm_O::static_className()));
#endif
#ifdef Use_CorePkg
core::UserData_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::UserData_O::static_packageName(),core::UserData_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Vector_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Vector_O::static_packageName(),core::Vector_O::static_className()));
#endif
#ifdef Use_CorePkg
core::WeakKeyHashTable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::WeakKeyHashTable_O::static_packageName(),core::WeakKeyHashTable_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::BoundingBox_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::BoundingBox_O::static_packageName(),geom::BoundingBox_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::Color_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::Color_O::static_packageName(),geom::Color_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::OMatrix_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::OMatrix_O::static_packageName(),geom::OMatrix_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::OVector2_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::OVector2_O::static_packageName(),geom::OVector2_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::OVector3_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::OVector3_O::static_packageName(),geom::OVector3_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::APFloat_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::APFloat_O::static_packageName(),llvmo::APFloat_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::APInt_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::APInt_O::static_packageName(),llvmo::APInt_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIArray_O::static_packageName(),llvmo::DIArray_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIBasicType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIBasicType_O::static_packageName(),llvmo::DIBasicType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIBuilder_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIBuilder_O::static_packageName(),llvmo::DIBuilder_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DICompileUnit_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DICompileUnit_O::static_packageName(),llvmo::DICompileUnit_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DICompositeType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DICompositeType_O::static_packageName(),llvmo::DICompositeType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIDerivedType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIDerivedType_O::static_packageName(),llvmo::DIDerivedType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIDescriptor_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIDescriptor_O::static_packageName(),llvmo::DIDescriptor_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIFile_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIFile_O::static_packageName(),llvmo::DIFile_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DILexicalBlock_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DILexicalBlock_O::static_packageName(),llvmo::DILexicalBlock_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIScope_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIScope_O::static_packageName(),llvmo::DIScope_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DISubprogram_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DISubprogram_O::static_packageName(),llvmo::DISubprogram_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DISubroutineType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DISubroutineType_O::static_packageName(),llvmo::DISubroutineType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DITypeArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DITypeArray_O::static_packageName(),llvmo::DITypeArray_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DIType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DIType_O::static_packageName(),llvmo::DIType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DataLayout_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DataLayout_O::static_packageName(),llvmo::DataLayout_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::EngineBuilder_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::EngineBuilder_O::static_packageName(),llvmo::EngineBuilder_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ExecutionEngine_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ExecutionEngine_O::static_packageName(),llvmo::ExecutionEngine_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::IRBuilderBase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::IRBuilderBase_O::static_packageName(),llvmo::IRBuilderBase_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::LLVMContext_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::LLVMContext_O::static_packageName(),llvmo::LLVMContext_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Linker_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Linker_O::static_packageName(),llvmo::Linker_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::MCSubtargetInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::MCSubtargetInfo_O::static_packageName(),llvmo::MCSubtargetInfo_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Metadata_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Metadata_O::static_packageName(),llvmo::Metadata_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Module_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Module_O::static_packageName(),llvmo::Module_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::NamedMDNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::NamedMDNode_O::static_packageName(),llvmo::NamedMDNode_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::PassManagerBase_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::PassManagerBase_O::static_packageName(),llvmo::PassManagerBase_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::PassManagerBuilder_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::PassManagerBuilder_O::static_packageName(),llvmo::PassManagerBuilder_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Pass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Pass_O::static_packageName(),llvmo::Pass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::TargetMachine_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::TargetMachine_O::static_packageName(),llvmo::TargetMachine_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::TargetOptions_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::TargetOptions_O::static_packageName(),llvmo::TargetOptions_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Target_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Target_O::static_packageName(),llvmo::Target_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Triple_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Triple_O::static_packageName(),llvmo::Triple_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Type_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Type_O::static_packageName(),llvmo::Type_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Value_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Value_O::static_packageName(),llvmo::Value_O::static_className()));
#endif
#ifdef Use_UnitsPkg
units::Dimension_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::Dimension_O::static_packageName(),units::Dimension_O::static_className()));
#endif
#ifdef Use_UnitsPkg
units::Quantity_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::Quantity_O::static_packageName(),units::Quantity_O::static_className()));
#endif
#ifdef Use_UnitsPkg
units::Unit_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::Unit_O::static_packageName(),units::Unit_O::static_className()));
#endif
#ifdef Use_CorePkg
core::BitVector_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::BitVector_O::static_packageName(),core::BitVector_O::static_className()));
#endif
#ifdef Use_CorePkg
core::BroadcastStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::BroadcastStream_O::static_packageName(),core::BroadcastStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::CompileTimeEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::CompileTimeEnvironment_O::static_packageName(),core::CompileTimeEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ConcatenatedStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ConcatenatedStream_O::static_packageName(),core::ConcatenatedStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::EchoStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::EchoStream_O::static_packageName(),core::EchoStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::FileStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::FileStream_O::static_packageName(),core::FileStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Float_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Float_O::static_packageName(),core::Float_O::static_className()));
#endif
#ifdef Use_CorePkg
core::FunctionFrame_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::FunctionFrame_O::static_packageName(),core::FunctionFrame_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Rational_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Rational_O::static_packageName(),core::Rational_O::static_className()));
#endif
#ifdef Use_CorePkg
core::RuntimeVisibleEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::RuntimeVisibleEnvironment_O::static_packageName(),core::RuntimeVisibleEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SexpLoadArchive_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SexpLoadArchive_O::static_packageName(),core::SexpLoadArchive_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SexpSaveArchive_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SexpSaveArchive_O::static_packageName(),core::SexpSaveArchive_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Specializer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Specializer_O::static_packageName(),core::Specializer_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StringStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StringStream_O::static_packageName(),core::StringStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::String_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::String_O::static_packageName(),core::String_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SynonymStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SynonymStream_O::static_packageName(),core::SynonymStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::TagbodyFrame_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::TagbodyFrame_O::static_packageName(),core::TagbodyFrame_O::static_className()));
#endif
#ifdef Use_CorePkg
core::TwoWayStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::TwoWayStream_O::static_packageName(),core::TwoWayStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ValueFrame_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ValueFrame_O::static_packageName(),core::ValueFrame_O::static_className()));
#endif
#ifdef Use_CorePkg
core::VectorObjects_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::VectorObjects_O::static_packageName(),core::VectorObjects_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::CoordinateArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::CoordinateArray_O::static_packageName(),geom::CoordinateArray_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Argument_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Argument_O::static_packageName(),llvmo::Argument_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::BasicBlock_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::BasicBlock_O::static_packageName(),llvmo::BasicBlock_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::CompositeType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::CompositeType_O::static_packageName(),llvmo::CompositeType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::FunctionPassManager_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::FunctionPassManager_O::static_packageName(),llvmo::FunctionPassManager_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::FunctionPass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::FunctionPass_O::static_packageName(),llvmo::FunctionPass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::FunctionType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::FunctionType_O::static_packageName(),llvmo::FunctionType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::IRBuilder_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::IRBuilder_O::static_packageName(),llvmo::IRBuilder_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::IntegerType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::IntegerType_O::static_packageName(),llvmo::IntegerType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::LLVMTargetMachine_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::LLVMTargetMachine_O::static_packageName(),llvmo::LLVMTargetMachine_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::MDNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::MDNode_O::static_packageName(),llvmo::MDNode_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::MDString_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::MDString_O::static_packageName(),llvmo::MDString_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ModulePass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ModulePass_O::static_packageName(),llvmo::ModulePass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::PassManager_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::PassManager_O::static_packageName(),llvmo::PassManager_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::TargetSubtargetInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::TargetSubtargetInfo_O::static_packageName(),llvmo::TargetSubtargetInfo_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::User_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::User_O::static_packageName(),llvmo::User_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ValueAsMetadata_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ValueAsMetadata_O::static_packageName(),llvmo::ValueAsMetadata_O::static_className()));
#endif
#ifdef Use_UnitsPkg
units::NamedUnit_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::NamedUnit_O::static_packageName(),units::NamedUnit_O::static_className()));
#endif
#ifdef Use_CorePkg
core::BlockEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::BlockEnvironment_O::static_packageName(),core::BlockEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::CatchEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::CatchEnvironment_O::static_packageName(),core::CatchEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Class_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Class_O::static_packageName(),core::Class_O::static_className()));
#endif
#ifdef Use_CorePkg
core::DoubleFloat_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::DoubleFloat_O::static_packageName(),core::DoubleFloat_O::static_className()));
#endif
#ifdef Use_CorePkg
core::FunctionContainerEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::FunctionContainerEnvironment_O::static_packageName(),core::FunctionContainerEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::FunctionValueEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::FunctionValueEnvironment_O::static_packageName(),core::FunctionValueEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::IOFileStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::IOFileStream_O::static_packageName(),core::IOFileStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::IOStreamStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::IOStreamStream_O::static_packageName(),core::IOStreamStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Integer_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Integer_O::static_packageName(),core::Integer_O::static_className()));
#endif
#ifdef Use_CorePkg
core::LongFloat_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::LongFloat_O::static_packageName(),core::LongFloat_O::static_className()));
#endif
#ifdef Use_CorePkg
core::MacroletEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::MacroletEnvironment_O::static_packageName(),core::MacroletEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Ratio_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Ratio_O::static_packageName(),core::Ratio_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ShortFloat_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ShortFloat_O::static_packageName(),core::ShortFloat_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SimpleBitVector_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SimpleBitVector_O::static_packageName(),core::SimpleBitVector_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SingleFloat_dummy_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SingleFloat_dummy_O::static_packageName(),core::SingleFloat_dummy_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StackValueEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StackValueEnvironment_O::static_packageName(),core::StackValueEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Str_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Str_O::static_packageName(),core::Str_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StringInputStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StringInputStream_O::static_packageName(),core::StringInputStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StringOutputStream_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StringOutputStream_O::static_packageName(),core::StringOutputStream_O::static_className()));
#endif
#ifdef Use_CorePkg
core::SymbolMacroletEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::SymbolMacroletEnvironment_O::static_packageName(),core::SymbolMacroletEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::TagbodyEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::TagbodyEnvironment_O::static_packageName(),core::TagbodyEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::UnwindProtectEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::UnwindProtectEnvironment_O::static_packageName(),core::UnwindProtectEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ValueEnvironment_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ValueEnvironment_O::static_packageName(),core::ValueEnvironment_O::static_className()));
#endif
#ifdef Use_CorePkg
core::VectorObjectsWithFillPtr_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::VectorObjectsWithFillPtr_O::static_packageName(),core::VectorObjectsWithFillPtr_O::static_className()));
#endif
#ifdef Use_GeomPkg
geom::CoordinateArrayWithHash_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::CoordinateArrayWithHash_O::static_packageName(),geom::CoordinateArrayWithHash_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Constant_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Constant_O::static_packageName(),llvmo::Constant_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ImmutablePass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ImmutablePass_O::static_packageName(),llvmo::ImmutablePass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Instruction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Instruction_O::static_packageName(),llvmo::Instruction_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::SequentialType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::SequentialType_O::static_packageName(),llvmo::SequentialType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::StructType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::StructType_O::static_packageName(),llvmo::StructType_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Bignum_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Bignum_O::static_packageName(),core::Bignum_O::static_className()));
#endif
#ifdef Use_CorePkg
core::BuiltInClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::BuiltInClass_O::static_packageName(),core::BuiltInClass_O::static_className()));
#endif
#ifdef Use_CorePkg
core::CxxClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::CxxClass_O::static_packageName(),core::CxxClass_O::static_className()));
#endif
#ifdef Use_CorePkg
core::Fixnum_dummy_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::Fixnum_dummy_O::static_packageName(),core::Fixnum_dummy_O::static_className()));
#endif
#ifdef Use_CorePkg
core::ForwardReferencedClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::ForwardReferencedClass_O::static_packageName(),core::ForwardReferencedClass_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StdClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StdClass_O::static_packageName(),core::StdClass_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StrWithFillPtr_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StrWithFillPtr_O::static_packageName(),core::StrWithFillPtr_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StructureClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StructureClass_O::static_packageName(),core::StructureClass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ArrayType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ArrayType_O::static_packageName(),llvmo::ArrayType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::AtomicCmpXchgInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::AtomicCmpXchgInst_O::static_packageName(),llvmo::AtomicCmpXchgInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::AtomicRMWInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::AtomicRMWInst_O::static_packageName(),llvmo::AtomicRMWInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::BlockAddress_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::BlockAddress_O::static_packageName(),llvmo::BlockAddress_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::CallInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::CallInst_O::static_packageName(),llvmo::CallInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantArray_O::static_packageName(),llvmo::ConstantArray_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantDataSequential_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantDataSequential_O::static_packageName(),llvmo::ConstantDataSequential_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantExpr_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantExpr_O::static_packageName(),llvmo::ConstantExpr_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantFP_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantFP_O::static_packageName(),llvmo::ConstantFP_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantInt_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantInt_O::static_packageName(),llvmo::ConstantInt_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantPointerNull_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantPointerNull_O::static_packageName(),llvmo::ConstantPointerNull_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantStruct_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantStruct_O::static_packageName(),llvmo::ConstantStruct_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::DataLayoutPass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::DataLayoutPass_O::static_packageName(),llvmo::DataLayoutPass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::FenceInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::FenceInst_O::static_packageName(),llvmo::FenceInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::GlobalValue_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::GlobalValue_O::static_packageName(),llvmo::GlobalValue_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::LandingPadInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::LandingPadInst_O::static_packageName(),llvmo::LandingPadInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::PHINode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::PHINode_O::static_packageName(),llvmo::PHINode_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::PointerType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::PointerType_O::static_packageName(),llvmo::PointerType_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::StoreInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::StoreInst_O::static_packageName(),llvmo::StoreInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::TargetLibraryInfo_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::TargetLibraryInfo_O::static_packageName(),llvmo::TargetLibraryInfo_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::TerminatorInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::TerminatorInst_O::static_packageName(),llvmo::TerminatorInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::UnaryInstruction_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::UnaryInstruction_O::static_packageName(),llvmo::UnaryInstruction_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::UndefValue_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::UndefValue_O::static_packageName(),llvmo::UndefValue_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::VectorType_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::VectorType_O::static_packageName(),llvmo::VectorType_O::static_className()));
#endif
#ifdef Use_ClbindPkg
clbind::ClassRep_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(clbind::ClassRep_O::static_packageName(),clbind::ClassRep_O::static_className()));
#endif
#ifdef Use_CorePkg
core::FuncallableStandardClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::FuncallableStandardClass_O::static_packageName(),core::FuncallableStandardClass_O::static_className()));
#endif
#ifdef Use_CorePkg
core::StandardClass_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(core::StandardClass_O::static_packageName(),core::StandardClass_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::AllocaInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::AllocaInst_O::static_packageName(),llvmo::AllocaInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::BranchInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::BranchInst_O::static_packageName(),llvmo::BranchInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ConstantDataArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ConstantDataArray_O::static_packageName(),llvmo::ConstantDataArray_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::Function_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::Function_O::static_packageName(),llvmo::Function_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::GlobalVariable_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::GlobalVariable_O::static_packageName(),llvmo::GlobalVariable_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::IndirectBrInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::IndirectBrInst_O::static_packageName(),llvmo::IndirectBrInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::InvokeInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::InvokeInst_O::static_packageName(),llvmo::InvokeInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::LoadInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::LoadInst_O::static_packageName(),llvmo::LoadInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ResumeInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ResumeInst_O::static_packageName(),llvmo::ResumeInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::ReturnInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::ReturnInst_O::static_packageName(),llvmo::ReturnInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::SwitchInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::SwitchInst_O::static_packageName(),llvmo::SwitchInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::UnreachableInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::UnreachableInst_O::static_packageName(),llvmo::UnreachableInst_O::static_className()));
#endif
#ifdef Use_LlvmoPkg
llvmo::VAArgInst_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(llvmo::VAArgInst_O::static_packageName(),llvmo::VAArgInst_O::static_className()));
#endif
#endif // SET_SYMBOL
#undef SET_SYMBOL
#if defined(CREATE_CLASS) || defined(ALL_STAGES)
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__ObjectSet_Oval]"));
    core::BuiltInClass_sp classadapt__ObjectSet_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__ObjectSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__ObjectSet_Oval,_lisp,adapt::ObjectSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::ObjectSet_O>::id,adapt::ObjectSet_O::static_classSymbol());
    adapt::ObjectSet_O::___staticClass = classadapt__ObjectSet_Oval;
#ifdef USE_MPS
    adapt::ObjectSet_O::static_Kind = gctools::GCKind<adapt::ObjectSet_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__ObjectSet_Oval,adapt::ObjectSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::ObjectSet_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::ObjectSet_O>>::allocateClass();
        adapt::ObjectSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::ObjectSet_O::static_className() % (void*)(adapt::ObjectSet_O::static_allocator) );
    classadapt__ObjectSet_Oval->setCreator(adapt::ObjectSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::ObjectSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__ObjectSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__StringList_Oval]"));
    core::BuiltInClass_sp classadapt__StringList_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__StringList_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__StringList_Oval,_lisp,adapt::StringList_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::StringList_O>::id,adapt::StringList_O::static_classSymbol());
    adapt::StringList_O::___staticClass = classadapt__StringList_Oval;
#ifdef USE_MPS
    adapt::StringList_O::static_Kind = gctools::GCKind<adapt::StringList_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__StringList_Oval,adapt::StringList_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::StringList_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::StringList_O>>::allocateClass();
        adapt::StringList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::StringList_O::static_className() % (void*)(adapt::StringList_O::static_allocator) );
    classadapt__StringList_Oval->setCreator(adapt::StringList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::StringList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__StringList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__SymbolSet_Oval]"));
    core::BuiltInClass_sp classadapt__SymbolSet_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__SymbolSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__SymbolSet_Oval,_lisp,adapt::SymbolSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::SymbolSet_O>::id,adapt::SymbolSet_O::static_classSymbol());
    adapt::SymbolSet_O::___staticClass = classadapt__SymbolSet_Oval;
#ifdef USE_MPS
    adapt::SymbolSet_O::static_Kind = gctools::GCKind<adapt::SymbolSet_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__SymbolSet_Oval,adapt::SymbolSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::SymbolSet_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::SymbolSet_O>>::allocateClass();
        adapt::SymbolSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::SymbolSet_O::static_className() % (void*)(adapt::SymbolSet_O::static_allocator) );
    classadapt__SymbolSet_Oval->setCreator(adapt::SymbolSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::SymbolSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__SymbolSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__T_Oval]"));
    core::BuiltInClass_sp classcore__T_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__T_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__T_Oval,_lisp,core::T_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::T_O>::id,core::T_O::static_classSymbol());
    core::T_O::___staticClass = classcore__T_Oval;
#ifdef USE_MPS
    core::T_O::static_Kind = gctools::GCKind<core::T_O>::Kind;
#endif
    core::af_setf_findClass(classcore__T_Oval,core::T_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::T_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::T_O>>::allocateClass();
        core::T_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::T_O::static_className() % (void*)(core::T_O::static_allocator) );
    classcore__T_Oval->setCreator(core::T_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::T_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__T_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AstToolingPkg
    LOG(BF("Creating class[classasttooling__AstVisitor_Oval]"));
    core::BuiltInClass_sp classasttooling__AstVisitor_Oval = core::BuiltInClass_O::createUncollectable();
    classasttooling__AstVisitor_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classasttooling__AstVisitor_Oval,_lisp,asttooling::AstVisitor_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<asttooling::AstVisitor_O>::id,asttooling::AstVisitor_O::static_classSymbol());
    asttooling::AstVisitor_O::___staticClass = classasttooling__AstVisitor_Oval;
#ifdef USE_MPS
    asttooling::AstVisitor_O::static_Kind = gctools::GCKind<asttooling::AstVisitor_O>::Kind;
#endif
    core::af_setf_findClass(classasttooling__AstVisitor_Oval,asttooling::AstVisitor_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<asttooling::AstVisitor_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<asttooling::AstVisitor_O>>::allocateClass();
        asttooling::AstVisitor_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% asttooling::AstVisitor_O::static_className() % (void*)(asttooling::AstVisitor_O::static_allocator) );
    classasttooling__AstVisitor_Oval->setCreator(asttooling::AstVisitor_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % asttooling::AstVisitor_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classasttooling__AstVisitor_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CffiPkg
    LOG(BF("Creating class[classcffi__Pointer_Oval]"));
    core::BuiltInClass_sp classcffi__Pointer_Oval = core::BuiltInClass_O::createUncollectable();
    classcffi__Pointer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcffi__Pointer_Oval,_lisp,cffi::Pointer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<cffi::Pointer_O>::id,cffi::Pointer_O::static_classSymbol());
    cffi::Pointer_O::___staticClass = classcffi__Pointer_Oval;
#ifdef USE_MPS
    cffi::Pointer_O::static_Kind = gctools::GCKind<cffi::Pointer_O>::Kind;
#endif
    core::af_setf_findClass(classcffi__Pointer_Oval,cffi::Pointer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<cffi::Pointer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<cffi::Pointer_O>>::allocateClass();
        cffi::Pointer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% cffi::Pointer_O::static_className() % (void*)(cffi::Pointer_O::static_allocator) );
    classcffi__Pointer_Oval->setCreator(cffi::Pointer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % cffi::Pointer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcffi__Pointer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_ClbindPkg
    LOG(BF("Creating class[classclbind__ClassRegistry_Oval]"));
    core::BuiltInClass_sp classclbind__ClassRegistry_Oval = core::BuiltInClass_O::createUncollectable();
    classclbind__ClassRegistry_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classclbind__ClassRegistry_Oval,_lisp,clbind::ClassRegistry_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<clbind::ClassRegistry_O>::id,clbind::ClassRegistry_O::static_classSymbol());
    clbind::ClassRegistry_O::___staticClass = classclbind__ClassRegistry_Oval;
#ifdef USE_MPS
    clbind::ClassRegistry_O::static_Kind = gctools::GCKind<clbind::ClassRegistry_O>::Kind;
#endif
    core::af_setf_findClass(classclbind__ClassRegistry_Oval,clbind::ClassRegistry_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<clbind::ClassRegistry_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<clbind::ClassRegistry_O>>::allocateClass();
        clbind::ClassRegistry_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% clbind::ClassRegistry_O::static_className() % (void*)(clbind::ClassRegistry_O::static_allocator) );
    classclbind__ClassRegistry_Oval->setCreator(clbind::ClassRegistry_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % clbind::ClassRegistry_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classclbind__ClassRegistry_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Archive_Oval]"));
    core::BuiltInClass_sp classcore__Archive_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Archive_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Archive_Oval,_lisp,core::Archive_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Archive_O>::id,core::Archive_O::static_classSymbol());
    core::Archive_O::___staticClass = classcore__Archive_Oval;
#ifdef USE_MPS
    core::Archive_O::static_Kind = gctools::GCKind<core::Archive_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Archive_Oval,core::Archive_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Archive_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Archive_O>>::allocateClass();
        core::Archive_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Archive_O::static_className() % (void*)(core::Archive_O::static_allocator) );
    classcore__Archive_Oval->setCreator(core::Archive_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Archive_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Archive_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Array_Oval]"));
    core::BuiltInClass_sp classcore__Array_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Array_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Array_Oval,_lisp,core::Array_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Array_O>::id,core::Array_O::static_classSymbol());
    core::Array_O::___staticClass = classcore__Array_Oval;
#ifdef USE_MPS
    core::Array_O::static_Kind = gctools::GCKind<core::Array_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Array_Oval,core::Array_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Array_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Array_O>>::allocateClass();
        core::Array_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Array_O::static_className() % (void*)(core::Array_O::static_allocator) );
    classcore__Array_Oval->setCreator(core::Array_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Array_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Array_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Binder_Oval]"));
    core::BuiltInClass_sp classcore__Binder_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Binder_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Binder_Oval,_lisp,core::Binder_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Binder_O>::id,core::Binder_O::static_classSymbol());
    core::Binder_O::___staticClass = classcore__Binder_Oval;
#ifdef USE_MPS
    core::Binder_O::static_Kind = gctools::GCKind<core::Binder_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Binder_Oval,core::Binder_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Binder_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Binder_O>>::allocateClass();
        core::Binder_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Binder_O::static_className() % (void*)(core::Binder_O::static_allocator) );
    classcore__Binder_Oval->setCreator(core::Binder_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Binder_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Binder_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__CandoException_Oval]"));
    core::BuiltInClass_sp classcore__CandoException_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__CandoException_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__CandoException_Oval,_lisp,core::CandoException_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::CandoException_O>::id,core::CandoException_O::static_classSymbol());
    core::CandoException_O::___staticClass = classcore__CandoException_Oval;
#ifdef USE_MPS
    core::CandoException_O::static_Kind = gctools::GCKind<core::CandoException_O>::Kind;
#endif
    core::af_setf_findClass(classcore__CandoException_Oval,core::CandoException_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::CandoException_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::CandoException_O>>::allocateClass();
        core::CandoException_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::CandoException_O::static_className() % (void*)(core::CandoException_O::static_allocator) );
    classcore__CandoException_Oval->setCreator(core::CandoException_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::CandoException_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__CandoException_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Character_dummy_Oval]"));
    core::BuiltInClass_sp classcore__Character_dummy_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Character_dummy_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Character_dummy_Oval,_lisp,core::Character_dummy_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Character_dummy_O>::id,core::Character_dummy_O::static_classSymbol());
    core::Character_dummy_O::___staticClass = classcore__Character_dummy_Oval;
#ifdef USE_MPS
    core::Character_dummy_O::static_Kind = gctools::GCKind<core::Character_dummy_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Character_dummy_Oval,core::Character_dummy_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Character_dummy_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Character_dummy_O>>::allocateClass();
        core::Character_dummy_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Character_dummy_O::static_className() % (void*)(core::Character_dummy_O::static_allocator) );
    classcore__Character_dummy_Oval->setCreator(core::Character_dummy_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Character_dummy_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Character_dummy_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Cons_Oval]"));
    core::BuiltInClass_sp classcore__Cons_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Cons_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Cons_Oval,_lisp,core::Cons_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Cons_O>::id,core::Cons_O::static_classSymbol());
    core::Cons_O::___staticClass = classcore__Cons_Oval;
#ifdef USE_MPS
    core::Cons_O::static_Kind = gctools::GCKind<core::Cons_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Cons_Oval,core::Cons_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Cons_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Cons_O>>::allocateClass();
        core::Cons_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Cons_O::static_className() % (void*)(core::Cons_O::static_allocator) );
    classcore__Cons_Oval->setCreator(core::Cons_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Cons_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Cons_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__CxxObject_Oval]"));
    core::BuiltInClass_sp classcore__CxxObject_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__CxxObject_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__CxxObject_Oval,_lisp,core::CxxObject_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::CxxObject_O>::id,core::CxxObject_O::static_classSymbol());
    core::CxxObject_O::___staticClass = classcore__CxxObject_Oval;
#ifdef USE_MPS
    core::CxxObject_O::static_Kind = gctools::GCKind<core::CxxObject_O>::Kind;
#endif
    core::af_setf_findClass(classcore__CxxObject_Oval,core::CxxObject_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::CxxObject_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::CxxObject_O>>::allocateClass();
        core::CxxObject_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::CxxObject_O::static_className() % (void*)(core::CxxObject_O::static_allocator) );
    classcore__CxxObject_Oval->setCreator(core::CxxObject_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::CxxObject_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__CxxObject_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__DirectoryEntry_Oval]"));
    core::BuiltInClass_sp classcore__DirectoryEntry_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__DirectoryEntry_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__DirectoryEntry_Oval,_lisp,core::DirectoryEntry_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::DirectoryEntry_O>::id,core::DirectoryEntry_O::static_classSymbol());
    core::DirectoryEntry_O::___staticClass = classcore__DirectoryEntry_Oval;
#ifdef USE_MPS
    core::DirectoryEntry_O::static_Kind = gctools::GCKind<core::DirectoryEntry_O>::Kind;
#endif
    core::af_setf_findClass(classcore__DirectoryEntry_Oval,core::DirectoryEntry_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::DirectoryEntry_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::DirectoryEntry_O>>::allocateClass();
        core::DirectoryEntry_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::DirectoryEntry_O::static_className() % (void*)(core::DirectoryEntry_O::static_allocator) );
    classcore__DirectoryEntry_Oval->setCreator(core::DirectoryEntry_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::DirectoryEntry_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__DirectoryEntry_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Environment_Oval]"));
    core::BuiltInClass_sp classcore__Environment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Environment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Environment_Oval,_lisp,core::Environment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Environment_O>::id,core::Environment_O::static_classSymbol());
    core::Environment_O::___staticClass = classcore__Environment_Oval;
#ifdef USE_MPS
    core::Environment_O::static_Kind = gctools::GCKind<core::Environment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Environment_Oval,core::Environment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Environment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Environment_O>>::allocateClass();
        core::Environment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Environment_O::static_className() % (void*)(core::Environment_O::static_allocator) );
    classcore__Environment_Oval->setCreator(core::Environment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Environment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Environment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ExternalObject_Oval]"));
    core::BuiltInClass_sp classcore__ExternalObject_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ExternalObject_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ExternalObject_Oval,_lisp,core::ExternalObject_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ExternalObject_O>::id,core::ExternalObject_O::static_classSymbol());
    core::ExternalObject_O::___staticClass = classcore__ExternalObject_Oval;
#ifdef USE_MPS
    core::ExternalObject_O::static_Kind = gctools::GCKind<core::ExternalObject_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ExternalObject_Oval,core::ExternalObject_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ExternalObject_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ExternalObject_O>>::allocateClass();
        core::ExternalObject_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ExternalObject_O::static_className() % (void*)(core::ExternalObject_O::static_allocator) );
    classcore__ExternalObject_Oval->setCreator(core::ExternalObject_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ExternalObject_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ExternalObject_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__FileStatus_Oval]"));
    core::BuiltInClass_sp classcore__FileStatus_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__FileStatus_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__FileStatus_Oval,_lisp,core::FileStatus_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::FileStatus_O>::id,core::FileStatus_O::static_classSymbol());
    core::FileStatus_O::___staticClass = classcore__FileStatus_Oval;
#ifdef USE_MPS
    core::FileStatus_O::static_Kind = gctools::GCKind<core::FileStatus_O>::Kind;
#endif
    core::af_setf_findClass(classcore__FileStatus_Oval,core::FileStatus_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::FileStatus_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::FileStatus_O>>::allocateClass();
        core::FileStatus_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::FileStatus_O::static_className() % (void*)(core::FileStatus_O::static_allocator) );
    classcore__FileStatus_Oval->setCreator(core::FileStatus_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::FileStatus_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__FileStatus_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Function_Oval]"));
    core::BuiltInClass_sp classcore__Function_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Function_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Function_Oval,_lisp,core::Function_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Function_O>::id,core::Function_O::static_classSymbol());
    core::Function_O::___staticClass = classcore__Function_Oval;
#ifdef USE_MPS
    core::Function_O::static_Kind = gctools::GCKind<core::Function_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Function_Oval,core::Function_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Function_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Function_O>>::allocateClass();
        core::Function_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Function_O::static_className() % (void*)(core::Function_O::static_allocator) );
    classcore__Function_Oval->setCreator(core::Function_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Function_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Function_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__HashTable_Oval]"));
    core::BuiltInClass_sp classcore__HashTable_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__HashTable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__HashTable_Oval,_lisp,core::HashTable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::HashTable_O>::id,core::HashTable_O::static_classSymbol());
    core::HashTable_O::___staticClass = classcore__HashTable_Oval;
#ifdef USE_MPS
    core::HashTable_O::static_Kind = gctools::GCKind<core::HashTable_O>::Kind;
#endif
    core::af_setf_findClass(classcore__HashTable_Oval,core::HashTable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::HashTable_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::HashTable_O>>::allocateClass();
        core::HashTable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::HashTable_O::static_className() % (void*)(core::HashTable_O::static_allocator) );
    classcore__HashTable_Oval->setCreator(core::HashTable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::HashTable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__HashTable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__IntArray_Oval]"));
    core::BuiltInClass_sp classcore__IntArray_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__IntArray_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__IntArray_Oval,_lisp,core::IntArray_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::IntArray_O>::id,core::IntArray_O::static_classSymbol());
    core::IntArray_O::___staticClass = classcore__IntArray_Oval;
#ifdef USE_MPS
    core::IntArray_O::static_Kind = gctools::GCKind<core::IntArray_O>::Kind;
#endif
    core::af_setf_findClass(classcore__IntArray_Oval,core::IntArray_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::IntArray_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::IntArray_O>>::allocateClass();
        core::IntArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::IntArray_O::static_className() % (void*)(core::IntArray_O::static_allocator) );
    classcore__IntArray_Oval->setCreator(core::IntArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::IntArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__IntArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__InvocationHistoryFrameIterator_Oval]"));
    core::BuiltInClass_sp classcore__InvocationHistoryFrameIterator_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__InvocationHistoryFrameIterator_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__InvocationHistoryFrameIterator_Oval,_lisp,core::InvocationHistoryFrameIterator_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::InvocationHistoryFrameIterator_O>::id,core::InvocationHistoryFrameIterator_O::static_classSymbol());
    core::InvocationHistoryFrameIterator_O::___staticClass = classcore__InvocationHistoryFrameIterator_Oval;
#ifdef USE_MPS
    core::InvocationHistoryFrameIterator_O::static_Kind = gctools::GCKind<core::InvocationHistoryFrameIterator_O>::Kind;
#endif
    core::af_setf_findClass(classcore__InvocationHistoryFrameIterator_Oval,core::InvocationHistoryFrameIterator_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::InvocationHistoryFrameIterator_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::InvocationHistoryFrameIterator_O>>::allocateClass();
        core::InvocationHistoryFrameIterator_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::InvocationHistoryFrameIterator_O::static_className() % (void*)(core::InvocationHistoryFrameIterator_O::static_allocator) );
    classcore__InvocationHistoryFrameIterator_Oval->setCreator(core::InvocationHistoryFrameIterator_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::InvocationHistoryFrameIterator_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__InvocationHistoryFrameIterator_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Iterator_Oval]"));
    core::BuiltInClass_sp classcore__Iterator_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Iterator_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Iterator_Oval,_lisp,core::Iterator_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Iterator_O>::id,core::Iterator_O::static_classSymbol());
    core::Iterator_O::___staticClass = classcore__Iterator_Oval;
#ifdef USE_MPS
    core::Iterator_O::static_Kind = gctools::GCKind<core::Iterator_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Iterator_Oval,core::Iterator_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Iterator_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Iterator_O>>::allocateClass();
        core::Iterator_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Iterator_O::static_className() % (void*)(core::Iterator_O::static_allocator) );
    classcore__Iterator_Oval->setCreator(core::Iterator_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Iterator_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Iterator_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LambdaListHandler_Oval]"));
    core::BuiltInClass_sp classcore__LambdaListHandler_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LambdaListHandler_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LambdaListHandler_Oval,_lisp,core::LambdaListHandler_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LambdaListHandler_O>::id,core::LambdaListHandler_O::static_classSymbol());
    core::LambdaListHandler_O::___staticClass = classcore__LambdaListHandler_Oval;
#ifdef USE_MPS
    core::LambdaListHandler_O::static_Kind = gctools::GCKind<core::LambdaListHandler_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LambdaListHandler_Oval,core::LambdaListHandler_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LambdaListHandler_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LambdaListHandler_O>>::allocateClass();
        core::LambdaListHandler_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LambdaListHandler_O::static_className() % (void*)(core::LambdaListHandler_O::static_allocator) );
    classcore__LambdaListHandler_Oval->setCreator(core::LambdaListHandler_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LambdaListHandler_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LambdaListHandler_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LightUserData_Oval]"));
    core::BuiltInClass_sp classcore__LightUserData_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LightUserData_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LightUserData_Oval,_lisp,core::LightUserData_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LightUserData_O>::id,core::LightUserData_O::static_classSymbol());
    core::LightUserData_O::___staticClass = classcore__LightUserData_Oval;
#ifdef USE_MPS
    core::LightUserData_O::static_Kind = gctools::GCKind<core::LightUserData_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LightUserData_Oval,core::LightUserData_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LightUserData_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LightUserData_O>>::allocateClass();
        core::LightUserData_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LightUserData_O::static_className() % (void*)(core::LightUserData_O::static_allocator) );
    classcore__LightUserData_Oval->setCreator(core::LightUserData_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LightUserData_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LightUserData_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LoadTimeValues_Oval]"));
    core::BuiltInClass_sp classcore__LoadTimeValues_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LoadTimeValues_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LoadTimeValues_Oval,_lisp,core::LoadTimeValues_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LoadTimeValues_O>::id,core::LoadTimeValues_O::static_classSymbol());
    core::LoadTimeValues_O::___staticClass = classcore__LoadTimeValues_Oval;
#ifdef USE_MPS
    core::LoadTimeValues_O::static_Kind = gctools::GCKind<core::LoadTimeValues_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LoadTimeValues_Oval,core::LoadTimeValues_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LoadTimeValues_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LoadTimeValues_O>>::allocateClass();
        core::LoadTimeValues_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LoadTimeValues_O::static_className() % (void*)(core::LoadTimeValues_O::static_allocator) );
    classcore__LoadTimeValues_Oval->setCreator(core::LoadTimeValues_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LoadTimeValues_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LoadTimeValues_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__MultiStringBuffer_Oval]"));
    core::BuiltInClass_sp classcore__MultiStringBuffer_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__MultiStringBuffer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__MultiStringBuffer_Oval,_lisp,core::MultiStringBuffer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::MultiStringBuffer_O>::id,core::MultiStringBuffer_O::static_classSymbol());
    core::MultiStringBuffer_O::___staticClass = classcore__MultiStringBuffer_Oval;
#ifdef USE_MPS
    core::MultiStringBuffer_O::static_Kind = gctools::GCKind<core::MultiStringBuffer_O>::Kind;
#endif
    core::af_setf_findClass(classcore__MultiStringBuffer_Oval,core::MultiStringBuffer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::MultiStringBuffer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::MultiStringBuffer_O>>::allocateClass();
        core::MultiStringBuffer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::MultiStringBuffer_O::static_className() % (void*)(core::MultiStringBuffer_O::static_allocator) );
    classcore__MultiStringBuffer_Oval->setCreator(core::MultiStringBuffer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::MultiStringBuffer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__MultiStringBuffer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Number_Oval]"));
    core::BuiltInClass_sp classcore__Number_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Number_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Number_Oval,_lisp,core::Number_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Number_O>::id,core::Number_O::static_classSymbol());
    core::Number_O::___staticClass = classcore__Number_Oval;
#ifdef USE_MPS
    core::Number_O::static_Kind = gctools::GCKind<core::Number_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Number_Oval,core::Number_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Number_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Number_O>>::allocateClass();
        core::Number_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Number_O::static_className() % (void*)(core::Number_O::static_allocator) );
    classcore__Number_Oval->setCreator(core::Number_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Number_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Number_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Package_Oval]"));
    core::BuiltInClass_sp classcore__Package_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Package_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Package_Oval,_lisp,core::Package_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Package_O>::id,core::Package_O::static_classSymbol());
    core::Package_O::___staticClass = classcore__Package_Oval;
#ifdef USE_MPS
    core::Package_O::static_Kind = gctools::GCKind<core::Package_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Package_Oval,core::Package_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Package_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Package_O>>::allocateClass();
        core::Package_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Package_O::static_className() % (void*)(core::Package_O::static_allocator) );
    classcore__Package_Oval->setCreator(core::Package_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Package_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Package_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Path_Oval]"));
    core::BuiltInClass_sp classcore__Path_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Path_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Path_Oval,_lisp,core::Path_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Path_O>::id,core::Path_O::static_classSymbol());
    core::Path_O::___staticClass = classcore__Path_Oval;
#ifdef USE_MPS
    core::Path_O::static_Kind = gctools::GCKind<core::Path_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Path_Oval,core::Path_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Path_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Path_O>>::allocateClass();
        core::Path_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Path_O::static_className() % (void*)(core::Path_O::static_allocator) );
    classcore__Path_Oval->setCreator(core::Path_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Path_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Path_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Pathname_Oval]"));
    core::BuiltInClass_sp classcore__Pathname_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Pathname_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Pathname_Oval,_lisp,core::Pathname_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Pathname_O>::id,core::Pathname_O::static_classSymbol());
    core::Pathname_O::___staticClass = classcore__Pathname_Oval;
#ifdef USE_MPS
    core::Pathname_O::static_Kind = gctools::GCKind<core::Pathname_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Pathname_Oval,core::Pathname_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Pathname_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Pathname_O>>::allocateClass();
        core::Pathname_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Pathname_O::static_className() % (void*)(core::Pathname_O::static_allocator) );
    classcore__Pathname_Oval->setCreator(core::Pathname_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Pathname_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Pathname_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Pointer_Oval]"));
    core::BuiltInClass_sp classcore__Pointer_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Pointer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Pointer_Oval,_lisp,core::Pointer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Pointer_O>::id,core::Pointer_O::static_classSymbol());
    core::Pointer_O::___staticClass = classcore__Pointer_Oval;
#ifdef USE_MPS
    core::Pointer_O::static_Kind = gctools::GCKind<core::Pointer_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Pointer_Oval,core::Pointer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Pointer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Pointer_O>>::allocateClass();
        core::Pointer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Pointer_O::static_className() % (void*)(core::Pointer_O::static_allocator) );
    classcore__Pointer_Oval->setCreator(core::Pointer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Pointer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Pointer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__PosixTimeDuration_Oval]"));
    core::BuiltInClass_sp classcore__PosixTimeDuration_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__PosixTimeDuration_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__PosixTimeDuration_Oval,_lisp,core::PosixTimeDuration_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::PosixTimeDuration_O>::id,core::PosixTimeDuration_O::static_classSymbol());
    core::PosixTimeDuration_O::___staticClass = classcore__PosixTimeDuration_Oval;
#ifdef USE_MPS
    core::PosixTimeDuration_O::static_Kind = gctools::GCKind<core::PosixTimeDuration_O>::Kind;
#endif
    core::af_setf_findClass(classcore__PosixTimeDuration_Oval,core::PosixTimeDuration_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::PosixTimeDuration_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::PosixTimeDuration_O>>::allocateClass();
        core::PosixTimeDuration_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::PosixTimeDuration_O::static_className() % (void*)(core::PosixTimeDuration_O::static_allocator) );
    classcore__PosixTimeDuration_Oval->setCreator(core::PosixTimeDuration_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::PosixTimeDuration_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__PosixTimeDuration_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__PosixTime_Oval]"));
    core::BuiltInClass_sp classcore__PosixTime_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__PosixTime_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__PosixTime_Oval,_lisp,core::PosixTime_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::PosixTime_O>::id,core::PosixTime_O::static_classSymbol());
    core::PosixTime_O::___staticClass = classcore__PosixTime_Oval;
#ifdef USE_MPS
    core::PosixTime_O::static_Kind = gctools::GCKind<core::PosixTime_O>::Kind;
#endif
    core::af_setf_findClass(classcore__PosixTime_Oval,core::PosixTime_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::PosixTime_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::PosixTime_O>>::allocateClass();
        core::PosixTime_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::PosixTime_O::static_className() % (void*)(core::PosixTime_O::static_allocator) );
    classcore__PosixTime_Oval->setCreator(core::PosixTime_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::PosixTime_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__PosixTime_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__RandomState_Oval]"));
    core::BuiltInClass_sp classcore__RandomState_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__RandomState_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__RandomState_Oval,_lisp,core::RandomState_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::RandomState_O>::id,core::RandomState_O::static_classSymbol());
    core::RandomState_O::___staticClass = classcore__RandomState_Oval;
#ifdef USE_MPS
    core::RandomState_O::static_Kind = gctools::GCKind<core::RandomState_O>::Kind;
#endif
    core::af_setf_findClass(classcore__RandomState_Oval,core::RandomState_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::RandomState_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::RandomState_O>>::allocateClass();
        core::RandomState_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::RandomState_O::static_className() % (void*)(core::RandomState_O::static_allocator) );
    classcore__RandomState_Oval->setCreator(core::RandomState_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::RandomState_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__RandomState_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ReadTable_Oval]"));
    core::BuiltInClass_sp classcore__ReadTable_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ReadTable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ReadTable_Oval,_lisp,core::ReadTable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ReadTable_O>::id,core::ReadTable_O::static_classSymbol());
    core::ReadTable_O::___staticClass = classcore__ReadTable_Oval;
#ifdef USE_MPS
    core::ReadTable_O::static_Kind = gctools::GCKind<core::ReadTable_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ReadTable_Oval,core::ReadTable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ReadTable_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ReadTable_O>>::allocateClass();
        core::ReadTable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ReadTable_O::static_className() % (void*)(core::ReadTable_O::static_allocator) );
    classcore__ReadTable_Oval->setCreator(core::ReadTable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ReadTable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ReadTable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Reader_Oval]"));
    core::BuiltInClass_sp classcore__Reader_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Reader_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Reader_Oval,_lisp,core::Reader_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Reader_O>::id,core::Reader_O::static_classSymbol());
    core::Reader_O::___staticClass = classcore__Reader_Oval;
#ifdef USE_MPS
    core::Reader_O::static_Kind = gctools::GCKind<core::Reader_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Reader_Oval,core::Reader_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Reader_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Reader_O>>::allocateClass();
        core::Reader_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Reader_O::static_className() % (void*)(core::Reader_O::static_allocator) );
    classcore__Reader_Oval->setCreator(core::Reader_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Reader_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Reader_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Record_Oval]"));
    core::BuiltInClass_sp classcore__Record_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Record_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Record_Oval,_lisp,core::Record_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Record_O>::id,core::Record_O::static_classSymbol());
    core::Record_O::___staticClass = classcore__Record_Oval;
#ifdef USE_MPS
    core::Record_O::static_Kind = gctools::GCKind<core::Record_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Record_Oval,core::Record_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Record_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Record_O>>::allocateClass();
        core::Record_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Record_O::static_className() % (void*)(core::Record_O::static_allocator) );
    classcore__Record_Oval->setCreator(core::Record_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Record_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Record_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__RegexMatch_Oval]"));
    core::BuiltInClass_sp classcore__RegexMatch_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__RegexMatch_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__RegexMatch_Oval,_lisp,core::RegexMatch_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::RegexMatch_O>::id,core::RegexMatch_O::static_classSymbol());
    core::RegexMatch_O::___staticClass = classcore__RegexMatch_Oval;
#ifdef USE_MPS
    core::RegexMatch_O::static_Kind = gctools::GCKind<core::RegexMatch_O>::Kind;
#endif
    core::af_setf_findClass(classcore__RegexMatch_Oval,core::RegexMatch_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::RegexMatch_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::RegexMatch_O>>::allocateClass();
        core::RegexMatch_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::RegexMatch_O::static_className() % (void*)(core::RegexMatch_O::static_allocator) );
    classcore__RegexMatch_Oval->setCreator(core::RegexMatch_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::RegexMatch_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__RegexMatch_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Regex_Oval]"));
    core::BuiltInClass_sp classcore__Regex_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Regex_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Regex_Oval,_lisp,core::Regex_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Regex_O>::id,core::Regex_O::static_classSymbol());
    core::Regex_O::___staticClass = classcore__Regex_Oval;
#ifdef USE_MPS
    core::Regex_O::static_Kind = gctools::GCKind<core::Regex_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Regex_Oval,core::Regex_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Regex_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Regex_O>>::allocateClass();
        core::Regex_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Regex_O::static_className() % (void*)(core::Regex_O::static_allocator) );
    classcore__Regex_Oval->setCreator(core::Regex_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Regex_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Regex_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SNode_Oval]"));
    core::BuiltInClass_sp classcore__SNode_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SNode_Oval,_lisp,core::SNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SNode_O>::id,core::SNode_O::static_classSymbol());
    core::SNode_O::___staticClass = classcore__SNode_Oval;
#ifdef USE_MPS
    core::SNode_O::static_Kind = gctools::GCKind<core::SNode_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SNode_Oval,core::SNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SNode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SNode_O>>::allocateClass();
        core::SNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SNode_O::static_className() % (void*)(core::SNode_O::static_allocator) );
    classcore__SNode_Oval->setCreator(core::SNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SingleDispatchMethod_Oval]"));
    core::BuiltInClass_sp classcore__SingleDispatchMethod_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SingleDispatchMethod_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SingleDispatchMethod_Oval,_lisp,core::SingleDispatchMethod_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SingleDispatchMethod_O>::id,core::SingleDispatchMethod_O::static_classSymbol());
    core::SingleDispatchMethod_O::___staticClass = classcore__SingleDispatchMethod_Oval;
#ifdef USE_MPS
    core::SingleDispatchMethod_O::static_Kind = gctools::GCKind<core::SingleDispatchMethod_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SingleDispatchMethod_Oval,core::SingleDispatchMethod_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SingleDispatchMethod_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SingleDispatchMethod_O>>::allocateClass();
        core::SingleDispatchMethod_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SingleDispatchMethod_O::static_className() % (void*)(core::SingleDispatchMethod_O::static_allocator) );
    classcore__SingleDispatchMethod_Oval->setCreator(core::SingleDispatchMethod_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SingleDispatchMethod_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SingleDispatchMethod_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SmallMap_Oval]"));
    core::BuiltInClass_sp classcore__SmallMap_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SmallMap_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SmallMap_Oval,_lisp,core::SmallMap_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SmallMap_O>::id,core::SmallMap_O::static_classSymbol());
    core::SmallMap_O::___staticClass = classcore__SmallMap_Oval;
#ifdef USE_MPS
    core::SmallMap_O::static_Kind = gctools::GCKind<core::SmallMap_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SmallMap_Oval,core::SmallMap_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SmallMap_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SmallMap_O>>::allocateClass();
        core::SmallMap_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SmallMap_O::static_className() % (void*)(core::SmallMap_O::static_allocator) );
    classcore__SmallMap_Oval->setCreator(core::SmallMap_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SmallMap_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SmallMap_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SmallMultimap_Oval]"));
    core::BuiltInClass_sp classcore__SmallMultimap_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SmallMultimap_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SmallMultimap_Oval,_lisp,core::SmallMultimap_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SmallMultimap_O>::id,core::SmallMultimap_O::static_classSymbol());
    core::SmallMultimap_O::___staticClass = classcore__SmallMultimap_Oval;
#ifdef USE_MPS
    core::SmallMultimap_O::static_Kind = gctools::GCKind<core::SmallMultimap_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SmallMultimap_Oval,core::SmallMultimap_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SmallMultimap_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SmallMultimap_O>>::allocateClass();
        core::SmallMultimap_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SmallMultimap_O::static_className() % (void*)(core::SmallMultimap_O::static_allocator) );
    classcore__SmallMultimap_Oval->setCreator(core::SmallMultimap_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SmallMultimap_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SmallMultimap_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SourceFileInfo_Oval]"));
    core::BuiltInClass_sp classcore__SourceFileInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SourceFileInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SourceFileInfo_Oval,_lisp,core::SourceFileInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SourceFileInfo_O>::id,core::SourceFileInfo_O::static_classSymbol());
    core::SourceFileInfo_O::___staticClass = classcore__SourceFileInfo_Oval;
#ifdef USE_MPS
    core::SourceFileInfo_O::static_Kind = gctools::GCKind<core::SourceFileInfo_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SourceFileInfo_Oval,core::SourceFileInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SourceFileInfo_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SourceFileInfo_O>>::allocateClass();
        core::SourceFileInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SourceFileInfo_O::static_className() % (void*)(core::SourceFileInfo_O::static_allocator) );
    classcore__SourceFileInfo_Oval->setCreator(core::SourceFileInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SourceFileInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SourceFileInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SourceManager_Oval]"));
    core::BuiltInClass_sp classcore__SourceManager_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SourceManager_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SourceManager_Oval,_lisp,core::SourceManager_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SourceManager_O>::id,core::SourceManager_O::static_classSymbol());
    core::SourceManager_O::___staticClass = classcore__SourceManager_Oval;
#ifdef USE_MPS
    core::SourceManager_O::static_Kind = gctools::GCKind<core::SourceManager_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SourceManager_Oval,core::SourceManager_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SourceManager_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SourceManager_O>>::allocateClass();
        core::SourceManager_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SourceManager_O::static_className() % (void*)(core::SourceManager_O::static_allocator) );
    classcore__SourceManager_Oval->setCreator(core::SourceManager_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SourceManager_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SourceManager_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SourcePosInfo_Oval]"));
    core::BuiltInClass_sp classcore__SourcePosInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SourcePosInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SourcePosInfo_Oval,_lisp,core::SourcePosInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SourcePosInfo_O>::id,core::SourcePosInfo_O::static_classSymbol());
    core::SourcePosInfo_O::___staticClass = classcore__SourcePosInfo_Oval;
#ifdef USE_MPS
    core::SourcePosInfo_O::static_Kind = gctools::GCKind<core::SourcePosInfo_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SourcePosInfo_Oval,core::SourcePosInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SourcePosInfo_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SourcePosInfo_O>>::allocateClass();
        core::SourcePosInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SourcePosInfo_O::static_className() % (void*)(core::SourcePosInfo_O::static_allocator) );
    classcore__SourcePosInfo_Oval->setCreator(core::SourcePosInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SourcePosInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SourcePosInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StandardObject_Oval]"));
    StandardClass_sp classcore__StandardObject_Oval = StandardClass_O::createUncollectable();
    classcore__StandardObject_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StandardObject_Oval,_lisp,core::StandardObject_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StandardObject_O>::id,core::StandardObject_O::static_classSymbol());
    core::StandardObject_O::___staticClass = classcore__StandardObject_Oval;
#ifdef USE_MPS
    core::StandardObject_O::static_Kind = gctools::GCKind<core::StandardObject_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StandardObject_Oval,core::StandardObject_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StandardObject_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StandardObject_O>>::allocateClass();
        core::StandardObject_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StandardObject_O::static_className() % (void*)(core::StandardObject_O::static_allocator) );
    classcore__StandardObject_Oval->setCreator(core::StandardObject_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StandardObject_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StandardObject_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Stream_Oval]"));
    core::BuiltInClass_sp classcore__Stream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Stream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Stream_Oval,_lisp,core::Stream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Stream_O>::id,core::Stream_O::static_classSymbol());
    core::Stream_O::___staticClass = classcore__Stream_Oval;
#ifdef USE_MPS
    core::Stream_O::static_Kind = gctools::GCKind<core::Stream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Stream_Oval,core::Stream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Stream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Stream_O>>::allocateClass();
        core::Stream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Stream_O::static_className() % (void*)(core::Stream_O::static_allocator) );
    classcore__Stream_Oval->setCreator(core::Stream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Stream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Stream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StructureObject_Oval]"));
    StructureClass_sp classcore__StructureObject_Oval = StructureClass_O::createUncollectable();
    classcore__StructureObject_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StructureObject_Oval,_lisp,core::StructureObject_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StructureObject_O>::id,core::StructureObject_O::static_classSymbol());
    core::StructureObject_O::___staticClass = classcore__StructureObject_Oval;
#ifdef USE_MPS
    core::StructureObject_O::static_Kind = gctools::GCKind<core::StructureObject_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StructureObject_Oval,core::StructureObject_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StructureObject_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StructureObject_O>>::allocateClass();
        core::StructureObject_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StructureObject_O::static_className() % (void*)(core::StructureObject_O::static_allocator) );
    classcore__StructureObject_Oval->setCreator(core::StructureObject_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StructureObject_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StructureObject_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SymbolToEnumConverter_Oval]"));
    core::BuiltInClass_sp classcore__SymbolToEnumConverter_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SymbolToEnumConverter_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SymbolToEnumConverter_Oval,_lisp,core::SymbolToEnumConverter_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SymbolToEnumConverter_O>::id,core::SymbolToEnumConverter_O::static_classSymbol());
    core::SymbolToEnumConverter_O::___staticClass = classcore__SymbolToEnumConverter_Oval;
#ifdef USE_MPS
    core::SymbolToEnumConverter_O::static_Kind = gctools::GCKind<core::SymbolToEnumConverter_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SymbolToEnumConverter_Oval,core::SymbolToEnumConverter_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SymbolToEnumConverter_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SymbolToEnumConverter_O>>::allocateClass();
        core::SymbolToEnumConverter_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SymbolToEnumConverter_O::static_className() % (void*)(core::SymbolToEnumConverter_O::static_allocator) );
    classcore__SymbolToEnumConverter_Oval->setCreator(core::SymbolToEnumConverter_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SymbolToEnumConverter_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SymbolToEnumConverter_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Symbol_Oval]"));
    core::BuiltInClass_sp classcore__Symbol_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Symbol_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Symbol_Oval,_lisp,core::Symbol_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Symbol_O>::id,core::Symbol_O::static_classSymbol());
    core::Symbol_O::___staticClass = classcore__Symbol_Oval;
#ifdef USE_MPS
    core::Symbol_O::static_Kind = gctools::GCKind<core::Symbol_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Symbol_Oval,core::Symbol_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Symbol_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Symbol_O>>::allocateClass();
        core::Symbol_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Symbol_O::static_className() % (void*)(core::Symbol_O::static_allocator) );
    classcore__Symbol_Oval->setCreator(core::Symbol_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Symbol_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Symbol_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__VaList_dummy_Oval]"));
    core::BuiltInClass_sp classcore__VaList_dummy_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__VaList_dummy_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__VaList_dummy_Oval,_lisp,core::VaList_dummy_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::VaList_dummy_O>::id,core::VaList_dummy_O::static_classSymbol());
    core::VaList_dummy_O::___staticClass = classcore__VaList_dummy_Oval;
#ifdef USE_MPS
    core::VaList_dummy_O::static_Kind = gctools::GCKind<core::VaList_dummy_O>::Kind;
#endif
    core::af_setf_findClass(classcore__VaList_dummy_Oval,core::VaList_dummy_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::VaList_dummy_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::VaList_dummy_O>>::allocateClass();
        core::VaList_dummy_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::VaList_dummy_O::static_className() % (void*)(core::VaList_dummy_O::static_allocator) );
    classcore__VaList_dummy_Oval->setCreator(core::VaList_dummy_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::VaList_dummy_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__VaList_dummy_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__WeakHashTable_Oval]"));
    core::BuiltInClass_sp classcore__WeakHashTable_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__WeakHashTable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__WeakHashTable_Oval,_lisp,core::WeakHashTable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::WeakHashTable_O>::id,core::WeakHashTable_O::static_classSymbol());
    core::WeakHashTable_O::___staticClass = classcore__WeakHashTable_Oval;
#ifdef USE_MPS
    core::WeakHashTable_O::static_Kind = gctools::GCKind<core::WeakHashTable_O>::Kind;
#endif
    core::af_setf_findClass(classcore__WeakHashTable_Oval,core::WeakHashTable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::WeakHashTable_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::WeakHashTable_O>>::allocateClass();
        core::WeakHashTable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::WeakHashTable_O::static_className() % (void*)(core::WeakHashTable_O::static_allocator) );
    classcore__WeakHashTable_Oval->setCreator(core::WeakHashTable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::WeakHashTable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__WeakHashTable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__WeakKeyMapping_Oval]"));
    core::BuiltInClass_sp classcore__WeakKeyMapping_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__WeakKeyMapping_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__WeakKeyMapping_Oval,_lisp,core::WeakKeyMapping_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::WeakKeyMapping_O>::id,core::WeakKeyMapping_O::static_classSymbol());
    core::WeakKeyMapping_O::___staticClass = classcore__WeakKeyMapping_Oval;
#ifdef USE_MPS
    core::WeakKeyMapping_O::static_Kind = gctools::GCKind<core::WeakKeyMapping_O>::Kind;
#endif
    core::af_setf_findClass(classcore__WeakKeyMapping_Oval,core::WeakKeyMapping_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::WeakKeyMapping_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::WeakKeyMapping_O>>::allocateClass();
        core::WeakKeyMapping_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::WeakKeyMapping_O::static_className() % (void*)(core::WeakKeyMapping_O::static_allocator) );
    classcore__WeakKeyMapping_Oval->setCreator(core::WeakKeyMapping_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::WeakKeyMapping_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__WeakKeyMapping_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__WeakPointer_Oval]"));
    core::BuiltInClass_sp classcore__WeakPointer_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__WeakPointer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__WeakPointer_Oval,_lisp,core::WeakPointer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::WeakPointer_O>::id,core::WeakPointer_O::static_classSymbol());
    core::WeakPointer_O::___staticClass = classcore__WeakPointer_Oval;
#ifdef USE_MPS
    core::WeakPointer_O::static_Kind = gctools::GCKind<core::WeakPointer_O>::Kind;
#endif
    core::af_setf_findClass(classcore__WeakPointer_Oval,core::WeakPointer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::WeakPointer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::WeakPointer_O>>::allocateClass();
        core::WeakPointer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::WeakPointer_O::static_className() % (void*)(core::WeakPointer_O::static_allocator) );
    classcore__WeakPointer_Oval->setCreator(core::WeakPointer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::WeakPointer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__WeakPointer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__WrappedPointer_Oval]"));
    core::BuiltInClass_sp classcore__WrappedPointer_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__WrappedPointer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__WrappedPointer_Oval,_lisp,core::WrappedPointer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::WrappedPointer_O>::id,core::WrappedPointer_O::static_classSymbol());
    core::WrappedPointer_O::___staticClass = classcore__WrappedPointer_Oval;
#ifdef USE_MPS
    core::WrappedPointer_O::static_Kind = gctools::GCKind<core::WrappedPointer_O>::Kind;
#endif
    core::af_setf_findClass(classcore__WrappedPointer_Oval,core::WrappedPointer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::WrappedPointer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::WrappedPointer_O>>::allocateClass();
        core::WrappedPointer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::WrappedPointer_O::static_className() % (void*)(core::WrappedPointer_O::static_allocator) );
    classcore__WrappedPointer_Oval->setCreator(core::WrappedPointer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::WrappedPointer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__WrappedPointer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__AttributeSet_Oval]"));
    core::BuiltInClass_sp classllvmo__AttributeSet_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__AttributeSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__AttributeSet_Oval,_lisp,llvmo::AttributeSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::AttributeSet_O>::id,llvmo::AttributeSet_O::static_classSymbol());
    llvmo::AttributeSet_O::___staticClass = classllvmo__AttributeSet_Oval;
#ifdef USE_MPS
    llvmo::AttributeSet_O::static_Kind = gctools::GCKind<llvmo::AttributeSet_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__AttributeSet_Oval,llvmo::AttributeSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::AttributeSet_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::AttributeSet_O>>::allocateClass();
        llvmo::AttributeSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::AttributeSet_O::static_className() % (void*)(llvmo::AttributeSet_O::static_allocator) );
    classllvmo__AttributeSet_Oval->setCreator(llvmo::AttributeSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::AttributeSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__AttributeSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Attribute_Oval]"));
    core::BuiltInClass_sp classllvmo__Attribute_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Attribute_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Attribute_Oval,_lisp,llvmo::Attribute_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Attribute_O>::id,llvmo::Attribute_O::static_classSymbol());
    llvmo::Attribute_O::___staticClass = classllvmo__Attribute_Oval;
#ifdef USE_MPS
    llvmo::Attribute_O::static_Kind = gctools::GCKind<llvmo::Attribute_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Attribute_Oval,llvmo::Attribute_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Attribute_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Attribute_O>>::allocateClass();
        llvmo::Attribute_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Attribute_O::static_className() % (void*)(llvmo::Attribute_O::static_allocator) );
    classllvmo__Attribute_Oval->setCreator(llvmo::Attribute_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Attribute_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Attribute_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DebugInfo_Oval]"));
    core::BuiltInClass_sp classllvmo__DebugInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DebugInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DebugInfo_Oval,_lisp,llvmo::DebugInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DebugInfo_O>::id,llvmo::DebugInfo_O::static_classSymbol());
    llvmo::DebugInfo_O::___staticClass = classllvmo__DebugInfo_Oval;
#ifdef USE_MPS
    llvmo::DebugInfo_O::static_Kind = gctools::GCKind<llvmo::DebugInfo_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DebugInfo_Oval,llvmo::DebugInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DebugInfo_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DebugInfo_O>>::allocateClass();
        llvmo::DebugInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DebugInfo_O::static_className() % (void*)(llvmo::DebugInfo_O::static_allocator) );
    classllvmo__DebugInfo_Oval->setCreator(llvmo::DebugInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DebugInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DebugInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DebugLoc_Oval]"));
    core::BuiltInClass_sp classllvmo__DebugLoc_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DebugLoc_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DebugLoc_Oval,_lisp,llvmo::DebugLoc_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DebugLoc_O>::id,llvmo::DebugLoc_O::static_classSymbol());
    llvmo::DebugLoc_O::___staticClass = classllvmo__DebugLoc_Oval;
#ifdef USE_MPS
    llvmo::DebugLoc_O::static_Kind = gctools::GCKind<llvmo::DebugLoc_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DebugLoc_Oval,llvmo::DebugLoc_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DebugLoc_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DebugLoc_O>>::allocateClass();
        llvmo::DebugLoc_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DebugLoc_O::static_className() % (void*)(llvmo::DebugLoc_O::static_allocator) );
    classllvmo__DebugLoc_Oval->setCreator(llvmo::DebugLoc_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DebugLoc_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DebugLoc_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__InsertPoint_Oval]"));
    core::BuiltInClass_sp classllvmo__InsertPoint_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__InsertPoint_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__InsertPoint_Oval,_lisp,llvmo::InsertPoint_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::InsertPoint_O>::id,llvmo::InsertPoint_O::static_classSymbol());
    llvmo::InsertPoint_O::___staticClass = classllvmo__InsertPoint_Oval;
#ifdef USE_MPS
    llvmo::InsertPoint_O::static_Kind = gctools::GCKind<llvmo::InsertPoint_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__InsertPoint_Oval,llvmo::InsertPoint_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::InsertPoint_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::InsertPoint_O>>::allocateClass();
        llvmo::InsertPoint_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::InsertPoint_O::static_className() % (void*)(llvmo::InsertPoint_O::static_allocator) );
    classllvmo__InsertPoint_Oval->setCreator(llvmo::InsertPoint_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::InsertPoint_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__InsertPoint_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__IterateCons_Oval]"));
    core::BuiltInClass_sp classadapt__IterateCons_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__IterateCons_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__IterateCons_Oval,_lisp,adapt::IterateCons_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::IterateCons_O>::id,adapt::IterateCons_O::static_classSymbol());
    adapt::IterateCons_O::___staticClass = classadapt__IterateCons_Oval;
#ifdef USE_MPS
    adapt::IterateCons_O::static_Kind = gctools::GCKind<adapt::IterateCons_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__IterateCons_Oval,adapt::IterateCons_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::IterateCons_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::IterateCons_O>>::allocateClass();
        adapt::IterateCons_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::IterateCons_O::static_className() % (void*)(adapt::IterateCons_O::static_allocator) );
    classadapt__IterateCons_Oval->setCreator(adapt::IterateCons_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::IterateCons_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__IterateCons_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__QDomNode_Oval]"));
    core::BuiltInClass_sp classadapt__QDomNode_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__QDomNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__QDomNode_Oval,_lisp,adapt::QDomNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::QDomNode_O>::id,adapt::QDomNode_O::static_classSymbol());
    adapt::QDomNode_O::___staticClass = classadapt__QDomNode_Oval;
#ifdef USE_MPS
    adapt::QDomNode_O::static_Kind = gctools::GCKind<adapt::QDomNode_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__QDomNode_Oval,adapt::QDomNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::QDomNode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::QDomNode_O>>::allocateClass();
        adapt::QDomNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::QDomNode_O::static_className() % (void*)(adapt::QDomNode_O::static_allocator) );
    classadapt__QDomNode_Oval->setCreator(adapt::QDomNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::QDomNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__QDomNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__StringSet_Oval]"));
    core::BuiltInClass_sp classadapt__StringSet_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__StringSet_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__StringSet_Oval,_lisp,adapt::StringSet_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::StringSet_O>::id,adapt::StringSet_O::static_classSymbol());
    adapt::StringSet_O::___staticClass = classadapt__StringSet_Oval;
#ifdef USE_MPS
    adapt::StringSet_O::static_Kind = gctools::GCKind<adapt::StringSet_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__StringSet_Oval,adapt::StringSet_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::StringSet_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::StringSet_O>>::allocateClass();
        adapt::StringSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::StringSet_O::static_className() % (void*)(adapt::StringSet_O::static_allocator) );
    classadapt__StringSet_Oval->setCreator(adapt::StringSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::StringSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__StringSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_AdaptPkg
    LOG(BF("Creating class[classadapt__SymbolList_Oval]"));
    core::BuiltInClass_sp classadapt__SymbolList_Oval = core::BuiltInClass_O::createUncollectable();
    classadapt__SymbolList_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classadapt__SymbolList_Oval,_lisp,adapt::SymbolList_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<adapt::SymbolList_O>::id,adapt::SymbolList_O::static_classSymbol());
    adapt::SymbolList_O::___staticClass = classadapt__SymbolList_Oval;
#ifdef USE_MPS
    adapt::SymbolList_O::static_Kind = gctools::GCKind<adapt::SymbolList_O>::Kind;
#endif
    core::af_setf_findClass(classadapt__SymbolList_Oval,adapt::SymbolList_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<adapt::SymbolList_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::SymbolList_O>>::allocateClass();
        adapt::SymbolList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::SymbolList_O::static_className() % (void*)(adapt::SymbolList_O::static_allocator) );
    classadapt__SymbolList_Oval->setCreator(adapt::SymbolList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::SymbolList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__SymbolList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ActivationFrame_Oval]"));
    core::BuiltInClass_sp classcore__ActivationFrame_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ActivationFrame_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ActivationFrame_Oval,_lisp,core::ActivationFrame_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ActivationFrame_O>::id,core::ActivationFrame_O::static_classSymbol());
    core::ActivationFrame_O::___staticClass = classcore__ActivationFrame_Oval;
#ifdef USE_MPS
    core::ActivationFrame_O::static_Kind = gctools::GCKind<core::ActivationFrame_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ActivationFrame_Oval,core::ActivationFrame_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ActivationFrame_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ActivationFrame_O>>::allocateClass();
        core::ActivationFrame_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ActivationFrame_O::static_className() % (void*)(core::ActivationFrame_O::static_allocator) );
    classcore__ActivationFrame_Oval->setCreator(core::ActivationFrame_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ActivationFrame_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ActivationFrame_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__AnsiStream_Oval]"));
    core::BuiltInClass_sp classcore__AnsiStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__AnsiStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__AnsiStream_Oval,_lisp,core::AnsiStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::AnsiStream_O>::id,core::AnsiStream_O::static_classSymbol());
    core::AnsiStream_O::___staticClass = classcore__AnsiStream_Oval;
#ifdef USE_MPS
    core::AnsiStream_O::static_Kind = gctools::GCKind<core::AnsiStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__AnsiStream_Oval,core::AnsiStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::AnsiStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::AnsiStream_O>>::allocateClass();
        core::AnsiStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::AnsiStream_O::static_className() % (void*)(core::AnsiStream_O::static_allocator) );
    classcore__AnsiStream_Oval->setCreator(core::AnsiStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::AnsiStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__AnsiStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ArrayDisplaced_Oval]"));
    core::BuiltInClass_sp classcore__ArrayDisplaced_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ArrayDisplaced_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ArrayDisplaced_Oval,_lisp,core::ArrayDisplaced_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ArrayDisplaced_O>::id,core::ArrayDisplaced_O::static_classSymbol());
    core::ArrayDisplaced_O::___staticClass = classcore__ArrayDisplaced_Oval;
#ifdef USE_MPS
    core::ArrayDisplaced_O::static_Kind = gctools::GCKind<core::ArrayDisplaced_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ArrayDisplaced_Oval,core::ArrayDisplaced_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ArrayDisplaced_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ArrayDisplaced_O>>::allocateClass();
        core::ArrayDisplaced_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ArrayDisplaced_O::static_className() % (void*)(core::ArrayDisplaced_O::static_allocator) );
    classcore__ArrayDisplaced_Oval->setCreator(core::ArrayDisplaced_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ArrayDisplaced_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ArrayDisplaced_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ArrayObjects_Oval]"));
    core::BuiltInClass_sp classcore__ArrayObjects_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ArrayObjects_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ArrayObjects_Oval,_lisp,core::ArrayObjects_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ArrayObjects_O>::id,core::ArrayObjects_O::static_classSymbol());
    core::ArrayObjects_O::___staticClass = classcore__ArrayObjects_Oval;
#ifdef USE_MPS
    core::ArrayObjects_O::static_Kind = gctools::GCKind<core::ArrayObjects_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ArrayObjects_Oval,core::ArrayObjects_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ArrayObjects_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ArrayObjects_O>>::allocateClass();
        core::ArrayObjects_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ArrayObjects_O::static_className() % (void*)(core::ArrayObjects_O::static_allocator) );
    classcore__ArrayObjects_Oval->setCreator(core::ArrayObjects_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ArrayObjects_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ArrayObjects_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__BranchSNode_Oval]"));
    core::BuiltInClass_sp classcore__BranchSNode_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__BranchSNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__BranchSNode_Oval,_lisp,core::BranchSNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::BranchSNode_O>::id,core::BranchSNode_O::static_classSymbol());
    core::BranchSNode_O::___staticClass = classcore__BranchSNode_Oval;
#ifdef USE_MPS
    core::BranchSNode_O::static_Kind = gctools::GCKind<core::BranchSNode_O>::Kind;
#endif
    core::af_setf_findClass(classcore__BranchSNode_Oval,core::BranchSNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::BranchSNode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::BranchSNode_O>>::allocateClass();
        core::BranchSNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::BranchSNode_O::static_className() % (void*)(core::BranchSNode_O::static_allocator) );
    classcore__BranchSNode_Oval->setCreator(core::BranchSNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::BranchSNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__BranchSNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__CompiledFunction_Oval]"));
    core::BuiltInClass_sp classcore__CompiledFunction_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__CompiledFunction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__CompiledFunction_Oval,_lisp,core::CompiledFunction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::CompiledFunction_O>::id,core::CompiledFunction_O::static_classSymbol());
    core::CompiledFunction_O::___staticClass = classcore__CompiledFunction_Oval;
#ifdef USE_MPS
    core::CompiledFunction_O::static_Kind = gctools::GCKind<core::CompiledFunction_O>::Kind;
#endif
    core::af_setf_findClass(classcore__CompiledFunction_Oval,core::CompiledFunction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::CompiledFunction_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::CompiledFunction_O>>::allocateClass();
        core::CompiledFunction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::CompiledFunction_O::static_className() % (void*)(core::CompiledFunction_O::static_allocator) );
    classcore__CompiledFunction_Oval->setCreator(core::CompiledFunction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::CompiledFunction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__CompiledFunction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Complex_Oval]"));
    core::BuiltInClass_sp classcore__Complex_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Complex_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Complex_Oval,_lisp,core::Complex_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Complex_O>::id,core::Complex_O::static_classSymbol());
    core::Complex_O::___staticClass = classcore__Complex_Oval;
#ifdef USE_MPS
    core::Complex_O::static_Kind = gctools::GCKind<core::Complex_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Complex_Oval,core::Complex_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Complex_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Complex_O>>::allocateClass();
        core::Complex_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Complex_O::static_className() % (void*)(core::Complex_O::static_allocator) );
    classcore__Complex_Oval->setCreator(core::Complex_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Complex_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Complex_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__DirectoryIterator_Oval]"));
    core::BuiltInClass_sp classcore__DirectoryIterator_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__DirectoryIterator_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__DirectoryIterator_Oval,_lisp,core::DirectoryIterator_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::DirectoryIterator_O>::id,core::DirectoryIterator_O::static_classSymbol());
    core::DirectoryIterator_O::___staticClass = classcore__DirectoryIterator_Oval;
#ifdef USE_MPS
    core::DirectoryIterator_O::static_Kind = gctools::GCKind<core::DirectoryIterator_O>::Kind;
#endif
    core::af_setf_findClass(classcore__DirectoryIterator_Oval,core::DirectoryIterator_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::DirectoryIterator_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::DirectoryIterator_O>>::allocateClass();
        core::DirectoryIterator_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::DirectoryIterator_O::static_className() % (void*)(core::DirectoryIterator_O::static_allocator) );
    classcore__DirectoryIterator_Oval->setCreator(core::DirectoryIterator_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::DirectoryIterator_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__DirectoryIterator_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ForeignData_Oval]"));
    core::BuiltInClass_sp classcore__ForeignData_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ForeignData_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ForeignData_Oval,_lisp,core::ForeignData_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ForeignData_O>::id,core::ForeignData_O::static_classSymbol());
    core::ForeignData_O::___staticClass = classcore__ForeignData_Oval;
#ifdef USE_MPS
    core::ForeignData_O::static_Kind = gctools::GCKind<core::ForeignData_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ForeignData_Oval,core::ForeignData_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ForeignData_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ForeignData_O>>::allocateClass();
        core::ForeignData_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ForeignData_O::static_className() % (void*)(core::ForeignData_O::static_allocator) );
    classcore__ForeignData_Oval->setCreator(core::ForeignData_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ForeignData_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ForeignData_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__GlueEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__GlueEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__GlueEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__GlueEnvironment_Oval,_lisp,core::GlueEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::GlueEnvironment_O>::id,core::GlueEnvironment_O::static_classSymbol());
    core::GlueEnvironment_O::___staticClass = classcore__GlueEnvironment_Oval;
#ifdef USE_MPS
    core::GlueEnvironment_O::static_Kind = gctools::GCKind<core::GlueEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__GlueEnvironment_Oval,core::GlueEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::GlueEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::GlueEnvironment_O>>::allocateClass();
        core::GlueEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::GlueEnvironment_O::static_className() % (void*)(core::GlueEnvironment_O::static_allocator) );
    classcore__GlueEnvironment_Oval->setCreator(core::GlueEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::GlueEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__GlueEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__HashTableEq_Oval]"));
    core::BuiltInClass_sp classcore__HashTableEq_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__HashTableEq_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__HashTableEq_Oval,_lisp,core::HashTableEq_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::HashTableEq_O>::id,core::HashTableEq_O::static_classSymbol());
    core::HashTableEq_O::___staticClass = classcore__HashTableEq_Oval;
#ifdef USE_MPS
    core::HashTableEq_O::static_Kind = gctools::GCKind<core::HashTableEq_O>::Kind;
#endif
    core::af_setf_findClass(classcore__HashTableEq_Oval,core::HashTableEq_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::HashTableEq_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::HashTableEq_O>>::allocateClass();
        core::HashTableEq_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::HashTableEq_O::static_className() % (void*)(core::HashTableEq_O::static_allocator) );
    classcore__HashTableEq_Oval->setCreator(core::HashTableEq_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::HashTableEq_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__HashTableEq_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__HashTableEql_Oval]"));
    core::BuiltInClass_sp classcore__HashTableEql_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__HashTableEql_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__HashTableEql_Oval,_lisp,core::HashTableEql_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::HashTableEql_O>::id,core::HashTableEql_O::static_classSymbol());
    core::HashTableEql_O::___staticClass = classcore__HashTableEql_Oval;
#ifdef USE_MPS
    core::HashTableEql_O::static_Kind = gctools::GCKind<core::HashTableEql_O>::Kind;
#endif
    core::af_setf_findClass(classcore__HashTableEql_Oval,core::HashTableEql_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::HashTableEql_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::HashTableEql_O>>::allocateClass();
        core::HashTableEql_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::HashTableEql_O::static_className() % (void*)(core::HashTableEql_O::static_allocator) );
    classcore__HashTableEql_Oval->setCreator(core::HashTableEql_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::HashTableEql_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__HashTableEql_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__HashTableEqual_Oval]"));
    core::BuiltInClass_sp classcore__HashTableEqual_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__HashTableEqual_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__HashTableEqual_Oval,_lisp,core::HashTableEqual_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::HashTableEqual_O>::id,core::HashTableEqual_O::static_classSymbol());
    core::HashTableEqual_O::___staticClass = classcore__HashTableEqual_Oval;
#ifdef USE_MPS
    core::HashTableEqual_O::static_Kind = gctools::GCKind<core::HashTableEqual_O>::Kind;
#endif
    core::af_setf_findClass(classcore__HashTableEqual_Oval,core::HashTableEqual_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::HashTableEqual_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::HashTableEqual_O>>::allocateClass();
        core::HashTableEqual_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::HashTableEqual_O::static_className() % (void*)(core::HashTableEqual_O::static_allocator) );
    classcore__HashTableEqual_Oval->setCreator(core::HashTableEqual_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::HashTableEqual_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__HashTableEqual_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__HashTableEqualp_Oval]"));
    core::BuiltInClass_sp classcore__HashTableEqualp_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__HashTableEqualp_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__HashTableEqualp_Oval,_lisp,core::HashTableEqualp_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::HashTableEqualp_O>::id,core::HashTableEqualp_O::static_classSymbol());
    core::HashTableEqualp_O::___staticClass = classcore__HashTableEqualp_Oval;
#ifdef USE_MPS
    core::HashTableEqualp_O::static_Kind = gctools::GCKind<core::HashTableEqualp_O>::Kind;
#endif
    core::af_setf_findClass(classcore__HashTableEqualp_Oval,core::HashTableEqualp_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::HashTableEqualp_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::HashTableEqualp_O>>::allocateClass();
        core::HashTableEqualp_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::HashTableEqualp_O::static_className() % (void*)(core::HashTableEqualp_O::static_allocator) );
    classcore__HashTableEqualp_Oval->setCreator(core::HashTableEqualp_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::HashTableEqualp_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__HashTableEqualp_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Instance_Oval]"));
    core::BuiltInClass_sp classcore__Instance_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Instance_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Instance_Oval,_lisp,core::Instance_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Instance_O>::id,core::Instance_O::static_classSymbol());
    core::Instance_O::___staticClass = classcore__Instance_Oval;
#ifdef USE_MPS
    core::Instance_O::static_Kind = gctools::GCKind<core::Instance_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Instance_Oval,core::Instance_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Instance_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Instance_O>>::allocateClass();
        core::Instance_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Instance_O::static_className() % (void*)(core::Instance_O::static_allocator) );
    classcore__Instance_Oval->setCreator(core::Instance_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Instance_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Instance_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LeafSNode_Oval]"));
    core::BuiltInClass_sp classcore__LeafSNode_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LeafSNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LeafSNode_Oval,_lisp,core::LeafSNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LeafSNode_O>::id,core::LeafSNode_O::static_classSymbol());
    core::LeafSNode_O::___staticClass = classcore__LeafSNode_Oval;
#ifdef USE_MPS
    core::LeafSNode_O::static_Kind = gctools::GCKind<core::LeafSNode_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LeafSNode_Oval,core::LeafSNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LeafSNode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LeafSNode_O>>::allocateClass();
        core::LeafSNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LeafSNode_O::static_className() % (void*)(core::LeafSNode_O::static_allocator) );
    classcore__LeafSNode_Oval->setCreator(core::LeafSNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LeafSNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LeafSNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LexicalEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__LexicalEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LexicalEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LexicalEnvironment_Oval,_lisp,core::LexicalEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LexicalEnvironment_O>::id,core::LexicalEnvironment_O::static_classSymbol());
    core::LexicalEnvironment_O::___staticClass = classcore__LexicalEnvironment_Oval;
#ifdef USE_MPS
    core::LexicalEnvironment_O::static_Kind = gctools::GCKind<core::LexicalEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LexicalEnvironment_Oval,core::LexicalEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LexicalEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LexicalEnvironment_O>>::allocateClass();
        core::LexicalEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LexicalEnvironment_O::static_className() % (void*)(core::LexicalEnvironment_O::static_allocator) );
    classcore__LexicalEnvironment_Oval->setCreator(core::LexicalEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LexicalEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LexicalEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LoadArchive_Oval]"));
    core::BuiltInClass_sp classcore__LoadArchive_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LoadArchive_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LoadArchive_Oval,_lisp,core::LoadArchive_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LoadArchive_O>::id,core::LoadArchive_O::static_classSymbol());
    core::LoadArchive_O::___staticClass = classcore__LoadArchive_Oval;
#ifdef USE_MPS
    core::LoadArchive_O::static_Kind = gctools::GCKind<core::LoadArchive_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LoadArchive_Oval,core::LoadArchive_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LoadArchive_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LoadArchive_O>>::allocateClass();
        core::LoadArchive_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LoadArchive_O::static_className() % (void*)(core::LoadArchive_O::static_allocator) );
    classcore__LoadArchive_Oval->setCreator(core::LoadArchive_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LoadArchive_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LoadArchive_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LogicalPathname_Oval]"));
    core::BuiltInClass_sp classcore__LogicalPathname_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LogicalPathname_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LogicalPathname_Oval,_lisp,core::LogicalPathname_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LogicalPathname_O>::id,core::LogicalPathname_O::static_classSymbol());
    core::LogicalPathname_O::___staticClass = classcore__LogicalPathname_Oval;
#ifdef USE_MPS
    core::LogicalPathname_O::static_Kind = gctools::GCKind<core::LogicalPathname_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LogicalPathname_Oval,core::LogicalPathname_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LogicalPathname_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LogicalPathname_O>>::allocateClass();
        core::LogicalPathname_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LogicalPathname_O::static_className() % (void*)(core::LogicalPathname_O::static_allocator) );
    classcore__LogicalPathname_Oval->setCreator(core::LogicalPathname_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LogicalPathname_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LogicalPathname_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Metaobject_Oval]"));
    StandardClass_sp classcore__Metaobject_Oval = StandardClass_O::createUncollectable();
    classcore__Metaobject_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Metaobject_Oval,_lisp,core::Metaobject_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Metaobject_O>::id,core::Metaobject_O::static_classSymbol());
    core::Metaobject_O::___staticClass = classcore__Metaobject_Oval;
#ifdef USE_MPS
    core::Metaobject_O::static_Kind = gctools::GCKind<core::Metaobject_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Metaobject_Oval,core::Metaobject_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Metaobject_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Metaobject_O>>::allocateClass();
        core::Metaobject_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Metaobject_O::static_className() % (void*)(core::Metaobject_O::static_allocator) );
    classcore__Metaobject_Oval->setCreator(core::Metaobject_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Metaobject_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Metaobject_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Null_Oval]"));
    core::BuiltInClass_sp classcore__Null_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Null_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Null_Oval,_lisp,core::Null_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Null_O>::id,core::Null_O::static_classSymbol());
    core::Null_O::___staticClass = classcore__Null_Oval;
#ifdef USE_MPS
    core::Null_O::static_Kind = gctools::GCKind<core::Null_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Null_Oval,core::Null_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Null_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Null_O>>::allocateClass();
        core::Null_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Null_O::static_className() % (void*)(core::Null_O::static_allocator) );
    classcore__Null_Oval->setCreator(core::Null_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Null_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Null_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Real_Oval]"));
    core::BuiltInClass_sp classcore__Real_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Real_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Real_Oval,_lisp,core::Real_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Real_O>::id,core::Real_O::static_classSymbol());
    core::Real_O::___staticClass = classcore__Real_Oval;
#ifdef USE_MPS
    core::Real_O::static_Kind = gctools::GCKind<core::Real_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Real_Oval,core::Real_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Real_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Real_O>>::allocateClass();
        core::Real_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Real_O::static_className() % (void*)(core::Real_O::static_allocator) );
    classcore__Real_Oval->setCreator(core::Real_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Real_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Real_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__RecursiveDirectoryIterator_Oval]"));
    core::BuiltInClass_sp classcore__RecursiveDirectoryIterator_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__RecursiveDirectoryIterator_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__RecursiveDirectoryIterator_Oval,_lisp,core::RecursiveDirectoryIterator_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::RecursiveDirectoryIterator_O>::id,core::RecursiveDirectoryIterator_O::static_classSymbol());
    core::RecursiveDirectoryIterator_O::___staticClass = classcore__RecursiveDirectoryIterator_Oval;
#ifdef USE_MPS
    core::RecursiveDirectoryIterator_O::static_Kind = gctools::GCKind<core::RecursiveDirectoryIterator_O>::Kind;
#endif
    core::af_setf_findClass(classcore__RecursiveDirectoryIterator_Oval,core::RecursiveDirectoryIterator_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::RecursiveDirectoryIterator_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::RecursiveDirectoryIterator_O>>::allocateClass();
        core::RecursiveDirectoryIterator_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::RecursiveDirectoryIterator_O::static_className() % (void*)(core::RecursiveDirectoryIterator_O::static_allocator) );
    classcore__RecursiveDirectoryIterator_Oval->setCreator(core::RecursiveDirectoryIterator_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::RecursiveDirectoryIterator_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__RecursiveDirectoryIterator_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SaveArchive_Oval]"));
    core::BuiltInClass_sp classcore__SaveArchive_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SaveArchive_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SaveArchive_Oval,_lisp,core::SaveArchive_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SaveArchive_O>::id,core::SaveArchive_O::static_classSymbol());
    core::SaveArchive_O::___staticClass = classcore__SaveArchive_Oval;
#ifdef USE_MPS
    core::SaveArchive_O::static_Kind = gctools::GCKind<core::SaveArchive_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SaveArchive_Oval,core::SaveArchive_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SaveArchive_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SaveArchive_O>>::allocateClass();
        core::SaveArchive_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SaveArchive_O::static_className() % (void*)(core::SaveArchive_O::static_allocator) );
    classcore__SaveArchive_Oval->setCreator(core::SaveArchive_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SaveArchive_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SaveArchive_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SingleDispatchEffectiveMethodFunction_Oval]"));
    core::BuiltInClass_sp classcore__SingleDispatchEffectiveMethodFunction_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SingleDispatchEffectiveMethodFunction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SingleDispatchEffectiveMethodFunction_Oval,_lisp,core::SingleDispatchEffectiveMethodFunction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SingleDispatchEffectiveMethodFunction_O>::id,core::SingleDispatchEffectiveMethodFunction_O::static_classSymbol());
    core::SingleDispatchEffectiveMethodFunction_O::___staticClass = classcore__SingleDispatchEffectiveMethodFunction_Oval;
#ifdef USE_MPS
    core::SingleDispatchEffectiveMethodFunction_O::static_Kind = gctools::GCKind<core::SingleDispatchEffectiveMethodFunction_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SingleDispatchEffectiveMethodFunction_Oval,core::SingleDispatchEffectiveMethodFunction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SingleDispatchEffectiveMethodFunction_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SingleDispatchEffectiveMethodFunction_O>>::allocateClass();
        core::SingleDispatchEffectiveMethodFunction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SingleDispatchEffectiveMethodFunction_O::static_className() % (void*)(core::SingleDispatchEffectiveMethodFunction_O::static_allocator) );
    classcore__SingleDispatchEffectiveMethodFunction_Oval->setCreator(core::SingleDispatchEffectiveMethodFunction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SingleDispatchEffectiveMethodFunction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SingleDispatchEffectiveMethodFunction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SingleDispatchGenericFunction_Oval]"));
    core::BuiltInClass_sp classcore__SingleDispatchGenericFunction_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SingleDispatchGenericFunction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SingleDispatchGenericFunction_Oval,_lisp,core::SingleDispatchGenericFunction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SingleDispatchGenericFunction_O>::id,core::SingleDispatchGenericFunction_O::static_classSymbol());
    core::SingleDispatchGenericFunction_O::___staticClass = classcore__SingleDispatchGenericFunction_Oval;
#ifdef USE_MPS
    core::SingleDispatchGenericFunction_O::static_Kind = gctools::GCKind<core::SingleDispatchGenericFunction_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SingleDispatchGenericFunction_Oval,core::SingleDispatchGenericFunction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SingleDispatchGenericFunction_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SingleDispatchGenericFunction_O>>::allocateClass();
        core::SingleDispatchGenericFunction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SingleDispatchGenericFunction_O::static_className() % (void*)(core::SingleDispatchGenericFunction_O::static_allocator) );
    classcore__SingleDispatchGenericFunction_Oval->setCreator(core::SingleDispatchGenericFunction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SingleDispatchGenericFunction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SingleDispatchGenericFunction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SpecialForm_Oval]"));
    core::BuiltInClass_sp classcore__SpecialForm_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SpecialForm_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SpecialForm_Oval,_lisp,core::SpecialForm_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SpecialForm_O>::id,core::SpecialForm_O::static_classSymbol());
    core::SpecialForm_O::___staticClass = classcore__SpecialForm_Oval;
#ifdef USE_MPS
    core::SpecialForm_O::static_Kind = gctools::GCKind<core::SpecialForm_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SpecialForm_Oval,core::SpecialForm_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SpecialForm_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SpecialForm_O>>::allocateClass();
        core::SpecialForm_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SpecialForm_O::static_className() % (void*)(core::SpecialForm_O::static_allocator) );
    classcore__SpecialForm_Oval->setCreator(core::SpecialForm_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SpecialForm_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SpecialForm_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__UserData_Oval]"));
    core::BuiltInClass_sp classcore__UserData_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__UserData_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__UserData_Oval,_lisp,core::UserData_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::UserData_O>::id,core::UserData_O::static_classSymbol());
    core::UserData_O::___staticClass = classcore__UserData_Oval;
#ifdef USE_MPS
    core::UserData_O::static_Kind = gctools::GCKind<core::UserData_O>::Kind;
#endif
    core::af_setf_findClass(classcore__UserData_Oval,core::UserData_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::UserData_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::UserData_O>>::allocateClass();
        core::UserData_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::UserData_O::static_className() % (void*)(core::UserData_O::static_allocator) );
    classcore__UserData_Oval->setCreator(core::UserData_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::UserData_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__UserData_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Vector_Oval]"));
    core::BuiltInClass_sp classcore__Vector_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Vector_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Vector_Oval,_lisp,core::Vector_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Vector_O>::id,core::Vector_O::static_classSymbol());
    core::Vector_O::___staticClass = classcore__Vector_Oval;
#ifdef USE_MPS
    core::Vector_O::static_Kind = gctools::GCKind<core::Vector_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Vector_Oval,core::Vector_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Vector_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Vector_O>>::allocateClass();
        core::Vector_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Vector_O::static_className() % (void*)(core::Vector_O::static_allocator) );
    classcore__Vector_Oval->setCreator(core::Vector_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Vector_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Vector_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__WeakKeyHashTable_Oval]"));
    core::BuiltInClass_sp classcore__WeakKeyHashTable_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__WeakKeyHashTable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__WeakKeyHashTable_Oval,_lisp,core::WeakKeyHashTable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::WeakKeyHashTable_O>::id,core::WeakKeyHashTable_O::static_classSymbol());
    core::WeakKeyHashTable_O::___staticClass = classcore__WeakKeyHashTable_Oval;
#ifdef USE_MPS
    core::WeakKeyHashTable_O::static_Kind = gctools::GCKind<core::WeakKeyHashTable_O>::Kind;
#endif
    core::af_setf_findClass(classcore__WeakKeyHashTable_Oval,core::WeakKeyHashTable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::WeakKeyHashTable_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::WeakKeyHashTable_O>>::allocateClass();
        core::WeakKeyHashTable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::WeakKeyHashTable_O::static_className() % (void*)(core::WeakKeyHashTable_O::static_allocator) );
    classcore__WeakKeyHashTable_Oval->setCreator(core::WeakKeyHashTable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::WeakKeyHashTable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__WeakKeyHashTable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__BoundingBox_Oval]"));
    core::BuiltInClass_sp classgeom__BoundingBox_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__BoundingBox_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__BoundingBox_Oval,_lisp,geom::BoundingBox_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::BoundingBox_O>::id,geom::BoundingBox_O::static_classSymbol());
    geom::BoundingBox_O::___staticClass = classgeom__BoundingBox_Oval;
#ifdef USE_MPS
    geom::BoundingBox_O::static_Kind = gctools::GCKind<geom::BoundingBox_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__BoundingBox_Oval,geom::BoundingBox_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::BoundingBox_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::BoundingBox_O>>::allocateClass();
        geom::BoundingBox_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::BoundingBox_O::static_className() % (void*)(geom::BoundingBox_O::static_allocator) );
    classgeom__BoundingBox_Oval->setCreator(geom::BoundingBox_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::BoundingBox_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__BoundingBox_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__Color_Oval]"));
    core::BuiltInClass_sp classgeom__Color_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__Color_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__Color_Oval,_lisp,geom::Color_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::Color_O>::id,geom::Color_O::static_classSymbol());
    geom::Color_O::___staticClass = classgeom__Color_Oval;
#ifdef USE_MPS
    geom::Color_O::static_Kind = gctools::GCKind<geom::Color_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__Color_Oval,geom::Color_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::Color_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::Color_O>>::allocateClass();
        geom::Color_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::Color_O::static_className() % (void*)(geom::Color_O::static_allocator) );
    classgeom__Color_Oval->setCreator(geom::Color_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::Color_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__Color_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__OMatrix_Oval]"));
    core::BuiltInClass_sp classgeom__OMatrix_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__OMatrix_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__OMatrix_Oval,_lisp,geom::OMatrix_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::OMatrix_O>::id,geom::OMatrix_O::static_classSymbol());
    geom::OMatrix_O::___staticClass = classgeom__OMatrix_Oval;
#ifdef USE_MPS
    geom::OMatrix_O::static_Kind = gctools::GCKind<geom::OMatrix_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__OMatrix_Oval,geom::OMatrix_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::OMatrix_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::OMatrix_O>>::allocateClass();
        geom::OMatrix_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::OMatrix_O::static_className() % (void*)(geom::OMatrix_O::static_allocator) );
    classgeom__OMatrix_Oval->setCreator(geom::OMatrix_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::OMatrix_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__OMatrix_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__OVector2_Oval]"));
    core::BuiltInClass_sp classgeom__OVector2_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__OVector2_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__OVector2_Oval,_lisp,geom::OVector2_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::OVector2_O>::id,geom::OVector2_O::static_classSymbol());
    geom::OVector2_O::___staticClass = classgeom__OVector2_Oval;
#ifdef USE_MPS
    geom::OVector2_O::static_Kind = gctools::GCKind<geom::OVector2_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__OVector2_Oval,geom::OVector2_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::OVector2_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::OVector2_O>>::allocateClass();
        geom::OVector2_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::OVector2_O::static_className() % (void*)(geom::OVector2_O::static_allocator) );
    classgeom__OVector2_Oval->setCreator(geom::OVector2_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::OVector2_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__OVector2_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__OVector3_Oval]"));
    core::BuiltInClass_sp classgeom__OVector3_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__OVector3_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__OVector3_Oval,_lisp,geom::OVector3_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::OVector3_O>::id,geom::OVector3_O::static_classSymbol());
    geom::OVector3_O::___staticClass = classgeom__OVector3_Oval;
#ifdef USE_MPS
    geom::OVector3_O::static_Kind = gctools::GCKind<geom::OVector3_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__OVector3_Oval,geom::OVector3_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::OVector3_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::OVector3_O>>::allocateClass();
        geom::OVector3_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::OVector3_O::static_className() % (void*)(geom::OVector3_O::static_allocator) );
    classgeom__OVector3_Oval->setCreator(geom::OVector3_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::OVector3_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__OVector3_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__APFloat_Oval]"));
    core::BuiltInClass_sp classllvmo__APFloat_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__APFloat_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__APFloat_Oval,_lisp,llvmo::APFloat_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::APFloat_O>::id,llvmo::APFloat_O::static_classSymbol());
    llvmo::APFloat_O::___staticClass = classllvmo__APFloat_Oval;
#ifdef USE_MPS
    llvmo::APFloat_O::static_Kind = gctools::GCKind<llvmo::APFloat_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__APFloat_Oval,llvmo::APFloat_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::APFloat_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::APFloat_O>>::allocateClass();
        llvmo::APFloat_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::APFloat_O::static_className() % (void*)(llvmo::APFloat_O::static_allocator) );
    classllvmo__APFloat_Oval->setCreator(llvmo::APFloat_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::APFloat_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__APFloat_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__APInt_Oval]"));
    core::BuiltInClass_sp classllvmo__APInt_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__APInt_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__APInt_Oval,_lisp,llvmo::APInt_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::APInt_O>::id,llvmo::APInt_O::static_classSymbol());
    llvmo::APInt_O::___staticClass = classllvmo__APInt_Oval;
#ifdef USE_MPS
    llvmo::APInt_O::static_Kind = gctools::GCKind<llvmo::APInt_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__APInt_Oval,llvmo::APInt_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::APInt_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::APInt_O>>::allocateClass();
        llvmo::APInt_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::APInt_O::static_className() % (void*)(llvmo::APInt_O::static_allocator) );
    classllvmo__APInt_Oval->setCreator(llvmo::APInt_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::APInt_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__APInt_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIArray_Oval]"));
    core::BuiltInClass_sp classllvmo__DIArray_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIArray_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIArray_Oval,_lisp,llvmo::DIArray_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIArray_O>::id,llvmo::DIArray_O::static_classSymbol());
    llvmo::DIArray_O::___staticClass = classllvmo__DIArray_Oval;
#ifdef USE_MPS
    llvmo::DIArray_O::static_Kind = gctools::GCKind<llvmo::DIArray_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIArray_Oval,llvmo::DIArray_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIArray_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIArray_O>>::allocateClass();
        llvmo::DIArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIArray_O::static_className() % (void*)(llvmo::DIArray_O::static_allocator) );
    classllvmo__DIArray_Oval->setCreator(llvmo::DIArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIBasicType_Oval]"));
    core::BuiltInClass_sp classllvmo__DIBasicType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIBasicType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIBasicType_Oval,_lisp,llvmo::DIBasicType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIBasicType_O>::id,llvmo::DIBasicType_O::static_classSymbol());
    llvmo::DIBasicType_O::___staticClass = classllvmo__DIBasicType_Oval;
#ifdef USE_MPS
    llvmo::DIBasicType_O::static_Kind = gctools::GCKind<llvmo::DIBasicType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIBasicType_Oval,llvmo::DIBasicType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIBasicType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIBasicType_O>>::allocateClass();
        llvmo::DIBasicType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIBasicType_O::static_className() % (void*)(llvmo::DIBasicType_O::static_allocator) );
    classllvmo__DIBasicType_Oval->setCreator(llvmo::DIBasicType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIBasicType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIBasicType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIBuilder_Oval]"));
    core::BuiltInClass_sp classllvmo__DIBuilder_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIBuilder_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIBuilder_Oval,_lisp,llvmo::DIBuilder_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIBuilder_O>::id,llvmo::DIBuilder_O::static_classSymbol());
    llvmo::DIBuilder_O::___staticClass = classllvmo__DIBuilder_Oval;
#ifdef USE_MPS
    llvmo::DIBuilder_O::static_Kind = gctools::GCKind<llvmo::DIBuilder_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIBuilder_Oval,llvmo::DIBuilder_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIBuilder_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIBuilder_O>>::allocateClass();
        llvmo::DIBuilder_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIBuilder_O::static_className() % (void*)(llvmo::DIBuilder_O::static_allocator) );
    classllvmo__DIBuilder_Oval->setCreator(llvmo::DIBuilder_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIBuilder_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIBuilder_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DICompileUnit_Oval]"));
    core::BuiltInClass_sp classllvmo__DICompileUnit_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DICompileUnit_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DICompileUnit_Oval,_lisp,llvmo::DICompileUnit_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DICompileUnit_O>::id,llvmo::DICompileUnit_O::static_classSymbol());
    llvmo::DICompileUnit_O::___staticClass = classllvmo__DICompileUnit_Oval;
#ifdef USE_MPS
    llvmo::DICompileUnit_O::static_Kind = gctools::GCKind<llvmo::DICompileUnit_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DICompileUnit_Oval,llvmo::DICompileUnit_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DICompileUnit_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DICompileUnit_O>>::allocateClass();
        llvmo::DICompileUnit_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DICompileUnit_O::static_className() % (void*)(llvmo::DICompileUnit_O::static_allocator) );
    classllvmo__DICompileUnit_Oval->setCreator(llvmo::DICompileUnit_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DICompileUnit_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DICompileUnit_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DICompositeType_Oval]"));
    core::BuiltInClass_sp classllvmo__DICompositeType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DICompositeType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DICompositeType_Oval,_lisp,llvmo::DICompositeType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DICompositeType_O>::id,llvmo::DICompositeType_O::static_classSymbol());
    llvmo::DICompositeType_O::___staticClass = classllvmo__DICompositeType_Oval;
#ifdef USE_MPS
    llvmo::DICompositeType_O::static_Kind = gctools::GCKind<llvmo::DICompositeType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DICompositeType_Oval,llvmo::DICompositeType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DICompositeType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DICompositeType_O>>::allocateClass();
        llvmo::DICompositeType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DICompositeType_O::static_className() % (void*)(llvmo::DICompositeType_O::static_allocator) );
    classllvmo__DICompositeType_Oval->setCreator(llvmo::DICompositeType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DICompositeType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DICompositeType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIDerivedType_Oval]"));
    core::BuiltInClass_sp classllvmo__DIDerivedType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIDerivedType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIDerivedType_Oval,_lisp,llvmo::DIDerivedType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIDerivedType_O>::id,llvmo::DIDerivedType_O::static_classSymbol());
    llvmo::DIDerivedType_O::___staticClass = classllvmo__DIDerivedType_Oval;
#ifdef USE_MPS
    llvmo::DIDerivedType_O::static_Kind = gctools::GCKind<llvmo::DIDerivedType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIDerivedType_Oval,llvmo::DIDerivedType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIDerivedType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIDerivedType_O>>::allocateClass();
        llvmo::DIDerivedType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIDerivedType_O::static_className() % (void*)(llvmo::DIDerivedType_O::static_allocator) );
    classllvmo__DIDerivedType_Oval->setCreator(llvmo::DIDerivedType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIDerivedType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIDerivedType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIDescriptor_Oval]"));
    core::BuiltInClass_sp classllvmo__DIDescriptor_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIDescriptor_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIDescriptor_Oval,_lisp,llvmo::DIDescriptor_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIDescriptor_O>::id,llvmo::DIDescriptor_O::static_classSymbol());
    llvmo::DIDescriptor_O::___staticClass = classllvmo__DIDescriptor_Oval;
#ifdef USE_MPS
    llvmo::DIDescriptor_O::static_Kind = gctools::GCKind<llvmo::DIDescriptor_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIDescriptor_Oval,llvmo::DIDescriptor_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIDescriptor_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIDescriptor_O>>::allocateClass();
        llvmo::DIDescriptor_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIDescriptor_O::static_className() % (void*)(llvmo::DIDescriptor_O::static_allocator) );
    classllvmo__DIDescriptor_Oval->setCreator(llvmo::DIDescriptor_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIDescriptor_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIDescriptor_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIFile_Oval]"));
    core::BuiltInClass_sp classllvmo__DIFile_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIFile_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIFile_Oval,_lisp,llvmo::DIFile_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIFile_O>::id,llvmo::DIFile_O::static_classSymbol());
    llvmo::DIFile_O::___staticClass = classllvmo__DIFile_Oval;
#ifdef USE_MPS
    llvmo::DIFile_O::static_Kind = gctools::GCKind<llvmo::DIFile_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIFile_Oval,llvmo::DIFile_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIFile_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIFile_O>>::allocateClass();
        llvmo::DIFile_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIFile_O::static_className() % (void*)(llvmo::DIFile_O::static_allocator) );
    classllvmo__DIFile_Oval->setCreator(llvmo::DIFile_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIFile_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIFile_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DILexicalBlock_Oval]"));
    core::BuiltInClass_sp classllvmo__DILexicalBlock_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DILexicalBlock_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DILexicalBlock_Oval,_lisp,llvmo::DILexicalBlock_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DILexicalBlock_O>::id,llvmo::DILexicalBlock_O::static_classSymbol());
    llvmo::DILexicalBlock_O::___staticClass = classllvmo__DILexicalBlock_Oval;
#ifdef USE_MPS
    llvmo::DILexicalBlock_O::static_Kind = gctools::GCKind<llvmo::DILexicalBlock_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DILexicalBlock_Oval,llvmo::DILexicalBlock_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DILexicalBlock_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DILexicalBlock_O>>::allocateClass();
        llvmo::DILexicalBlock_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DILexicalBlock_O::static_className() % (void*)(llvmo::DILexicalBlock_O::static_allocator) );
    classllvmo__DILexicalBlock_Oval->setCreator(llvmo::DILexicalBlock_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DILexicalBlock_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DILexicalBlock_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIScope_Oval]"));
    core::BuiltInClass_sp classllvmo__DIScope_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIScope_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIScope_Oval,_lisp,llvmo::DIScope_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIScope_O>::id,llvmo::DIScope_O::static_classSymbol());
    llvmo::DIScope_O::___staticClass = classllvmo__DIScope_Oval;
#ifdef USE_MPS
    llvmo::DIScope_O::static_Kind = gctools::GCKind<llvmo::DIScope_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIScope_Oval,llvmo::DIScope_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIScope_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIScope_O>>::allocateClass();
        llvmo::DIScope_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIScope_O::static_className() % (void*)(llvmo::DIScope_O::static_allocator) );
    classllvmo__DIScope_Oval->setCreator(llvmo::DIScope_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIScope_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIScope_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DISubprogram_Oval]"));
    core::BuiltInClass_sp classllvmo__DISubprogram_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DISubprogram_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DISubprogram_Oval,_lisp,llvmo::DISubprogram_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DISubprogram_O>::id,llvmo::DISubprogram_O::static_classSymbol());
    llvmo::DISubprogram_O::___staticClass = classllvmo__DISubprogram_Oval;
#ifdef USE_MPS
    llvmo::DISubprogram_O::static_Kind = gctools::GCKind<llvmo::DISubprogram_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DISubprogram_Oval,llvmo::DISubprogram_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DISubprogram_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DISubprogram_O>>::allocateClass();
        llvmo::DISubprogram_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DISubprogram_O::static_className() % (void*)(llvmo::DISubprogram_O::static_allocator) );
    classllvmo__DISubprogram_Oval->setCreator(llvmo::DISubprogram_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DISubprogram_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DISubprogram_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DISubroutineType_Oval]"));
    core::BuiltInClass_sp classllvmo__DISubroutineType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DISubroutineType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DISubroutineType_Oval,_lisp,llvmo::DISubroutineType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DISubroutineType_O>::id,llvmo::DISubroutineType_O::static_classSymbol());
    llvmo::DISubroutineType_O::___staticClass = classllvmo__DISubroutineType_Oval;
#ifdef USE_MPS
    llvmo::DISubroutineType_O::static_Kind = gctools::GCKind<llvmo::DISubroutineType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DISubroutineType_Oval,llvmo::DISubroutineType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DISubroutineType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DISubroutineType_O>>::allocateClass();
        llvmo::DISubroutineType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DISubroutineType_O::static_className() % (void*)(llvmo::DISubroutineType_O::static_allocator) );
    classllvmo__DISubroutineType_Oval->setCreator(llvmo::DISubroutineType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DISubroutineType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DISubroutineType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DITypeArray_Oval]"));
    core::BuiltInClass_sp classllvmo__DITypeArray_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DITypeArray_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DITypeArray_Oval,_lisp,llvmo::DITypeArray_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DITypeArray_O>::id,llvmo::DITypeArray_O::static_classSymbol());
    llvmo::DITypeArray_O::___staticClass = classllvmo__DITypeArray_Oval;
#ifdef USE_MPS
    llvmo::DITypeArray_O::static_Kind = gctools::GCKind<llvmo::DITypeArray_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DITypeArray_Oval,llvmo::DITypeArray_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DITypeArray_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DITypeArray_O>>::allocateClass();
        llvmo::DITypeArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DITypeArray_O::static_className() % (void*)(llvmo::DITypeArray_O::static_allocator) );
    classllvmo__DITypeArray_Oval->setCreator(llvmo::DITypeArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DITypeArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DITypeArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DIType_Oval]"));
    core::BuiltInClass_sp classllvmo__DIType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DIType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DIType_Oval,_lisp,llvmo::DIType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DIType_O>::id,llvmo::DIType_O::static_classSymbol());
    llvmo::DIType_O::___staticClass = classllvmo__DIType_Oval;
#ifdef USE_MPS
    llvmo::DIType_O::static_Kind = gctools::GCKind<llvmo::DIType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DIType_Oval,llvmo::DIType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DIType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DIType_O>>::allocateClass();
        llvmo::DIType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DIType_O::static_className() % (void*)(llvmo::DIType_O::static_allocator) );
    classllvmo__DIType_Oval->setCreator(llvmo::DIType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DIType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DIType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DataLayout_Oval]"));
    core::BuiltInClass_sp classllvmo__DataLayout_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DataLayout_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DataLayout_Oval,_lisp,llvmo::DataLayout_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DataLayout_O>::id,llvmo::DataLayout_O::static_classSymbol());
    llvmo::DataLayout_O::___staticClass = classllvmo__DataLayout_Oval;
#ifdef USE_MPS
    llvmo::DataLayout_O::static_Kind = gctools::GCKind<llvmo::DataLayout_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DataLayout_Oval,llvmo::DataLayout_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DataLayout_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DataLayout_O>>::allocateClass();
        llvmo::DataLayout_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DataLayout_O::static_className() % (void*)(llvmo::DataLayout_O::static_allocator) );
    classllvmo__DataLayout_Oval->setCreator(llvmo::DataLayout_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DataLayout_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DataLayout_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__EngineBuilder_Oval]"));
    core::BuiltInClass_sp classllvmo__EngineBuilder_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__EngineBuilder_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__EngineBuilder_Oval,_lisp,llvmo::EngineBuilder_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::EngineBuilder_O>::id,llvmo::EngineBuilder_O::static_classSymbol());
    llvmo::EngineBuilder_O::___staticClass = classllvmo__EngineBuilder_Oval;
#ifdef USE_MPS
    llvmo::EngineBuilder_O::static_Kind = gctools::GCKind<llvmo::EngineBuilder_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__EngineBuilder_Oval,llvmo::EngineBuilder_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::EngineBuilder_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::EngineBuilder_O>>::allocateClass();
        llvmo::EngineBuilder_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::EngineBuilder_O::static_className() % (void*)(llvmo::EngineBuilder_O::static_allocator) );
    classllvmo__EngineBuilder_Oval->setCreator(llvmo::EngineBuilder_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::EngineBuilder_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__EngineBuilder_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ExecutionEngine_Oval]"));
    core::BuiltInClass_sp classllvmo__ExecutionEngine_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ExecutionEngine_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ExecutionEngine_Oval,_lisp,llvmo::ExecutionEngine_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ExecutionEngine_O>::id,llvmo::ExecutionEngine_O::static_classSymbol());
    llvmo::ExecutionEngine_O::___staticClass = classllvmo__ExecutionEngine_Oval;
#ifdef USE_MPS
    llvmo::ExecutionEngine_O::static_Kind = gctools::GCKind<llvmo::ExecutionEngine_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ExecutionEngine_Oval,llvmo::ExecutionEngine_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ExecutionEngine_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ExecutionEngine_O>>::allocateClass();
        llvmo::ExecutionEngine_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ExecutionEngine_O::static_className() % (void*)(llvmo::ExecutionEngine_O::static_allocator) );
    classllvmo__ExecutionEngine_Oval->setCreator(llvmo::ExecutionEngine_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ExecutionEngine_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ExecutionEngine_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__IRBuilderBase_Oval]"));
    core::BuiltInClass_sp classllvmo__IRBuilderBase_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__IRBuilderBase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__IRBuilderBase_Oval,_lisp,llvmo::IRBuilderBase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::IRBuilderBase_O>::id,llvmo::IRBuilderBase_O::static_classSymbol());
    llvmo::IRBuilderBase_O::___staticClass = classllvmo__IRBuilderBase_Oval;
#ifdef USE_MPS
    llvmo::IRBuilderBase_O::static_Kind = gctools::GCKind<llvmo::IRBuilderBase_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__IRBuilderBase_Oval,llvmo::IRBuilderBase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::IRBuilderBase_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::IRBuilderBase_O>>::allocateClass();
        llvmo::IRBuilderBase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::IRBuilderBase_O::static_className() % (void*)(llvmo::IRBuilderBase_O::static_allocator) );
    classllvmo__IRBuilderBase_Oval->setCreator(llvmo::IRBuilderBase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::IRBuilderBase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__IRBuilderBase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__LLVMContext_Oval]"));
    core::BuiltInClass_sp classllvmo__LLVMContext_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__LLVMContext_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__LLVMContext_Oval,_lisp,llvmo::LLVMContext_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::LLVMContext_O>::id,llvmo::LLVMContext_O::static_classSymbol());
    llvmo::LLVMContext_O::___staticClass = classllvmo__LLVMContext_Oval;
#ifdef USE_MPS
    llvmo::LLVMContext_O::static_Kind = gctools::GCKind<llvmo::LLVMContext_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__LLVMContext_Oval,llvmo::LLVMContext_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::LLVMContext_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::LLVMContext_O>>::allocateClass();
        llvmo::LLVMContext_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::LLVMContext_O::static_className() % (void*)(llvmo::LLVMContext_O::static_allocator) );
    classllvmo__LLVMContext_Oval->setCreator(llvmo::LLVMContext_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::LLVMContext_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__LLVMContext_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Linker_Oval]"));
    core::BuiltInClass_sp classllvmo__Linker_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Linker_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Linker_Oval,_lisp,llvmo::Linker_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Linker_O>::id,llvmo::Linker_O::static_classSymbol());
    llvmo::Linker_O::___staticClass = classllvmo__Linker_Oval;
#ifdef USE_MPS
    llvmo::Linker_O::static_Kind = gctools::GCKind<llvmo::Linker_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Linker_Oval,llvmo::Linker_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Linker_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Linker_O>>::allocateClass();
        llvmo::Linker_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Linker_O::static_className() % (void*)(llvmo::Linker_O::static_allocator) );
    classllvmo__Linker_Oval->setCreator(llvmo::Linker_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Linker_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Linker_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__MCSubtargetInfo_Oval]"));
    core::BuiltInClass_sp classllvmo__MCSubtargetInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__MCSubtargetInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__MCSubtargetInfo_Oval,_lisp,llvmo::MCSubtargetInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::MCSubtargetInfo_O>::id,llvmo::MCSubtargetInfo_O::static_classSymbol());
    llvmo::MCSubtargetInfo_O::___staticClass = classllvmo__MCSubtargetInfo_Oval;
#ifdef USE_MPS
    llvmo::MCSubtargetInfo_O::static_Kind = gctools::GCKind<llvmo::MCSubtargetInfo_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__MCSubtargetInfo_Oval,llvmo::MCSubtargetInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::MCSubtargetInfo_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::MCSubtargetInfo_O>>::allocateClass();
        llvmo::MCSubtargetInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::MCSubtargetInfo_O::static_className() % (void*)(llvmo::MCSubtargetInfo_O::static_allocator) );
    classllvmo__MCSubtargetInfo_Oval->setCreator(llvmo::MCSubtargetInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::MCSubtargetInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__MCSubtargetInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Metadata_Oval]"));
    core::BuiltInClass_sp classllvmo__Metadata_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Metadata_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Metadata_Oval,_lisp,llvmo::Metadata_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Metadata_O>::id,llvmo::Metadata_O::static_classSymbol());
    llvmo::Metadata_O::___staticClass = classllvmo__Metadata_Oval;
#ifdef USE_MPS
    llvmo::Metadata_O::static_Kind = gctools::GCKind<llvmo::Metadata_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Metadata_Oval,llvmo::Metadata_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Metadata_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Metadata_O>>::allocateClass();
        llvmo::Metadata_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Metadata_O::static_className() % (void*)(llvmo::Metadata_O::static_allocator) );
    classllvmo__Metadata_Oval->setCreator(llvmo::Metadata_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Metadata_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Metadata_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Module_Oval]"));
    core::BuiltInClass_sp classllvmo__Module_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Module_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Module_Oval,_lisp,llvmo::Module_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Module_O>::id,llvmo::Module_O::static_classSymbol());
    llvmo::Module_O::___staticClass = classllvmo__Module_Oval;
#ifdef USE_MPS
    llvmo::Module_O::static_Kind = gctools::GCKind<llvmo::Module_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Module_Oval,llvmo::Module_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Module_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Module_O>>::allocateClass();
        llvmo::Module_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Module_O::static_className() % (void*)(llvmo::Module_O::static_allocator) );
    classllvmo__Module_Oval->setCreator(llvmo::Module_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Module_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Module_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__NamedMDNode_Oval]"));
    core::BuiltInClass_sp classllvmo__NamedMDNode_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__NamedMDNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__NamedMDNode_Oval,_lisp,llvmo::NamedMDNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::NamedMDNode_O>::id,llvmo::NamedMDNode_O::static_classSymbol());
    llvmo::NamedMDNode_O::___staticClass = classllvmo__NamedMDNode_Oval;
#ifdef USE_MPS
    llvmo::NamedMDNode_O::static_Kind = gctools::GCKind<llvmo::NamedMDNode_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__NamedMDNode_Oval,llvmo::NamedMDNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::NamedMDNode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::NamedMDNode_O>>::allocateClass();
        llvmo::NamedMDNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::NamedMDNode_O::static_className() % (void*)(llvmo::NamedMDNode_O::static_allocator) );
    classllvmo__NamedMDNode_Oval->setCreator(llvmo::NamedMDNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::NamedMDNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__NamedMDNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__PassManagerBase_Oval]"));
    core::BuiltInClass_sp classllvmo__PassManagerBase_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__PassManagerBase_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__PassManagerBase_Oval,_lisp,llvmo::PassManagerBase_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::PassManagerBase_O>::id,llvmo::PassManagerBase_O::static_classSymbol());
    llvmo::PassManagerBase_O::___staticClass = classllvmo__PassManagerBase_Oval;
#ifdef USE_MPS
    llvmo::PassManagerBase_O::static_Kind = gctools::GCKind<llvmo::PassManagerBase_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__PassManagerBase_Oval,llvmo::PassManagerBase_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::PassManagerBase_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::PassManagerBase_O>>::allocateClass();
        llvmo::PassManagerBase_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::PassManagerBase_O::static_className() % (void*)(llvmo::PassManagerBase_O::static_allocator) );
    classllvmo__PassManagerBase_Oval->setCreator(llvmo::PassManagerBase_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::PassManagerBase_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__PassManagerBase_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__PassManagerBuilder_Oval]"));
    core::BuiltInClass_sp classllvmo__PassManagerBuilder_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__PassManagerBuilder_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__PassManagerBuilder_Oval,_lisp,llvmo::PassManagerBuilder_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::PassManagerBuilder_O>::id,llvmo::PassManagerBuilder_O::static_classSymbol());
    llvmo::PassManagerBuilder_O::___staticClass = classllvmo__PassManagerBuilder_Oval;
#ifdef USE_MPS
    llvmo::PassManagerBuilder_O::static_Kind = gctools::GCKind<llvmo::PassManagerBuilder_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__PassManagerBuilder_Oval,llvmo::PassManagerBuilder_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::PassManagerBuilder_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::PassManagerBuilder_O>>::allocateClass();
        llvmo::PassManagerBuilder_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::PassManagerBuilder_O::static_className() % (void*)(llvmo::PassManagerBuilder_O::static_allocator) );
    classllvmo__PassManagerBuilder_Oval->setCreator(llvmo::PassManagerBuilder_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::PassManagerBuilder_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__PassManagerBuilder_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Pass_Oval]"));
    core::BuiltInClass_sp classllvmo__Pass_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Pass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Pass_Oval,_lisp,llvmo::Pass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Pass_O>::id,llvmo::Pass_O::static_classSymbol());
    llvmo::Pass_O::___staticClass = classllvmo__Pass_Oval;
#ifdef USE_MPS
    llvmo::Pass_O::static_Kind = gctools::GCKind<llvmo::Pass_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Pass_Oval,llvmo::Pass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Pass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Pass_O>>::allocateClass();
        llvmo::Pass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Pass_O::static_className() % (void*)(llvmo::Pass_O::static_allocator) );
    classllvmo__Pass_Oval->setCreator(llvmo::Pass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Pass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Pass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__TargetMachine_Oval]"));
    core::BuiltInClass_sp classllvmo__TargetMachine_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__TargetMachine_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__TargetMachine_Oval,_lisp,llvmo::TargetMachine_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::TargetMachine_O>::id,llvmo::TargetMachine_O::static_classSymbol());
    llvmo::TargetMachine_O::___staticClass = classllvmo__TargetMachine_Oval;
#ifdef USE_MPS
    llvmo::TargetMachine_O::static_Kind = gctools::GCKind<llvmo::TargetMachine_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__TargetMachine_Oval,llvmo::TargetMachine_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::TargetMachine_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::TargetMachine_O>>::allocateClass();
        llvmo::TargetMachine_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::TargetMachine_O::static_className() % (void*)(llvmo::TargetMachine_O::static_allocator) );
    classllvmo__TargetMachine_Oval->setCreator(llvmo::TargetMachine_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::TargetMachine_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__TargetMachine_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__TargetOptions_Oval]"));
    core::BuiltInClass_sp classllvmo__TargetOptions_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__TargetOptions_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__TargetOptions_Oval,_lisp,llvmo::TargetOptions_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::TargetOptions_O>::id,llvmo::TargetOptions_O::static_classSymbol());
    llvmo::TargetOptions_O::___staticClass = classllvmo__TargetOptions_Oval;
#ifdef USE_MPS
    llvmo::TargetOptions_O::static_Kind = gctools::GCKind<llvmo::TargetOptions_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__TargetOptions_Oval,llvmo::TargetOptions_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::TargetOptions_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::TargetOptions_O>>::allocateClass();
        llvmo::TargetOptions_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::TargetOptions_O::static_className() % (void*)(llvmo::TargetOptions_O::static_allocator) );
    classllvmo__TargetOptions_Oval->setCreator(llvmo::TargetOptions_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::TargetOptions_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__TargetOptions_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Target_Oval]"));
    core::BuiltInClass_sp classllvmo__Target_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Target_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Target_Oval,_lisp,llvmo::Target_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Target_O>::id,llvmo::Target_O::static_classSymbol());
    llvmo::Target_O::___staticClass = classllvmo__Target_Oval;
#ifdef USE_MPS
    llvmo::Target_O::static_Kind = gctools::GCKind<llvmo::Target_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Target_Oval,llvmo::Target_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Target_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Target_O>>::allocateClass();
        llvmo::Target_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Target_O::static_className() % (void*)(llvmo::Target_O::static_allocator) );
    classllvmo__Target_Oval->setCreator(llvmo::Target_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Target_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Target_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Triple_Oval]"));
    core::BuiltInClass_sp classllvmo__Triple_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Triple_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Triple_Oval,_lisp,llvmo::Triple_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Triple_O>::id,llvmo::Triple_O::static_classSymbol());
    llvmo::Triple_O::___staticClass = classllvmo__Triple_Oval;
#ifdef USE_MPS
    llvmo::Triple_O::static_Kind = gctools::GCKind<llvmo::Triple_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Triple_Oval,llvmo::Triple_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Triple_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Triple_O>>::allocateClass();
        llvmo::Triple_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Triple_O::static_className() % (void*)(llvmo::Triple_O::static_allocator) );
    classllvmo__Triple_Oval->setCreator(llvmo::Triple_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Triple_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Triple_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Type_Oval]"));
    core::BuiltInClass_sp classllvmo__Type_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Type_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Type_Oval,_lisp,llvmo::Type_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Type_O>::id,llvmo::Type_O::static_classSymbol());
    llvmo::Type_O::___staticClass = classllvmo__Type_Oval;
#ifdef USE_MPS
    llvmo::Type_O::static_Kind = gctools::GCKind<llvmo::Type_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Type_Oval,llvmo::Type_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Type_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Type_O>>::allocateClass();
        llvmo::Type_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Type_O::static_className() % (void*)(llvmo::Type_O::static_allocator) );
    classllvmo__Type_Oval->setCreator(llvmo::Type_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Type_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Type_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Value_Oval]"));
    core::BuiltInClass_sp classllvmo__Value_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Value_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Value_Oval,_lisp,llvmo::Value_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Value_O>::id,llvmo::Value_O::static_classSymbol());
    llvmo::Value_O::___staticClass = classllvmo__Value_Oval;
#ifdef USE_MPS
    llvmo::Value_O::static_Kind = gctools::GCKind<llvmo::Value_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Value_Oval,llvmo::Value_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Value_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Value_O>>::allocateClass();
        llvmo::Value_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Value_O::static_className() % (void*)(llvmo::Value_O::static_allocator) );
    classllvmo__Value_Oval->setCreator(llvmo::Value_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Value_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Value_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_UnitsPkg
    LOG(BF("Creating class[classunits__Dimension_Oval]"));
    core::BuiltInClass_sp classunits__Dimension_Oval = core::BuiltInClass_O::createUncollectable();
    classunits__Dimension_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__Dimension_Oval,_lisp,units::Dimension_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<units::Dimension_O>::id,units::Dimension_O::static_classSymbol());
    units::Dimension_O::___staticClass = classunits__Dimension_Oval;
#ifdef USE_MPS
    units::Dimension_O::static_Kind = gctools::GCKind<units::Dimension_O>::Kind;
#endif
    core::af_setf_findClass(classunits__Dimension_Oval,units::Dimension_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<units::Dimension_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<units::Dimension_O>>::allocateClass();
        units::Dimension_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::Dimension_O::static_className() % (void*)(units::Dimension_O::static_allocator) );
    classunits__Dimension_Oval->setCreator(units::Dimension_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::Dimension_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__Dimension_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_UnitsPkg
    LOG(BF("Creating class[classunits__Quantity_Oval]"));
    core::BuiltInClass_sp classunits__Quantity_Oval = core::BuiltInClass_O::createUncollectable();
    classunits__Quantity_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__Quantity_Oval,_lisp,units::Quantity_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<units::Quantity_O>::id,units::Quantity_O::static_classSymbol());
    units::Quantity_O::___staticClass = classunits__Quantity_Oval;
#ifdef USE_MPS
    units::Quantity_O::static_Kind = gctools::GCKind<units::Quantity_O>::Kind;
#endif
    core::af_setf_findClass(classunits__Quantity_Oval,units::Quantity_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<units::Quantity_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<units::Quantity_O>>::allocateClass();
        units::Quantity_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::Quantity_O::static_className() % (void*)(units::Quantity_O::static_allocator) );
    classunits__Quantity_Oval->setCreator(units::Quantity_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::Quantity_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__Quantity_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_UnitsPkg
    LOG(BF("Creating class[classunits__Unit_Oval]"));
    core::BuiltInClass_sp classunits__Unit_Oval = core::BuiltInClass_O::createUncollectable();
    classunits__Unit_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__Unit_Oval,_lisp,units::Unit_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<units::Unit_O>::id,units::Unit_O::static_classSymbol());
    units::Unit_O::___staticClass = classunits__Unit_Oval;
#ifdef USE_MPS
    units::Unit_O::static_Kind = gctools::GCKind<units::Unit_O>::Kind;
#endif
    core::af_setf_findClass(classunits__Unit_Oval,units::Unit_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<units::Unit_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<units::Unit_O>>::allocateClass();
        units::Unit_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::Unit_O::static_className() % (void*)(units::Unit_O::static_allocator) );
    classunits__Unit_Oval->setCreator(units::Unit_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::Unit_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__Unit_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__BitVector_Oval]"));
    core::BuiltInClass_sp classcore__BitVector_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__BitVector_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__BitVector_Oval,_lisp,core::BitVector_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::BitVector_O>::id,core::BitVector_O::static_classSymbol());
    core::BitVector_O::___staticClass = classcore__BitVector_Oval;
#ifdef USE_MPS
    core::BitVector_O::static_Kind = gctools::GCKind<core::BitVector_O>::Kind;
#endif
    core::af_setf_findClass(classcore__BitVector_Oval,core::BitVector_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::BitVector_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::BitVector_O>>::allocateClass();
        core::BitVector_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::BitVector_O::static_className() % (void*)(core::BitVector_O::static_allocator) );
    classcore__BitVector_Oval->setCreator(core::BitVector_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::BitVector_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__BitVector_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__BroadcastStream_Oval]"));
    core::BuiltInClass_sp classcore__BroadcastStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__BroadcastStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__BroadcastStream_Oval,_lisp,core::BroadcastStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::BroadcastStream_O>::id,core::BroadcastStream_O::static_classSymbol());
    core::BroadcastStream_O::___staticClass = classcore__BroadcastStream_Oval;
#ifdef USE_MPS
    core::BroadcastStream_O::static_Kind = gctools::GCKind<core::BroadcastStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__BroadcastStream_Oval,core::BroadcastStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::BroadcastStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::BroadcastStream_O>>::allocateClass();
        core::BroadcastStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::BroadcastStream_O::static_className() % (void*)(core::BroadcastStream_O::static_allocator) );
    classcore__BroadcastStream_Oval->setCreator(core::BroadcastStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::BroadcastStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__BroadcastStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__CompileTimeEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__CompileTimeEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__CompileTimeEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__CompileTimeEnvironment_Oval,_lisp,core::CompileTimeEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::CompileTimeEnvironment_O>::id,core::CompileTimeEnvironment_O::static_classSymbol());
    core::CompileTimeEnvironment_O::___staticClass = classcore__CompileTimeEnvironment_Oval;
#ifdef USE_MPS
    core::CompileTimeEnvironment_O::static_Kind = gctools::GCKind<core::CompileTimeEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__CompileTimeEnvironment_Oval,core::CompileTimeEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::CompileTimeEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::CompileTimeEnvironment_O>>::allocateClass();
        core::CompileTimeEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::CompileTimeEnvironment_O::static_className() % (void*)(core::CompileTimeEnvironment_O::static_allocator) );
    classcore__CompileTimeEnvironment_Oval->setCreator(core::CompileTimeEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::CompileTimeEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__CompileTimeEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ConcatenatedStream_Oval]"));
    core::BuiltInClass_sp classcore__ConcatenatedStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ConcatenatedStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ConcatenatedStream_Oval,_lisp,core::ConcatenatedStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ConcatenatedStream_O>::id,core::ConcatenatedStream_O::static_classSymbol());
    core::ConcatenatedStream_O::___staticClass = classcore__ConcatenatedStream_Oval;
#ifdef USE_MPS
    core::ConcatenatedStream_O::static_Kind = gctools::GCKind<core::ConcatenatedStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ConcatenatedStream_Oval,core::ConcatenatedStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ConcatenatedStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ConcatenatedStream_O>>::allocateClass();
        core::ConcatenatedStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ConcatenatedStream_O::static_className() % (void*)(core::ConcatenatedStream_O::static_allocator) );
    classcore__ConcatenatedStream_Oval->setCreator(core::ConcatenatedStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ConcatenatedStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ConcatenatedStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__EchoStream_Oval]"));
    core::BuiltInClass_sp classcore__EchoStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__EchoStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__EchoStream_Oval,_lisp,core::EchoStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::EchoStream_O>::id,core::EchoStream_O::static_classSymbol());
    core::EchoStream_O::___staticClass = classcore__EchoStream_Oval;
#ifdef USE_MPS
    core::EchoStream_O::static_Kind = gctools::GCKind<core::EchoStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__EchoStream_Oval,core::EchoStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::EchoStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::EchoStream_O>>::allocateClass();
        core::EchoStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::EchoStream_O::static_className() % (void*)(core::EchoStream_O::static_allocator) );
    classcore__EchoStream_Oval->setCreator(core::EchoStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::EchoStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__EchoStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__FileStream_Oval]"));
    core::BuiltInClass_sp classcore__FileStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__FileStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__FileStream_Oval,_lisp,core::FileStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::FileStream_O>::id,core::FileStream_O::static_classSymbol());
    core::FileStream_O::___staticClass = classcore__FileStream_Oval;
#ifdef USE_MPS
    core::FileStream_O::static_Kind = gctools::GCKind<core::FileStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__FileStream_Oval,core::FileStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::FileStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::FileStream_O>>::allocateClass();
        core::FileStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::FileStream_O::static_className() % (void*)(core::FileStream_O::static_allocator) );
    classcore__FileStream_Oval->setCreator(core::FileStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::FileStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__FileStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Float_Oval]"));
    core::BuiltInClass_sp classcore__Float_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Float_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Float_Oval,_lisp,core::Float_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Float_O>::id,core::Float_O::static_classSymbol());
    core::Float_O::___staticClass = classcore__Float_Oval;
#ifdef USE_MPS
    core::Float_O::static_Kind = gctools::GCKind<core::Float_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Float_Oval,core::Float_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Float_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Float_O>>::allocateClass();
        core::Float_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Float_O::static_className() % (void*)(core::Float_O::static_allocator) );
    classcore__Float_Oval->setCreator(core::Float_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Float_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Float_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__FunctionFrame_Oval]"));
    core::BuiltInClass_sp classcore__FunctionFrame_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__FunctionFrame_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__FunctionFrame_Oval,_lisp,core::FunctionFrame_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::FunctionFrame_O>::id,core::FunctionFrame_O::static_classSymbol());
    core::FunctionFrame_O::___staticClass = classcore__FunctionFrame_Oval;
#ifdef USE_MPS
    core::FunctionFrame_O::static_Kind = gctools::GCKind<core::FunctionFrame_O>::Kind;
#endif
    core::af_setf_findClass(classcore__FunctionFrame_Oval,core::FunctionFrame_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::FunctionFrame_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::FunctionFrame_O>>::allocateClass();
        core::FunctionFrame_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::FunctionFrame_O::static_className() % (void*)(core::FunctionFrame_O::static_allocator) );
    classcore__FunctionFrame_Oval->setCreator(core::FunctionFrame_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::FunctionFrame_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__FunctionFrame_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Rational_Oval]"));
    core::BuiltInClass_sp classcore__Rational_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Rational_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Rational_Oval,_lisp,core::Rational_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Rational_O>::id,core::Rational_O::static_classSymbol());
    core::Rational_O::___staticClass = classcore__Rational_Oval;
#ifdef USE_MPS
    core::Rational_O::static_Kind = gctools::GCKind<core::Rational_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Rational_Oval,core::Rational_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Rational_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Rational_O>>::allocateClass();
        core::Rational_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Rational_O::static_className() % (void*)(core::Rational_O::static_allocator) );
    classcore__Rational_Oval->setCreator(core::Rational_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Rational_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Rational_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__RuntimeVisibleEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__RuntimeVisibleEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__RuntimeVisibleEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__RuntimeVisibleEnvironment_Oval,_lisp,core::RuntimeVisibleEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::RuntimeVisibleEnvironment_O>::id,core::RuntimeVisibleEnvironment_O::static_classSymbol());
    core::RuntimeVisibleEnvironment_O::___staticClass = classcore__RuntimeVisibleEnvironment_Oval;
#ifdef USE_MPS
    core::RuntimeVisibleEnvironment_O::static_Kind = gctools::GCKind<core::RuntimeVisibleEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__RuntimeVisibleEnvironment_Oval,core::RuntimeVisibleEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::RuntimeVisibleEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::RuntimeVisibleEnvironment_O>>::allocateClass();
        core::RuntimeVisibleEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::RuntimeVisibleEnvironment_O::static_className() % (void*)(core::RuntimeVisibleEnvironment_O::static_allocator) );
    classcore__RuntimeVisibleEnvironment_Oval->setCreator(core::RuntimeVisibleEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::RuntimeVisibleEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__RuntimeVisibleEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SexpLoadArchive_Oval]"));
    core::BuiltInClass_sp classcore__SexpLoadArchive_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SexpLoadArchive_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SexpLoadArchive_Oval,_lisp,core::SexpLoadArchive_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SexpLoadArchive_O>::id,core::SexpLoadArchive_O::static_classSymbol());
    core::SexpLoadArchive_O::___staticClass = classcore__SexpLoadArchive_Oval;
#ifdef USE_MPS
    core::SexpLoadArchive_O::static_Kind = gctools::GCKind<core::SexpLoadArchive_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SexpLoadArchive_Oval,core::SexpLoadArchive_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SexpLoadArchive_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SexpLoadArchive_O>>::allocateClass();
        core::SexpLoadArchive_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SexpLoadArchive_O::static_className() % (void*)(core::SexpLoadArchive_O::static_allocator) );
    classcore__SexpLoadArchive_Oval->setCreator(core::SexpLoadArchive_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SexpLoadArchive_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SexpLoadArchive_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SexpSaveArchive_Oval]"));
    core::BuiltInClass_sp classcore__SexpSaveArchive_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SexpSaveArchive_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SexpSaveArchive_Oval,_lisp,core::SexpSaveArchive_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SexpSaveArchive_O>::id,core::SexpSaveArchive_O::static_classSymbol());
    core::SexpSaveArchive_O::___staticClass = classcore__SexpSaveArchive_Oval;
#ifdef USE_MPS
    core::SexpSaveArchive_O::static_Kind = gctools::GCKind<core::SexpSaveArchive_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SexpSaveArchive_Oval,core::SexpSaveArchive_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SexpSaveArchive_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SexpSaveArchive_O>>::allocateClass();
        core::SexpSaveArchive_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SexpSaveArchive_O::static_className() % (void*)(core::SexpSaveArchive_O::static_allocator) );
    classcore__SexpSaveArchive_Oval->setCreator(core::SexpSaveArchive_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SexpSaveArchive_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SexpSaveArchive_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Specializer_Oval]"));
    StandardClass_sp classcore__Specializer_Oval = StandardClass_O::createUncollectable();
    classcore__Specializer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Specializer_Oval,_lisp,core::Specializer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Specializer_O>::id,core::Specializer_O::static_classSymbol());
    core::Specializer_O::___staticClass = classcore__Specializer_Oval;
#ifdef USE_MPS
    core::Specializer_O::static_Kind = gctools::GCKind<core::Specializer_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Specializer_Oval,core::Specializer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Specializer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Specializer_O>>::allocateClass();
        core::Specializer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Specializer_O::static_className() % (void*)(core::Specializer_O::static_allocator) );
    classcore__Specializer_Oval->setCreator(core::Specializer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Specializer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Specializer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StringStream_Oval]"));
    core::BuiltInClass_sp classcore__StringStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__StringStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StringStream_Oval,_lisp,core::StringStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StringStream_O>::id,core::StringStream_O::static_classSymbol());
    core::StringStream_O::___staticClass = classcore__StringStream_Oval;
#ifdef USE_MPS
    core::StringStream_O::static_Kind = gctools::GCKind<core::StringStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StringStream_Oval,core::StringStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StringStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StringStream_O>>::allocateClass();
        core::StringStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StringStream_O::static_className() % (void*)(core::StringStream_O::static_allocator) );
    classcore__StringStream_Oval->setCreator(core::StringStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StringStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StringStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__String_Oval]"));
    core::BuiltInClass_sp classcore__String_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__String_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__String_Oval,_lisp,core::String_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::String_O>::id,core::String_O::static_classSymbol());
    core::String_O::___staticClass = classcore__String_Oval;
#ifdef USE_MPS
    core::String_O::static_Kind = gctools::GCKind<core::String_O>::Kind;
#endif
    core::af_setf_findClass(classcore__String_Oval,core::String_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::String_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::String_O>>::allocateClass();
        core::String_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::String_O::static_className() % (void*)(core::String_O::static_allocator) );
    classcore__String_Oval->setCreator(core::String_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::String_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__String_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SynonymStream_Oval]"));
    core::BuiltInClass_sp classcore__SynonymStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SynonymStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SynonymStream_Oval,_lisp,core::SynonymStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SynonymStream_O>::id,core::SynonymStream_O::static_classSymbol());
    core::SynonymStream_O::___staticClass = classcore__SynonymStream_Oval;
#ifdef USE_MPS
    core::SynonymStream_O::static_Kind = gctools::GCKind<core::SynonymStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SynonymStream_Oval,core::SynonymStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SynonymStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SynonymStream_O>>::allocateClass();
        core::SynonymStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SynonymStream_O::static_className() % (void*)(core::SynonymStream_O::static_allocator) );
    classcore__SynonymStream_Oval->setCreator(core::SynonymStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SynonymStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SynonymStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__TagbodyFrame_Oval]"));
    core::BuiltInClass_sp classcore__TagbodyFrame_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__TagbodyFrame_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__TagbodyFrame_Oval,_lisp,core::TagbodyFrame_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::TagbodyFrame_O>::id,core::TagbodyFrame_O::static_classSymbol());
    core::TagbodyFrame_O::___staticClass = classcore__TagbodyFrame_Oval;
#ifdef USE_MPS
    core::TagbodyFrame_O::static_Kind = gctools::GCKind<core::TagbodyFrame_O>::Kind;
#endif
    core::af_setf_findClass(classcore__TagbodyFrame_Oval,core::TagbodyFrame_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::TagbodyFrame_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::TagbodyFrame_O>>::allocateClass();
        core::TagbodyFrame_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::TagbodyFrame_O::static_className() % (void*)(core::TagbodyFrame_O::static_allocator) );
    classcore__TagbodyFrame_Oval->setCreator(core::TagbodyFrame_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::TagbodyFrame_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__TagbodyFrame_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__TwoWayStream_Oval]"));
    core::BuiltInClass_sp classcore__TwoWayStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__TwoWayStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__TwoWayStream_Oval,_lisp,core::TwoWayStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::TwoWayStream_O>::id,core::TwoWayStream_O::static_classSymbol());
    core::TwoWayStream_O::___staticClass = classcore__TwoWayStream_Oval;
#ifdef USE_MPS
    core::TwoWayStream_O::static_Kind = gctools::GCKind<core::TwoWayStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__TwoWayStream_Oval,core::TwoWayStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::TwoWayStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::TwoWayStream_O>>::allocateClass();
        core::TwoWayStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::TwoWayStream_O::static_className() % (void*)(core::TwoWayStream_O::static_allocator) );
    classcore__TwoWayStream_Oval->setCreator(core::TwoWayStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::TwoWayStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__TwoWayStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ValueFrame_Oval]"));
    core::BuiltInClass_sp classcore__ValueFrame_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ValueFrame_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ValueFrame_Oval,_lisp,core::ValueFrame_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ValueFrame_O>::id,core::ValueFrame_O::static_classSymbol());
    core::ValueFrame_O::___staticClass = classcore__ValueFrame_Oval;
#ifdef USE_MPS
    core::ValueFrame_O::static_Kind = gctools::GCKind<core::ValueFrame_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ValueFrame_Oval,core::ValueFrame_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ValueFrame_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ValueFrame_O>>::allocateClass();
        core::ValueFrame_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ValueFrame_O::static_className() % (void*)(core::ValueFrame_O::static_allocator) );
    classcore__ValueFrame_Oval->setCreator(core::ValueFrame_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ValueFrame_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ValueFrame_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__VectorObjects_Oval]"));
    core::BuiltInClass_sp classcore__VectorObjects_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__VectorObjects_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__VectorObjects_Oval,_lisp,core::VectorObjects_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::VectorObjects_O>::id,core::VectorObjects_O::static_classSymbol());
    core::VectorObjects_O::___staticClass = classcore__VectorObjects_Oval;
#ifdef USE_MPS
    core::VectorObjects_O::static_Kind = gctools::GCKind<core::VectorObjects_O>::Kind;
#endif
    core::af_setf_findClass(classcore__VectorObjects_Oval,core::VectorObjects_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::VectorObjects_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::VectorObjects_O>>::allocateClass();
        core::VectorObjects_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::VectorObjects_O::static_className() % (void*)(core::VectorObjects_O::static_allocator) );
    classcore__VectorObjects_Oval->setCreator(core::VectorObjects_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::VectorObjects_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__VectorObjects_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__CoordinateArray_Oval]"));
    core::BuiltInClass_sp classgeom__CoordinateArray_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__CoordinateArray_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__CoordinateArray_Oval,_lisp,geom::CoordinateArray_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::CoordinateArray_O>::id,geom::CoordinateArray_O::static_classSymbol());
    geom::CoordinateArray_O::___staticClass = classgeom__CoordinateArray_Oval;
#ifdef USE_MPS
    geom::CoordinateArray_O::static_Kind = gctools::GCKind<geom::CoordinateArray_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__CoordinateArray_Oval,geom::CoordinateArray_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::CoordinateArray_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::CoordinateArray_O>>::allocateClass();
        geom::CoordinateArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::CoordinateArray_O::static_className() % (void*)(geom::CoordinateArray_O::static_allocator) );
    classgeom__CoordinateArray_Oval->setCreator(geom::CoordinateArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::CoordinateArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__CoordinateArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Argument_Oval]"));
    core::BuiltInClass_sp classllvmo__Argument_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Argument_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Argument_Oval,_lisp,llvmo::Argument_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Argument_O>::id,llvmo::Argument_O::static_classSymbol());
    llvmo::Argument_O::___staticClass = classllvmo__Argument_Oval;
#ifdef USE_MPS
    llvmo::Argument_O::static_Kind = gctools::GCKind<llvmo::Argument_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Argument_Oval,llvmo::Argument_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Argument_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Argument_O>>::allocateClass();
        llvmo::Argument_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Argument_O::static_className() % (void*)(llvmo::Argument_O::static_allocator) );
    classllvmo__Argument_Oval->setCreator(llvmo::Argument_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Argument_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Argument_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__BasicBlock_Oval]"));
    core::BuiltInClass_sp classllvmo__BasicBlock_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__BasicBlock_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__BasicBlock_Oval,_lisp,llvmo::BasicBlock_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::BasicBlock_O>::id,llvmo::BasicBlock_O::static_classSymbol());
    llvmo::BasicBlock_O::___staticClass = classllvmo__BasicBlock_Oval;
#ifdef USE_MPS
    llvmo::BasicBlock_O::static_Kind = gctools::GCKind<llvmo::BasicBlock_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__BasicBlock_Oval,llvmo::BasicBlock_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::BasicBlock_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::BasicBlock_O>>::allocateClass();
        llvmo::BasicBlock_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::BasicBlock_O::static_className() % (void*)(llvmo::BasicBlock_O::static_allocator) );
    classllvmo__BasicBlock_Oval->setCreator(llvmo::BasicBlock_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::BasicBlock_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__BasicBlock_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__CompositeType_Oval]"));
    core::BuiltInClass_sp classllvmo__CompositeType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__CompositeType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__CompositeType_Oval,_lisp,llvmo::CompositeType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::CompositeType_O>::id,llvmo::CompositeType_O::static_classSymbol());
    llvmo::CompositeType_O::___staticClass = classllvmo__CompositeType_Oval;
#ifdef USE_MPS
    llvmo::CompositeType_O::static_Kind = gctools::GCKind<llvmo::CompositeType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__CompositeType_Oval,llvmo::CompositeType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::CompositeType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::CompositeType_O>>::allocateClass();
        llvmo::CompositeType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::CompositeType_O::static_className() % (void*)(llvmo::CompositeType_O::static_allocator) );
    classllvmo__CompositeType_Oval->setCreator(llvmo::CompositeType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::CompositeType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__CompositeType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__FunctionPassManager_Oval]"));
    core::BuiltInClass_sp classllvmo__FunctionPassManager_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__FunctionPassManager_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__FunctionPassManager_Oval,_lisp,llvmo::FunctionPassManager_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::FunctionPassManager_O>::id,llvmo::FunctionPassManager_O::static_classSymbol());
    llvmo::FunctionPassManager_O::___staticClass = classllvmo__FunctionPassManager_Oval;
#ifdef USE_MPS
    llvmo::FunctionPassManager_O::static_Kind = gctools::GCKind<llvmo::FunctionPassManager_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__FunctionPassManager_Oval,llvmo::FunctionPassManager_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::FunctionPassManager_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::FunctionPassManager_O>>::allocateClass();
        llvmo::FunctionPassManager_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::FunctionPassManager_O::static_className() % (void*)(llvmo::FunctionPassManager_O::static_allocator) );
    classllvmo__FunctionPassManager_Oval->setCreator(llvmo::FunctionPassManager_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::FunctionPassManager_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__FunctionPassManager_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__FunctionPass_Oval]"));
    core::BuiltInClass_sp classllvmo__FunctionPass_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__FunctionPass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__FunctionPass_Oval,_lisp,llvmo::FunctionPass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::FunctionPass_O>::id,llvmo::FunctionPass_O::static_classSymbol());
    llvmo::FunctionPass_O::___staticClass = classllvmo__FunctionPass_Oval;
#ifdef USE_MPS
    llvmo::FunctionPass_O::static_Kind = gctools::GCKind<llvmo::FunctionPass_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__FunctionPass_Oval,llvmo::FunctionPass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::FunctionPass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::FunctionPass_O>>::allocateClass();
        llvmo::FunctionPass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::FunctionPass_O::static_className() % (void*)(llvmo::FunctionPass_O::static_allocator) );
    classllvmo__FunctionPass_Oval->setCreator(llvmo::FunctionPass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::FunctionPass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__FunctionPass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__FunctionType_Oval]"));
    core::BuiltInClass_sp classllvmo__FunctionType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__FunctionType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__FunctionType_Oval,_lisp,llvmo::FunctionType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::FunctionType_O>::id,llvmo::FunctionType_O::static_classSymbol());
    llvmo::FunctionType_O::___staticClass = classllvmo__FunctionType_Oval;
#ifdef USE_MPS
    llvmo::FunctionType_O::static_Kind = gctools::GCKind<llvmo::FunctionType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__FunctionType_Oval,llvmo::FunctionType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::FunctionType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::FunctionType_O>>::allocateClass();
        llvmo::FunctionType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::FunctionType_O::static_className() % (void*)(llvmo::FunctionType_O::static_allocator) );
    classllvmo__FunctionType_Oval->setCreator(llvmo::FunctionType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::FunctionType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__FunctionType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__IRBuilder_Oval]"));
    core::BuiltInClass_sp classllvmo__IRBuilder_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__IRBuilder_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__IRBuilder_Oval,_lisp,llvmo::IRBuilder_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::IRBuilder_O>::id,llvmo::IRBuilder_O::static_classSymbol());
    llvmo::IRBuilder_O::___staticClass = classllvmo__IRBuilder_Oval;
#ifdef USE_MPS
    llvmo::IRBuilder_O::static_Kind = gctools::GCKind<llvmo::IRBuilder_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__IRBuilder_Oval,llvmo::IRBuilder_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::IRBuilder_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::IRBuilder_O>>::allocateClass();
        llvmo::IRBuilder_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::IRBuilder_O::static_className() % (void*)(llvmo::IRBuilder_O::static_allocator) );
    classllvmo__IRBuilder_Oval->setCreator(llvmo::IRBuilder_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::IRBuilder_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__IRBuilder_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__IntegerType_Oval]"));
    core::BuiltInClass_sp classllvmo__IntegerType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__IntegerType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__IntegerType_Oval,_lisp,llvmo::IntegerType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::IntegerType_O>::id,llvmo::IntegerType_O::static_classSymbol());
    llvmo::IntegerType_O::___staticClass = classllvmo__IntegerType_Oval;
#ifdef USE_MPS
    llvmo::IntegerType_O::static_Kind = gctools::GCKind<llvmo::IntegerType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__IntegerType_Oval,llvmo::IntegerType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::IntegerType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::IntegerType_O>>::allocateClass();
        llvmo::IntegerType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::IntegerType_O::static_className() % (void*)(llvmo::IntegerType_O::static_allocator) );
    classllvmo__IntegerType_Oval->setCreator(llvmo::IntegerType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::IntegerType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__IntegerType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__LLVMTargetMachine_Oval]"));
    core::BuiltInClass_sp classllvmo__LLVMTargetMachine_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__LLVMTargetMachine_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__LLVMTargetMachine_Oval,_lisp,llvmo::LLVMTargetMachine_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::LLVMTargetMachine_O>::id,llvmo::LLVMTargetMachine_O::static_classSymbol());
    llvmo::LLVMTargetMachine_O::___staticClass = classllvmo__LLVMTargetMachine_Oval;
#ifdef USE_MPS
    llvmo::LLVMTargetMachine_O::static_Kind = gctools::GCKind<llvmo::LLVMTargetMachine_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__LLVMTargetMachine_Oval,llvmo::LLVMTargetMachine_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::LLVMTargetMachine_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::LLVMTargetMachine_O>>::allocateClass();
        llvmo::LLVMTargetMachine_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::LLVMTargetMachine_O::static_className() % (void*)(llvmo::LLVMTargetMachine_O::static_allocator) );
    classllvmo__LLVMTargetMachine_Oval->setCreator(llvmo::LLVMTargetMachine_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::LLVMTargetMachine_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__LLVMTargetMachine_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__MDNode_Oval]"));
    core::BuiltInClass_sp classllvmo__MDNode_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__MDNode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__MDNode_Oval,_lisp,llvmo::MDNode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::MDNode_O>::id,llvmo::MDNode_O::static_classSymbol());
    llvmo::MDNode_O::___staticClass = classllvmo__MDNode_Oval;
#ifdef USE_MPS
    llvmo::MDNode_O::static_Kind = gctools::GCKind<llvmo::MDNode_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__MDNode_Oval,llvmo::MDNode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::MDNode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::MDNode_O>>::allocateClass();
        llvmo::MDNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::MDNode_O::static_className() % (void*)(llvmo::MDNode_O::static_allocator) );
    classllvmo__MDNode_Oval->setCreator(llvmo::MDNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::MDNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__MDNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__MDString_Oval]"));
    core::BuiltInClass_sp classllvmo__MDString_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__MDString_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__MDString_Oval,_lisp,llvmo::MDString_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::MDString_O>::id,llvmo::MDString_O::static_classSymbol());
    llvmo::MDString_O::___staticClass = classllvmo__MDString_Oval;
#ifdef USE_MPS
    llvmo::MDString_O::static_Kind = gctools::GCKind<llvmo::MDString_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__MDString_Oval,llvmo::MDString_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::MDString_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::MDString_O>>::allocateClass();
        llvmo::MDString_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::MDString_O::static_className() % (void*)(llvmo::MDString_O::static_allocator) );
    classllvmo__MDString_Oval->setCreator(llvmo::MDString_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::MDString_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__MDString_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ModulePass_Oval]"));
    core::BuiltInClass_sp classllvmo__ModulePass_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ModulePass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ModulePass_Oval,_lisp,llvmo::ModulePass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ModulePass_O>::id,llvmo::ModulePass_O::static_classSymbol());
    llvmo::ModulePass_O::___staticClass = classllvmo__ModulePass_Oval;
#ifdef USE_MPS
    llvmo::ModulePass_O::static_Kind = gctools::GCKind<llvmo::ModulePass_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ModulePass_Oval,llvmo::ModulePass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ModulePass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ModulePass_O>>::allocateClass();
        llvmo::ModulePass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ModulePass_O::static_className() % (void*)(llvmo::ModulePass_O::static_allocator) );
    classllvmo__ModulePass_Oval->setCreator(llvmo::ModulePass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ModulePass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ModulePass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__PassManager_Oval]"));
    core::BuiltInClass_sp classllvmo__PassManager_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__PassManager_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__PassManager_Oval,_lisp,llvmo::PassManager_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::PassManager_O>::id,llvmo::PassManager_O::static_classSymbol());
    llvmo::PassManager_O::___staticClass = classllvmo__PassManager_Oval;
#ifdef USE_MPS
    llvmo::PassManager_O::static_Kind = gctools::GCKind<llvmo::PassManager_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__PassManager_Oval,llvmo::PassManager_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::PassManager_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::PassManager_O>>::allocateClass();
        llvmo::PassManager_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::PassManager_O::static_className() % (void*)(llvmo::PassManager_O::static_allocator) );
    classllvmo__PassManager_Oval->setCreator(llvmo::PassManager_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::PassManager_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__PassManager_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__TargetSubtargetInfo_Oval]"));
    core::BuiltInClass_sp classllvmo__TargetSubtargetInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__TargetSubtargetInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__TargetSubtargetInfo_Oval,_lisp,llvmo::TargetSubtargetInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::TargetSubtargetInfo_O>::id,llvmo::TargetSubtargetInfo_O::static_classSymbol());
    llvmo::TargetSubtargetInfo_O::___staticClass = classllvmo__TargetSubtargetInfo_Oval;
#ifdef USE_MPS
    llvmo::TargetSubtargetInfo_O::static_Kind = gctools::GCKind<llvmo::TargetSubtargetInfo_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__TargetSubtargetInfo_Oval,llvmo::TargetSubtargetInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::TargetSubtargetInfo_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::TargetSubtargetInfo_O>>::allocateClass();
        llvmo::TargetSubtargetInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::TargetSubtargetInfo_O::static_className() % (void*)(llvmo::TargetSubtargetInfo_O::static_allocator) );
    classllvmo__TargetSubtargetInfo_Oval->setCreator(llvmo::TargetSubtargetInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::TargetSubtargetInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__TargetSubtargetInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__User_Oval]"));
    core::BuiltInClass_sp classllvmo__User_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__User_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__User_Oval,_lisp,llvmo::User_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::User_O>::id,llvmo::User_O::static_classSymbol());
    llvmo::User_O::___staticClass = classllvmo__User_Oval;
#ifdef USE_MPS
    llvmo::User_O::static_Kind = gctools::GCKind<llvmo::User_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__User_Oval,llvmo::User_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::User_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::User_O>>::allocateClass();
        llvmo::User_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::User_O::static_className() % (void*)(llvmo::User_O::static_allocator) );
    classllvmo__User_Oval->setCreator(llvmo::User_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::User_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__User_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ValueAsMetadata_Oval]"));
    core::BuiltInClass_sp classllvmo__ValueAsMetadata_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ValueAsMetadata_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ValueAsMetadata_Oval,_lisp,llvmo::ValueAsMetadata_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ValueAsMetadata_O>::id,llvmo::ValueAsMetadata_O::static_classSymbol());
    llvmo::ValueAsMetadata_O::___staticClass = classllvmo__ValueAsMetadata_Oval;
#ifdef USE_MPS
    llvmo::ValueAsMetadata_O::static_Kind = gctools::GCKind<llvmo::ValueAsMetadata_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ValueAsMetadata_Oval,llvmo::ValueAsMetadata_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ValueAsMetadata_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ValueAsMetadata_O>>::allocateClass();
        llvmo::ValueAsMetadata_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ValueAsMetadata_O::static_className() % (void*)(llvmo::ValueAsMetadata_O::static_allocator) );
    classllvmo__ValueAsMetadata_Oval->setCreator(llvmo::ValueAsMetadata_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ValueAsMetadata_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ValueAsMetadata_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_UnitsPkg
    LOG(BF("Creating class[classunits__NamedUnit_Oval]"));
    core::BuiltInClass_sp classunits__NamedUnit_Oval = core::BuiltInClass_O::createUncollectable();
    classunits__NamedUnit_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__NamedUnit_Oval,_lisp,units::NamedUnit_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<units::NamedUnit_O>::id,units::NamedUnit_O::static_classSymbol());
    units::NamedUnit_O::___staticClass = classunits__NamedUnit_Oval;
#ifdef USE_MPS
    units::NamedUnit_O::static_Kind = gctools::GCKind<units::NamedUnit_O>::Kind;
#endif
    core::af_setf_findClass(classunits__NamedUnit_Oval,units::NamedUnit_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<units::NamedUnit_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<units::NamedUnit_O>>::allocateClass();
        units::NamedUnit_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::NamedUnit_O::static_className() % (void*)(units::NamedUnit_O::static_allocator) );
    classunits__NamedUnit_Oval->setCreator(units::NamedUnit_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::NamedUnit_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__NamedUnit_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__BlockEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__BlockEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__BlockEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__BlockEnvironment_Oval,_lisp,core::BlockEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::BlockEnvironment_O>::id,core::BlockEnvironment_O::static_classSymbol());
    core::BlockEnvironment_O::___staticClass = classcore__BlockEnvironment_Oval;
#ifdef USE_MPS
    core::BlockEnvironment_O::static_Kind = gctools::GCKind<core::BlockEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__BlockEnvironment_Oval,core::BlockEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::BlockEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::BlockEnvironment_O>>::allocateClass();
        core::BlockEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::BlockEnvironment_O::static_className() % (void*)(core::BlockEnvironment_O::static_allocator) );
    classcore__BlockEnvironment_Oval->setCreator(core::BlockEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::BlockEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__BlockEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__CatchEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__CatchEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__CatchEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__CatchEnvironment_Oval,_lisp,core::CatchEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::CatchEnvironment_O>::id,core::CatchEnvironment_O::static_classSymbol());
    core::CatchEnvironment_O::___staticClass = classcore__CatchEnvironment_Oval;
#ifdef USE_MPS
    core::CatchEnvironment_O::static_Kind = gctools::GCKind<core::CatchEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__CatchEnvironment_Oval,core::CatchEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::CatchEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::CatchEnvironment_O>>::allocateClass();
        core::CatchEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::CatchEnvironment_O::static_className() % (void*)(core::CatchEnvironment_O::static_allocator) );
    classcore__CatchEnvironment_Oval->setCreator(core::CatchEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::CatchEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__CatchEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Class_Oval]"));
    StandardClass_sp classcore__Class_Oval = StandardClass_O::createUncollectable();
    classcore__Class_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Class_Oval,_lisp,core::Class_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Class_O>::id,core::Class_O::static_classSymbol());
    core::Class_O::___staticClass = classcore__Class_Oval;
#ifdef USE_MPS
    core::Class_O::static_Kind = gctools::GCKind<core::Class_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Class_Oval,core::Class_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Class_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Class_O>>::allocateClass();
        core::Class_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Class_O::static_className() % (void*)(core::Class_O::static_allocator) );
    classcore__Class_Oval->setCreator(core::Class_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Class_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Class_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__DoubleFloat_Oval]"));
    core::BuiltInClass_sp classcore__DoubleFloat_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__DoubleFloat_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__DoubleFloat_Oval,_lisp,core::DoubleFloat_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::DoubleFloat_O>::id,core::DoubleFloat_O::static_classSymbol());
    core::DoubleFloat_O::___staticClass = classcore__DoubleFloat_Oval;
#ifdef USE_MPS
    core::DoubleFloat_O::static_Kind = gctools::GCKind<core::DoubleFloat_O>::Kind;
#endif
    core::af_setf_findClass(classcore__DoubleFloat_Oval,core::DoubleFloat_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::DoubleFloat_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::DoubleFloat_O>>::allocateClass();
        core::DoubleFloat_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::DoubleFloat_O::static_className() % (void*)(core::DoubleFloat_O::static_allocator) );
    classcore__DoubleFloat_Oval->setCreator(core::DoubleFloat_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::DoubleFloat_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__DoubleFloat_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__FunctionContainerEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__FunctionContainerEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__FunctionContainerEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__FunctionContainerEnvironment_Oval,_lisp,core::FunctionContainerEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::FunctionContainerEnvironment_O>::id,core::FunctionContainerEnvironment_O::static_classSymbol());
    core::FunctionContainerEnvironment_O::___staticClass = classcore__FunctionContainerEnvironment_Oval;
#ifdef USE_MPS
    core::FunctionContainerEnvironment_O::static_Kind = gctools::GCKind<core::FunctionContainerEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__FunctionContainerEnvironment_Oval,core::FunctionContainerEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::FunctionContainerEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::FunctionContainerEnvironment_O>>::allocateClass();
        core::FunctionContainerEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::FunctionContainerEnvironment_O::static_className() % (void*)(core::FunctionContainerEnvironment_O::static_allocator) );
    classcore__FunctionContainerEnvironment_Oval->setCreator(core::FunctionContainerEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::FunctionContainerEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__FunctionContainerEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__FunctionValueEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__FunctionValueEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__FunctionValueEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__FunctionValueEnvironment_Oval,_lisp,core::FunctionValueEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::FunctionValueEnvironment_O>::id,core::FunctionValueEnvironment_O::static_classSymbol());
    core::FunctionValueEnvironment_O::___staticClass = classcore__FunctionValueEnvironment_Oval;
#ifdef USE_MPS
    core::FunctionValueEnvironment_O::static_Kind = gctools::GCKind<core::FunctionValueEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__FunctionValueEnvironment_Oval,core::FunctionValueEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::FunctionValueEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::FunctionValueEnvironment_O>>::allocateClass();
        core::FunctionValueEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::FunctionValueEnvironment_O::static_className() % (void*)(core::FunctionValueEnvironment_O::static_allocator) );
    classcore__FunctionValueEnvironment_Oval->setCreator(core::FunctionValueEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::FunctionValueEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__FunctionValueEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__IOFileStream_Oval]"));
    core::BuiltInClass_sp classcore__IOFileStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__IOFileStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__IOFileStream_Oval,_lisp,core::IOFileStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::IOFileStream_O>::id,core::IOFileStream_O::static_classSymbol());
    core::IOFileStream_O::___staticClass = classcore__IOFileStream_Oval;
#ifdef USE_MPS
    core::IOFileStream_O::static_Kind = gctools::GCKind<core::IOFileStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__IOFileStream_Oval,core::IOFileStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::IOFileStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::IOFileStream_O>>::allocateClass();
        core::IOFileStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::IOFileStream_O::static_className() % (void*)(core::IOFileStream_O::static_allocator) );
    classcore__IOFileStream_Oval->setCreator(core::IOFileStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::IOFileStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__IOFileStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__IOStreamStream_Oval]"));
    core::BuiltInClass_sp classcore__IOStreamStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__IOStreamStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__IOStreamStream_Oval,_lisp,core::IOStreamStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::IOStreamStream_O>::id,core::IOStreamStream_O::static_classSymbol());
    core::IOStreamStream_O::___staticClass = classcore__IOStreamStream_Oval;
#ifdef USE_MPS
    core::IOStreamStream_O::static_Kind = gctools::GCKind<core::IOStreamStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__IOStreamStream_Oval,core::IOStreamStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::IOStreamStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::IOStreamStream_O>>::allocateClass();
        core::IOStreamStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::IOStreamStream_O::static_className() % (void*)(core::IOStreamStream_O::static_allocator) );
    classcore__IOStreamStream_Oval->setCreator(core::IOStreamStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::IOStreamStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__IOStreamStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Integer_Oval]"));
    core::BuiltInClass_sp classcore__Integer_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Integer_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Integer_Oval,_lisp,core::Integer_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Integer_O>::id,core::Integer_O::static_classSymbol());
    core::Integer_O::___staticClass = classcore__Integer_Oval;
#ifdef USE_MPS
    core::Integer_O::static_Kind = gctools::GCKind<core::Integer_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Integer_Oval,core::Integer_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Integer_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Integer_O>>::allocateClass();
        core::Integer_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Integer_O::static_className() % (void*)(core::Integer_O::static_allocator) );
    classcore__Integer_Oval->setCreator(core::Integer_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Integer_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Integer_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__LongFloat_Oval]"));
    core::BuiltInClass_sp classcore__LongFloat_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__LongFloat_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__LongFloat_Oval,_lisp,core::LongFloat_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::LongFloat_O>::id,core::LongFloat_O::static_classSymbol());
    core::LongFloat_O::___staticClass = classcore__LongFloat_Oval;
#ifdef USE_MPS
    core::LongFloat_O::static_Kind = gctools::GCKind<core::LongFloat_O>::Kind;
#endif
    core::af_setf_findClass(classcore__LongFloat_Oval,core::LongFloat_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::LongFloat_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::LongFloat_O>>::allocateClass();
        core::LongFloat_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::LongFloat_O::static_className() % (void*)(core::LongFloat_O::static_allocator) );
    classcore__LongFloat_Oval->setCreator(core::LongFloat_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::LongFloat_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__LongFloat_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__MacroletEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__MacroletEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__MacroletEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__MacroletEnvironment_Oval,_lisp,core::MacroletEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::MacroletEnvironment_O>::id,core::MacroletEnvironment_O::static_classSymbol());
    core::MacroletEnvironment_O::___staticClass = classcore__MacroletEnvironment_Oval;
#ifdef USE_MPS
    core::MacroletEnvironment_O::static_Kind = gctools::GCKind<core::MacroletEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__MacroletEnvironment_Oval,core::MacroletEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::MacroletEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::MacroletEnvironment_O>>::allocateClass();
        core::MacroletEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::MacroletEnvironment_O::static_className() % (void*)(core::MacroletEnvironment_O::static_allocator) );
    classcore__MacroletEnvironment_Oval->setCreator(core::MacroletEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::MacroletEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__MacroletEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Ratio_Oval]"));
    core::BuiltInClass_sp classcore__Ratio_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Ratio_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Ratio_Oval,_lisp,core::Ratio_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Ratio_O>::id,core::Ratio_O::static_classSymbol());
    core::Ratio_O::___staticClass = classcore__Ratio_Oval;
#ifdef USE_MPS
    core::Ratio_O::static_Kind = gctools::GCKind<core::Ratio_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Ratio_Oval,core::Ratio_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Ratio_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Ratio_O>>::allocateClass();
        core::Ratio_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Ratio_O::static_className() % (void*)(core::Ratio_O::static_allocator) );
    classcore__Ratio_Oval->setCreator(core::Ratio_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Ratio_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Ratio_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ShortFloat_Oval]"));
    core::BuiltInClass_sp classcore__ShortFloat_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ShortFloat_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ShortFloat_Oval,_lisp,core::ShortFloat_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ShortFloat_O>::id,core::ShortFloat_O::static_classSymbol());
    core::ShortFloat_O::___staticClass = classcore__ShortFloat_Oval;
#ifdef USE_MPS
    core::ShortFloat_O::static_Kind = gctools::GCKind<core::ShortFloat_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ShortFloat_Oval,core::ShortFloat_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ShortFloat_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ShortFloat_O>>::allocateClass();
        core::ShortFloat_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ShortFloat_O::static_className() % (void*)(core::ShortFloat_O::static_allocator) );
    classcore__ShortFloat_Oval->setCreator(core::ShortFloat_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ShortFloat_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ShortFloat_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SimpleBitVector_Oval]"));
    core::BuiltInClass_sp classcore__SimpleBitVector_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SimpleBitVector_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SimpleBitVector_Oval,_lisp,core::SimpleBitVector_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SimpleBitVector_O>::id,core::SimpleBitVector_O::static_classSymbol());
    core::SimpleBitVector_O::___staticClass = classcore__SimpleBitVector_Oval;
#ifdef USE_MPS
    core::SimpleBitVector_O::static_Kind = gctools::GCKind<core::SimpleBitVector_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SimpleBitVector_Oval,core::SimpleBitVector_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SimpleBitVector_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SimpleBitVector_O>>::allocateClass();
        core::SimpleBitVector_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SimpleBitVector_O::static_className() % (void*)(core::SimpleBitVector_O::static_allocator) );
    classcore__SimpleBitVector_Oval->setCreator(core::SimpleBitVector_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SimpleBitVector_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SimpleBitVector_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SingleFloat_dummy_Oval]"));
    core::BuiltInClass_sp classcore__SingleFloat_dummy_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SingleFloat_dummy_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SingleFloat_dummy_Oval,_lisp,core::SingleFloat_dummy_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SingleFloat_dummy_O>::id,core::SingleFloat_dummy_O::static_classSymbol());
    core::SingleFloat_dummy_O::___staticClass = classcore__SingleFloat_dummy_Oval;
#ifdef USE_MPS
    core::SingleFloat_dummy_O::static_Kind = gctools::GCKind<core::SingleFloat_dummy_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SingleFloat_dummy_Oval,core::SingleFloat_dummy_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SingleFloat_dummy_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SingleFloat_dummy_O>>::allocateClass();
        core::SingleFloat_dummy_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SingleFloat_dummy_O::static_className() % (void*)(core::SingleFloat_dummy_O::static_allocator) );
    classcore__SingleFloat_dummy_Oval->setCreator(core::SingleFloat_dummy_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SingleFloat_dummy_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SingleFloat_dummy_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StackValueEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__StackValueEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__StackValueEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StackValueEnvironment_Oval,_lisp,core::StackValueEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StackValueEnvironment_O>::id,core::StackValueEnvironment_O::static_classSymbol());
    core::StackValueEnvironment_O::___staticClass = classcore__StackValueEnvironment_Oval;
#ifdef USE_MPS
    core::StackValueEnvironment_O::static_Kind = gctools::GCKind<core::StackValueEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StackValueEnvironment_Oval,core::StackValueEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StackValueEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StackValueEnvironment_O>>::allocateClass();
        core::StackValueEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StackValueEnvironment_O::static_className() % (void*)(core::StackValueEnvironment_O::static_allocator) );
    classcore__StackValueEnvironment_Oval->setCreator(core::StackValueEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StackValueEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StackValueEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Str_Oval]"));
    core::BuiltInClass_sp classcore__Str_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Str_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Str_Oval,_lisp,core::Str_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Str_O>::id,core::Str_O::static_classSymbol());
    core::Str_O::___staticClass = classcore__Str_Oval;
#ifdef USE_MPS
    core::Str_O::static_Kind = gctools::GCKind<core::Str_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Str_Oval,core::Str_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Str_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Str_O>>::allocateClass();
        core::Str_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Str_O::static_className() % (void*)(core::Str_O::static_allocator) );
    classcore__Str_Oval->setCreator(core::Str_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Str_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Str_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StringInputStream_Oval]"));
    core::BuiltInClass_sp classcore__StringInputStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__StringInputStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StringInputStream_Oval,_lisp,core::StringInputStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StringInputStream_O>::id,core::StringInputStream_O::static_classSymbol());
    core::StringInputStream_O::___staticClass = classcore__StringInputStream_Oval;
#ifdef USE_MPS
    core::StringInputStream_O::static_Kind = gctools::GCKind<core::StringInputStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StringInputStream_Oval,core::StringInputStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StringInputStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StringInputStream_O>>::allocateClass();
        core::StringInputStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StringInputStream_O::static_className() % (void*)(core::StringInputStream_O::static_allocator) );
    classcore__StringInputStream_Oval->setCreator(core::StringInputStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StringInputStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StringInputStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StringOutputStream_Oval]"));
    core::BuiltInClass_sp classcore__StringOutputStream_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__StringOutputStream_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StringOutputStream_Oval,_lisp,core::StringOutputStream_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StringOutputStream_O>::id,core::StringOutputStream_O::static_classSymbol());
    core::StringOutputStream_O::___staticClass = classcore__StringOutputStream_Oval;
#ifdef USE_MPS
    core::StringOutputStream_O::static_Kind = gctools::GCKind<core::StringOutputStream_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StringOutputStream_Oval,core::StringOutputStream_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StringOutputStream_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StringOutputStream_O>>::allocateClass();
        core::StringOutputStream_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StringOutputStream_O::static_className() % (void*)(core::StringOutputStream_O::static_allocator) );
    classcore__StringOutputStream_Oval->setCreator(core::StringOutputStream_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StringOutputStream_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StringOutputStream_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__SymbolMacroletEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__SymbolMacroletEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__SymbolMacroletEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__SymbolMacroletEnvironment_Oval,_lisp,core::SymbolMacroletEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::SymbolMacroletEnvironment_O>::id,core::SymbolMacroletEnvironment_O::static_classSymbol());
    core::SymbolMacroletEnvironment_O::___staticClass = classcore__SymbolMacroletEnvironment_Oval;
#ifdef USE_MPS
    core::SymbolMacroletEnvironment_O::static_Kind = gctools::GCKind<core::SymbolMacroletEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__SymbolMacroletEnvironment_Oval,core::SymbolMacroletEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::SymbolMacroletEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::SymbolMacroletEnvironment_O>>::allocateClass();
        core::SymbolMacroletEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::SymbolMacroletEnvironment_O::static_className() % (void*)(core::SymbolMacroletEnvironment_O::static_allocator) );
    classcore__SymbolMacroletEnvironment_Oval->setCreator(core::SymbolMacroletEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::SymbolMacroletEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__SymbolMacroletEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__TagbodyEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__TagbodyEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__TagbodyEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__TagbodyEnvironment_Oval,_lisp,core::TagbodyEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::TagbodyEnvironment_O>::id,core::TagbodyEnvironment_O::static_classSymbol());
    core::TagbodyEnvironment_O::___staticClass = classcore__TagbodyEnvironment_Oval;
#ifdef USE_MPS
    core::TagbodyEnvironment_O::static_Kind = gctools::GCKind<core::TagbodyEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__TagbodyEnvironment_Oval,core::TagbodyEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::TagbodyEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::TagbodyEnvironment_O>>::allocateClass();
        core::TagbodyEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::TagbodyEnvironment_O::static_className() % (void*)(core::TagbodyEnvironment_O::static_allocator) );
    classcore__TagbodyEnvironment_Oval->setCreator(core::TagbodyEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::TagbodyEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__TagbodyEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__UnwindProtectEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__UnwindProtectEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__UnwindProtectEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__UnwindProtectEnvironment_Oval,_lisp,core::UnwindProtectEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::UnwindProtectEnvironment_O>::id,core::UnwindProtectEnvironment_O::static_classSymbol());
    core::UnwindProtectEnvironment_O::___staticClass = classcore__UnwindProtectEnvironment_Oval;
#ifdef USE_MPS
    core::UnwindProtectEnvironment_O::static_Kind = gctools::GCKind<core::UnwindProtectEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__UnwindProtectEnvironment_Oval,core::UnwindProtectEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::UnwindProtectEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::UnwindProtectEnvironment_O>>::allocateClass();
        core::UnwindProtectEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::UnwindProtectEnvironment_O::static_className() % (void*)(core::UnwindProtectEnvironment_O::static_allocator) );
    classcore__UnwindProtectEnvironment_Oval->setCreator(core::UnwindProtectEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::UnwindProtectEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__UnwindProtectEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ValueEnvironment_Oval]"));
    core::BuiltInClass_sp classcore__ValueEnvironment_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__ValueEnvironment_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ValueEnvironment_Oval,_lisp,core::ValueEnvironment_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ValueEnvironment_O>::id,core::ValueEnvironment_O::static_classSymbol());
    core::ValueEnvironment_O::___staticClass = classcore__ValueEnvironment_Oval;
#ifdef USE_MPS
    core::ValueEnvironment_O::static_Kind = gctools::GCKind<core::ValueEnvironment_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ValueEnvironment_Oval,core::ValueEnvironment_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ValueEnvironment_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ValueEnvironment_O>>::allocateClass();
        core::ValueEnvironment_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ValueEnvironment_O::static_className() % (void*)(core::ValueEnvironment_O::static_allocator) );
    classcore__ValueEnvironment_Oval->setCreator(core::ValueEnvironment_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ValueEnvironment_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ValueEnvironment_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__VectorObjectsWithFillPtr_Oval]"));
    core::BuiltInClass_sp classcore__VectorObjectsWithFillPtr_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__VectorObjectsWithFillPtr_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__VectorObjectsWithFillPtr_Oval,_lisp,core::VectorObjectsWithFillPtr_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::VectorObjectsWithFillPtr_O>::id,core::VectorObjectsWithFillPtr_O::static_classSymbol());
    core::VectorObjectsWithFillPtr_O::___staticClass = classcore__VectorObjectsWithFillPtr_Oval;
#ifdef USE_MPS
    core::VectorObjectsWithFillPtr_O::static_Kind = gctools::GCKind<core::VectorObjectsWithFillPtr_O>::Kind;
#endif
    core::af_setf_findClass(classcore__VectorObjectsWithFillPtr_Oval,core::VectorObjectsWithFillPtr_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::VectorObjectsWithFillPtr_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::VectorObjectsWithFillPtr_O>>::allocateClass();
        core::VectorObjectsWithFillPtr_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::VectorObjectsWithFillPtr_O::static_className() % (void*)(core::VectorObjectsWithFillPtr_O::static_allocator) );
    classcore__VectorObjectsWithFillPtr_Oval->setCreator(core::VectorObjectsWithFillPtr_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::VectorObjectsWithFillPtr_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__VectorObjectsWithFillPtr_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_GeomPkg
    LOG(BF("Creating class[classgeom__CoordinateArrayWithHash_Oval]"));
    core::BuiltInClass_sp classgeom__CoordinateArrayWithHash_Oval = core::BuiltInClass_O::createUncollectable();
    classgeom__CoordinateArrayWithHash_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classgeom__CoordinateArrayWithHash_Oval,_lisp,geom::CoordinateArrayWithHash_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<geom::CoordinateArrayWithHash_O>::id,geom::CoordinateArrayWithHash_O::static_classSymbol());
    geom::CoordinateArrayWithHash_O::___staticClass = classgeom__CoordinateArrayWithHash_Oval;
#ifdef USE_MPS
    geom::CoordinateArrayWithHash_O::static_Kind = gctools::GCKind<geom::CoordinateArrayWithHash_O>::Kind;
#endif
    core::af_setf_findClass(classgeom__CoordinateArrayWithHash_Oval,geom::CoordinateArrayWithHash_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<geom::CoordinateArrayWithHash_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<geom::CoordinateArrayWithHash_O>>::allocateClass();
        geom::CoordinateArrayWithHash_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::CoordinateArrayWithHash_O::static_className() % (void*)(geom::CoordinateArrayWithHash_O::static_allocator) );
    classgeom__CoordinateArrayWithHash_Oval->setCreator(geom::CoordinateArrayWithHash_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::CoordinateArrayWithHash_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__CoordinateArrayWithHash_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Constant_Oval]"));
    core::BuiltInClass_sp classllvmo__Constant_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Constant_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Constant_Oval,_lisp,llvmo::Constant_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Constant_O>::id,llvmo::Constant_O::static_classSymbol());
    llvmo::Constant_O::___staticClass = classllvmo__Constant_Oval;
#ifdef USE_MPS
    llvmo::Constant_O::static_Kind = gctools::GCKind<llvmo::Constant_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Constant_Oval,llvmo::Constant_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Constant_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Constant_O>>::allocateClass();
        llvmo::Constant_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Constant_O::static_className() % (void*)(llvmo::Constant_O::static_allocator) );
    classllvmo__Constant_Oval->setCreator(llvmo::Constant_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Constant_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Constant_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ImmutablePass_Oval]"));
    core::BuiltInClass_sp classllvmo__ImmutablePass_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ImmutablePass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ImmutablePass_Oval,_lisp,llvmo::ImmutablePass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ImmutablePass_O>::id,llvmo::ImmutablePass_O::static_classSymbol());
    llvmo::ImmutablePass_O::___staticClass = classllvmo__ImmutablePass_Oval;
#ifdef USE_MPS
    llvmo::ImmutablePass_O::static_Kind = gctools::GCKind<llvmo::ImmutablePass_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ImmutablePass_Oval,llvmo::ImmutablePass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ImmutablePass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ImmutablePass_O>>::allocateClass();
        llvmo::ImmutablePass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ImmutablePass_O::static_className() % (void*)(llvmo::ImmutablePass_O::static_allocator) );
    classllvmo__ImmutablePass_Oval->setCreator(llvmo::ImmutablePass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ImmutablePass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ImmutablePass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Instruction_Oval]"));
    core::BuiltInClass_sp classllvmo__Instruction_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Instruction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Instruction_Oval,_lisp,llvmo::Instruction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Instruction_O>::id,llvmo::Instruction_O::static_classSymbol());
    llvmo::Instruction_O::___staticClass = classllvmo__Instruction_Oval;
#ifdef USE_MPS
    llvmo::Instruction_O::static_Kind = gctools::GCKind<llvmo::Instruction_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Instruction_Oval,llvmo::Instruction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Instruction_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Instruction_O>>::allocateClass();
        llvmo::Instruction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Instruction_O::static_className() % (void*)(llvmo::Instruction_O::static_allocator) );
    classllvmo__Instruction_Oval->setCreator(llvmo::Instruction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Instruction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Instruction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__SequentialType_Oval]"));
    core::BuiltInClass_sp classllvmo__SequentialType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__SequentialType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__SequentialType_Oval,_lisp,llvmo::SequentialType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::SequentialType_O>::id,llvmo::SequentialType_O::static_classSymbol());
    llvmo::SequentialType_O::___staticClass = classllvmo__SequentialType_Oval;
#ifdef USE_MPS
    llvmo::SequentialType_O::static_Kind = gctools::GCKind<llvmo::SequentialType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__SequentialType_Oval,llvmo::SequentialType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::SequentialType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::SequentialType_O>>::allocateClass();
        llvmo::SequentialType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::SequentialType_O::static_className() % (void*)(llvmo::SequentialType_O::static_allocator) );
    classllvmo__SequentialType_Oval->setCreator(llvmo::SequentialType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::SequentialType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__SequentialType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__StructType_Oval]"));
    core::BuiltInClass_sp classllvmo__StructType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__StructType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__StructType_Oval,_lisp,llvmo::StructType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::StructType_O>::id,llvmo::StructType_O::static_classSymbol());
    llvmo::StructType_O::___staticClass = classllvmo__StructType_Oval;
#ifdef USE_MPS
    llvmo::StructType_O::static_Kind = gctools::GCKind<llvmo::StructType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__StructType_Oval,llvmo::StructType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::StructType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::StructType_O>>::allocateClass();
        llvmo::StructType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::StructType_O::static_className() % (void*)(llvmo::StructType_O::static_allocator) );
    classllvmo__StructType_Oval->setCreator(llvmo::StructType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::StructType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__StructType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Bignum_Oval]"));
    core::BuiltInClass_sp classcore__Bignum_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Bignum_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Bignum_Oval,_lisp,core::Bignum_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Bignum_O>::id,core::Bignum_O::static_classSymbol());
    core::Bignum_O::___staticClass = classcore__Bignum_Oval;
#ifdef USE_MPS
    core::Bignum_O::static_Kind = gctools::GCKind<core::Bignum_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Bignum_Oval,core::Bignum_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Bignum_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Bignum_O>>::allocateClass();
        core::Bignum_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Bignum_O::static_className() % (void*)(core::Bignum_O::static_allocator) );
    classcore__Bignum_Oval->setCreator(core::Bignum_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Bignum_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Bignum_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__BuiltInClass_Oval]"));
    StandardClass_sp classcore__BuiltInClass_Oval = StandardClass_O::createUncollectable();
    classcore__BuiltInClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__BuiltInClass_Oval,_lisp,core::BuiltInClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::BuiltInClass_O>::id,core::BuiltInClass_O::static_classSymbol());
    core::BuiltInClass_O::___staticClass = classcore__BuiltInClass_Oval;
#ifdef USE_MPS
    core::BuiltInClass_O::static_Kind = gctools::GCKind<core::BuiltInClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__BuiltInClass_Oval,core::BuiltInClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::BuiltInClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::BuiltInClass_O>>::allocateClass();
        core::BuiltInClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::BuiltInClass_O::static_className() % (void*)(core::BuiltInClass_O::static_allocator) );
    classcore__BuiltInClass_Oval->setCreator(core::BuiltInClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::BuiltInClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__BuiltInClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__CxxClass_Oval]"));
    StandardClass_sp classcore__CxxClass_Oval = StandardClass_O::createUncollectable();
    classcore__CxxClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__CxxClass_Oval,_lisp,core::CxxClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::CxxClass_O>::id,core::CxxClass_O::static_classSymbol());
    core::CxxClass_O::___staticClass = classcore__CxxClass_Oval;
#ifdef USE_MPS
    core::CxxClass_O::static_Kind = gctools::GCKind<core::CxxClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__CxxClass_Oval,core::CxxClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::CxxClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::CxxClass_O>>::allocateClass();
        core::CxxClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::CxxClass_O::static_className() % (void*)(core::CxxClass_O::static_allocator) );
    classcore__CxxClass_Oval->setCreator(core::CxxClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::CxxClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__CxxClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__Fixnum_dummy_Oval]"));
    core::BuiltInClass_sp classcore__Fixnum_dummy_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__Fixnum_dummy_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__Fixnum_dummy_Oval,_lisp,core::Fixnum_dummy_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::Fixnum_dummy_O>::id,core::Fixnum_dummy_O::static_classSymbol());
    core::Fixnum_dummy_O::___staticClass = classcore__Fixnum_dummy_Oval;
#ifdef USE_MPS
    core::Fixnum_dummy_O::static_Kind = gctools::GCKind<core::Fixnum_dummy_O>::Kind;
#endif
    core::af_setf_findClass(classcore__Fixnum_dummy_Oval,core::Fixnum_dummy_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::Fixnum_dummy_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::Fixnum_dummy_O>>::allocateClass();
        core::Fixnum_dummy_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::Fixnum_dummy_O::static_className() % (void*)(core::Fixnum_dummy_O::static_allocator) );
    classcore__Fixnum_dummy_Oval->setCreator(core::Fixnum_dummy_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::Fixnum_dummy_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__Fixnum_dummy_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__ForwardReferencedClass_Oval]"));
    StandardClass_sp classcore__ForwardReferencedClass_Oval = StandardClass_O::createUncollectable();
    classcore__ForwardReferencedClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__ForwardReferencedClass_Oval,_lisp,core::ForwardReferencedClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::ForwardReferencedClass_O>::id,core::ForwardReferencedClass_O::static_classSymbol());
    core::ForwardReferencedClass_O::___staticClass = classcore__ForwardReferencedClass_Oval;
#ifdef USE_MPS
    core::ForwardReferencedClass_O::static_Kind = gctools::GCKind<core::ForwardReferencedClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__ForwardReferencedClass_Oval,core::ForwardReferencedClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::ForwardReferencedClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::ForwardReferencedClass_O>>::allocateClass();
        core::ForwardReferencedClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::ForwardReferencedClass_O::static_className() % (void*)(core::ForwardReferencedClass_O::static_allocator) );
    classcore__ForwardReferencedClass_Oval->setCreator(core::ForwardReferencedClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::ForwardReferencedClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__ForwardReferencedClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StdClass_Oval]"));
    StandardClass_sp classcore__StdClass_Oval = StandardClass_O::createUncollectable();
    classcore__StdClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StdClass_Oval,_lisp,core::StdClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StdClass_O>::id,core::StdClass_O::static_classSymbol());
    core::StdClass_O::___staticClass = classcore__StdClass_Oval;
#ifdef USE_MPS
    core::StdClass_O::static_Kind = gctools::GCKind<core::StdClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StdClass_Oval,core::StdClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StdClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StdClass_O>>::allocateClass();
        core::StdClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StdClass_O::static_className() % (void*)(core::StdClass_O::static_allocator) );
    classcore__StdClass_Oval->setCreator(core::StdClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StdClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StdClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StrWithFillPtr_Oval]"));
    core::BuiltInClass_sp classcore__StrWithFillPtr_Oval = core::BuiltInClass_O::createUncollectable();
    classcore__StrWithFillPtr_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StrWithFillPtr_Oval,_lisp,core::StrWithFillPtr_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StrWithFillPtr_O>::id,core::StrWithFillPtr_O::static_classSymbol());
    core::StrWithFillPtr_O::___staticClass = classcore__StrWithFillPtr_Oval;
#ifdef USE_MPS
    core::StrWithFillPtr_O::static_Kind = gctools::GCKind<core::StrWithFillPtr_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StrWithFillPtr_Oval,core::StrWithFillPtr_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StrWithFillPtr_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StrWithFillPtr_O>>::allocateClass();
        core::StrWithFillPtr_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StrWithFillPtr_O::static_className() % (void*)(core::StrWithFillPtr_O::static_allocator) );
    classcore__StrWithFillPtr_Oval->setCreator(core::StrWithFillPtr_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StrWithFillPtr_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StrWithFillPtr_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StructureClass_Oval]"));
    StandardClass_sp classcore__StructureClass_Oval = StandardClass_O::createUncollectable();
    classcore__StructureClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StructureClass_Oval,_lisp,core::StructureClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StructureClass_O>::id,core::StructureClass_O::static_classSymbol());
    core::StructureClass_O::___staticClass = classcore__StructureClass_Oval;
#ifdef USE_MPS
    core::StructureClass_O::static_Kind = gctools::GCKind<core::StructureClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StructureClass_Oval,core::StructureClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StructureClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StructureClass_O>>::allocateClass();
        core::StructureClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StructureClass_O::static_className() % (void*)(core::StructureClass_O::static_allocator) );
    classcore__StructureClass_Oval->setCreator(core::StructureClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StructureClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StructureClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ArrayType_Oval]"));
    core::BuiltInClass_sp classllvmo__ArrayType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ArrayType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ArrayType_Oval,_lisp,llvmo::ArrayType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ArrayType_O>::id,llvmo::ArrayType_O::static_classSymbol());
    llvmo::ArrayType_O::___staticClass = classllvmo__ArrayType_Oval;
#ifdef USE_MPS
    llvmo::ArrayType_O::static_Kind = gctools::GCKind<llvmo::ArrayType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ArrayType_Oval,llvmo::ArrayType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ArrayType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ArrayType_O>>::allocateClass();
        llvmo::ArrayType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ArrayType_O::static_className() % (void*)(llvmo::ArrayType_O::static_allocator) );
    classllvmo__ArrayType_Oval->setCreator(llvmo::ArrayType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ArrayType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ArrayType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__AtomicCmpXchgInst_Oval]"));
    core::BuiltInClass_sp classllvmo__AtomicCmpXchgInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__AtomicCmpXchgInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__AtomicCmpXchgInst_Oval,_lisp,llvmo::AtomicCmpXchgInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::AtomicCmpXchgInst_O>::id,llvmo::AtomicCmpXchgInst_O::static_classSymbol());
    llvmo::AtomicCmpXchgInst_O::___staticClass = classllvmo__AtomicCmpXchgInst_Oval;
#ifdef USE_MPS
    llvmo::AtomicCmpXchgInst_O::static_Kind = gctools::GCKind<llvmo::AtomicCmpXchgInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__AtomicCmpXchgInst_Oval,llvmo::AtomicCmpXchgInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::AtomicCmpXchgInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::AtomicCmpXchgInst_O>>::allocateClass();
        llvmo::AtomicCmpXchgInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::AtomicCmpXchgInst_O::static_className() % (void*)(llvmo::AtomicCmpXchgInst_O::static_allocator) );
    classllvmo__AtomicCmpXchgInst_Oval->setCreator(llvmo::AtomicCmpXchgInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::AtomicCmpXchgInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__AtomicCmpXchgInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__AtomicRMWInst_Oval]"));
    core::BuiltInClass_sp classllvmo__AtomicRMWInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__AtomicRMWInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__AtomicRMWInst_Oval,_lisp,llvmo::AtomicRMWInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::AtomicRMWInst_O>::id,llvmo::AtomicRMWInst_O::static_classSymbol());
    llvmo::AtomicRMWInst_O::___staticClass = classllvmo__AtomicRMWInst_Oval;
#ifdef USE_MPS
    llvmo::AtomicRMWInst_O::static_Kind = gctools::GCKind<llvmo::AtomicRMWInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__AtomicRMWInst_Oval,llvmo::AtomicRMWInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::AtomicRMWInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::AtomicRMWInst_O>>::allocateClass();
        llvmo::AtomicRMWInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::AtomicRMWInst_O::static_className() % (void*)(llvmo::AtomicRMWInst_O::static_allocator) );
    classllvmo__AtomicRMWInst_Oval->setCreator(llvmo::AtomicRMWInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::AtomicRMWInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__AtomicRMWInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__BlockAddress_Oval]"));
    core::BuiltInClass_sp classllvmo__BlockAddress_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__BlockAddress_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__BlockAddress_Oval,_lisp,llvmo::BlockAddress_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::BlockAddress_O>::id,llvmo::BlockAddress_O::static_classSymbol());
    llvmo::BlockAddress_O::___staticClass = classllvmo__BlockAddress_Oval;
#ifdef USE_MPS
    llvmo::BlockAddress_O::static_Kind = gctools::GCKind<llvmo::BlockAddress_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__BlockAddress_Oval,llvmo::BlockAddress_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::BlockAddress_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::BlockAddress_O>>::allocateClass();
        llvmo::BlockAddress_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::BlockAddress_O::static_className() % (void*)(llvmo::BlockAddress_O::static_allocator) );
    classllvmo__BlockAddress_Oval->setCreator(llvmo::BlockAddress_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::BlockAddress_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__BlockAddress_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__CallInst_Oval]"));
    core::BuiltInClass_sp classllvmo__CallInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__CallInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__CallInst_Oval,_lisp,llvmo::CallInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::CallInst_O>::id,llvmo::CallInst_O::static_classSymbol());
    llvmo::CallInst_O::___staticClass = classllvmo__CallInst_Oval;
#ifdef USE_MPS
    llvmo::CallInst_O::static_Kind = gctools::GCKind<llvmo::CallInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__CallInst_Oval,llvmo::CallInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::CallInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::CallInst_O>>::allocateClass();
        llvmo::CallInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::CallInst_O::static_className() % (void*)(llvmo::CallInst_O::static_allocator) );
    classllvmo__CallInst_Oval->setCreator(llvmo::CallInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::CallInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__CallInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantArray_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantArray_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantArray_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantArray_Oval,_lisp,llvmo::ConstantArray_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantArray_O>::id,llvmo::ConstantArray_O::static_classSymbol());
    llvmo::ConstantArray_O::___staticClass = classllvmo__ConstantArray_Oval;
#ifdef USE_MPS
    llvmo::ConstantArray_O::static_Kind = gctools::GCKind<llvmo::ConstantArray_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantArray_Oval,llvmo::ConstantArray_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantArray_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantArray_O>>::allocateClass();
        llvmo::ConstantArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantArray_O::static_className() % (void*)(llvmo::ConstantArray_O::static_allocator) );
    classllvmo__ConstantArray_Oval->setCreator(llvmo::ConstantArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantDataSequential_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantDataSequential_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantDataSequential_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantDataSequential_Oval,_lisp,llvmo::ConstantDataSequential_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantDataSequential_O>::id,llvmo::ConstantDataSequential_O::static_classSymbol());
    llvmo::ConstantDataSequential_O::___staticClass = classllvmo__ConstantDataSequential_Oval;
#ifdef USE_MPS
    llvmo::ConstantDataSequential_O::static_Kind = gctools::GCKind<llvmo::ConstantDataSequential_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantDataSequential_Oval,llvmo::ConstantDataSequential_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantDataSequential_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantDataSequential_O>>::allocateClass();
        llvmo::ConstantDataSequential_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantDataSequential_O::static_className() % (void*)(llvmo::ConstantDataSequential_O::static_allocator) );
    classllvmo__ConstantDataSequential_Oval->setCreator(llvmo::ConstantDataSequential_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantDataSequential_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantDataSequential_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantExpr_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantExpr_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantExpr_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantExpr_Oval,_lisp,llvmo::ConstantExpr_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantExpr_O>::id,llvmo::ConstantExpr_O::static_classSymbol());
    llvmo::ConstantExpr_O::___staticClass = classllvmo__ConstantExpr_Oval;
#ifdef USE_MPS
    llvmo::ConstantExpr_O::static_Kind = gctools::GCKind<llvmo::ConstantExpr_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantExpr_Oval,llvmo::ConstantExpr_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantExpr_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantExpr_O>>::allocateClass();
        llvmo::ConstantExpr_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantExpr_O::static_className() % (void*)(llvmo::ConstantExpr_O::static_allocator) );
    classllvmo__ConstantExpr_Oval->setCreator(llvmo::ConstantExpr_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantExpr_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantExpr_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantFP_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantFP_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantFP_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantFP_Oval,_lisp,llvmo::ConstantFP_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantFP_O>::id,llvmo::ConstantFP_O::static_classSymbol());
    llvmo::ConstantFP_O::___staticClass = classllvmo__ConstantFP_Oval;
#ifdef USE_MPS
    llvmo::ConstantFP_O::static_Kind = gctools::GCKind<llvmo::ConstantFP_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantFP_Oval,llvmo::ConstantFP_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantFP_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantFP_O>>::allocateClass();
        llvmo::ConstantFP_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantFP_O::static_className() % (void*)(llvmo::ConstantFP_O::static_allocator) );
    classllvmo__ConstantFP_Oval->setCreator(llvmo::ConstantFP_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantFP_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantFP_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantInt_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantInt_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantInt_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantInt_Oval,_lisp,llvmo::ConstantInt_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantInt_O>::id,llvmo::ConstantInt_O::static_classSymbol());
    llvmo::ConstantInt_O::___staticClass = classllvmo__ConstantInt_Oval;
#ifdef USE_MPS
    llvmo::ConstantInt_O::static_Kind = gctools::GCKind<llvmo::ConstantInt_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantInt_Oval,llvmo::ConstantInt_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantInt_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantInt_O>>::allocateClass();
        llvmo::ConstantInt_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantInt_O::static_className() % (void*)(llvmo::ConstantInt_O::static_allocator) );
    classllvmo__ConstantInt_Oval->setCreator(llvmo::ConstantInt_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantInt_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantInt_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantPointerNull_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantPointerNull_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantPointerNull_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantPointerNull_Oval,_lisp,llvmo::ConstantPointerNull_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantPointerNull_O>::id,llvmo::ConstantPointerNull_O::static_classSymbol());
    llvmo::ConstantPointerNull_O::___staticClass = classllvmo__ConstantPointerNull_Oval;
#ifdef USE_MPS
    llvmo::ConstantPointerNull_O::static_Kind = gctools::GCKind<llvmo::ConstantPointerNull_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantPointerNull_Oval,llvmo::ConstantPointerNull_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantPointerNull_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantPointerNull_O>>::allocateClass();
        llvmo::ConstantPointerNull_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantPointerNull_O::static_className() % (void*)(llvmo::ConstantPointerNull_O::static_allocator) );
    classllvmo__ConstantPointerNull_Oval->setCreator(llvmo::ConstantPointerNull_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantPointerNull_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantPointerNull_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantStruct_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantStruct_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantStruct_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantStruct_Oval,_lisp,llvmo::ConstantStruct_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantStruct_O>::id,llvmo::ConstantStruct_O::static_classSymbol());
    llvmo::ConstantStruct_O::___staticClass = classllvmo__ConstantStruct_Oval;
#ifdef USE_MPS
    llvmo::ConstantStruct_O::static_Kind = gctools::GCKind<llvmo::ConstantStruct_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantStruct_Oval,llvmo::ConstantStruct_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantStruct_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantStruct_O>>::allocateClass();
        llvmo::ConstantStruct_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantStruct_O::static_className() % (void*)(llvmo::ConstantStruct_O::static_allocator) );
    classllvmo__ConstantStruct_Oval->setCreator(llvmo::ConstantStruct_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantStruct_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantStruct_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__DataLayoutPass_Oval]"));
    core::BuiltInClass_sp classllvmo__DataLayoutPass_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__DataLayoutPass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__DataLayoutPass_Oval,_lisp,llvmo::DataLayoutPass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::DataLayoutPass_O>::id,llvmo::DataLayoutPass_O::static_classSymbol());
    llvmo::DataLayoutPass_O::___staticClass = classllvmo__DataLayoutPass_Oval;
#ifdef USE_MPS
    llvmo::DataLayoutPass_O::static_Kind = gctools::GCKind<llvmo::DataLayoutPass_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__DataLayoutPass_Oval,llvmo::DataLayoutPass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::DataLayoutPass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::DataLayoutPass_O>>::allocateClass();
        llvmo::DataLayoutPass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::DataLayoutPass_O::static_className() % (void*)(llvmo::DataLayoutPass_O::static_allocator) );
    classllvmo__DataLayoutPass_Oval->setCreator(llvmo::DataLayoutPass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::DataLayoutPass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__DataLayoutPass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__FenceInst_Oval]"));
    core::BuiltInClass_sp classllvmo__FenceInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__FenceInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__FenceInst_Oval,_lisp,llvmo::FenceInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::FenceInst_O>::id,llvmo::FenceInst_O::static_classSymbol());
    llvmo::FenceInst_O::___staticClass = classllvmo__FenceInst_Oval;
#ifdef USE_MPS
    llvmo::FenceInst_O::static_Kind = gctools::GCKind<llvmo::FenceInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__FenceInst_Oval,llvmo::FenceInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::FenceInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::FenceInst_O>>::allocateClass();
        llvmo::FenceInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::FenceInst_O::static_className() % (void*)(llvmo::FenceInst_O::static_allocator) );
    classllvmo__FenceInst_Oval->setCreator(llvmo::FenceInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::FenceInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__FenceInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__GlobalValue_Oval]"));
    core::BuiltInClass_sp classllvmo__GlobalValue_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__GlobalValue_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__GlobalValue_Oval,_lisp,llvmo::GlobalValue_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::GlobalValue_O>::id,llvmo::GlobalValue_O::static_classSymbol());
    llvmo::GlobalValue_O::___staticClass = classllvmo__GlobalValue_Oval;
#ifdef USE_MPS
    llvmo::GlobalValue_O::static_Kind = gctools::GCKind<llvmo::GlobalValue_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__GlobalValue_Oval,llvmo::GlobalValue_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::GlobalValue_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::GlobalValue_O>>::allocateClass();
        llvmo::GlobalValue_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::GlobalValue_O::static_className() % (void*)(llvmo::GlobalValue_O::static_allocator) );
    classllvmo__GlobalValue_Oval->setCreator(llvmo::GlobalValue_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::GlobalValue_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__GlobalValue_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__LandingPadInst_Oval]"));
    core::BuiltInClass_sp classllvmo__LandingPadInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__LandingPadInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__LandingPadInst_Oval,_lisp,llvmo::LandingPadInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::LandingPadInst_O>::id,llvmo::LandingPadInst_O::static_classSymbol());
    llvmo::LandingPadInst_O::___staticClass = classllvmo__LandingPadInst_Oval;
#ifdef USE_MPS
    llvmo::LandingPadInst_O::static_Kind = gctools::GCKind<llvmo::LandingPadInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__LandingPadInst_Oval,llvmo::LandingPadInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::LandingPadInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::LandingPadInst_O>>::allocateClass();
        llvmo::LandingPadInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::LandingPadInst_O::static_className() % (void*)(llvmo::LandingPadInst_O::static_allocator) );
    classllvmo__LandingPadInst_Oval->setCreator(llvmo::LandingPadInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::LandingPadInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__LandingPadInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__PHINode_Oval]"));
    core::BuiltInClass_sp classllvmo__PHINode_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__PHINode_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__PHINode_Oval,_lisp,llvmo::PHINode_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::PHINode_O>::id,llvmo::PHINode_O::static_classSymbol());
    llvmo::PHINode_O::___staticClass = classllvmo__PHINode_Oval;
#ifdef USE_MPS
    llvmo::PHINode_O::static_Kind = gctools::GCKind<llvmo::PHINode_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__PHINode_Oval,llvmo::PHINode_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::PHINode_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::PHINode_O>>::allocateClass();
        llvmo::PHINode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::PHINode_O::static_className() % (void*)(llvmo::PHINode_O::static_allocator) );
    classllvmo__PHINode_Oval->setCreator(llvmo::PHINode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::PHINode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__PHINode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__PointerType_Oval]"));
    core::BuiltInClass_sp classllvmo__PointerType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__PointerType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__PointerType_Oval,_lisp,llvmo::PointerType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::PointerType_O>::id,llvmo::PointerType_O::static_classSymbol());
    llvmo::PointerType_O::___staticClass = classllvmo__PointerType_Oval;
#ifdef USE_MPS
    llvmo::PointerType_O::static_Kind = gctools::GCKind<llvmo::PointerType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__PointerType_Oval,llvmo::PointerType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::PointerType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::PointerType_O>>::allocateClass();
        llvmo::PointerType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::PointerType_O::static_className() % (void*)(llvmo::PointerType_O::static_allocator) );
    classllvmo__PointerType_Oval->setCreator(llvmo::PointerType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::PointerType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__PointerType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__StoreInst_Oval]"));
    core::BuiltInClass_sp classllvmo__StoreInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__StoreInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__StoreInst_Oval,_lisp,llvmo::StoreInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::StoreInst_O>::id,llvmo::StoreInst_O::static_classSymbol());
    llvmo::StoreInst_O::___staticClass = classllvmo__StoreInst_Oval;
#ifdef USE_MPS
    llvmo::StoreInst_O::static_Kind = gctools::GCKind<llvmo::StoreInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__StoreInst_Oval,llvmo::StoreInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::StoreInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::StoreInst_O>>::allocateClass();
        llvmo::StoreInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::StoreInst_O::static_className() % (void*)(llvmo::StoreInst_O::static_allocator) );
    classllvmo__StoreInst_Oval->setCreator(llvmo::StoreInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::StoreInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__StoreInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__TargetLibraryInfo_Oval]"));
    core::BuiltInClass_sp classllvmo__TargetLibraryInfo_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__TargetLibraryInfo_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__TargetLibraryInfo_Oval,_lisp,llvmo::TargetLibraryInfo_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::TargetLibraryInfo_O>::id,llvmo::TargetLibraryInfo_O::static_classSymbol());
    llvmo::TargetLibraryInfo_O::___staticClass = classllvmo__TargetLibraryInfo_Oval;
#ifdef USE_MPS
    llvmo::TargetLibraryInfo_O::static_Kind = gctools::GCKind<llvmo::TargetLibraryInfo_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__TargetLibraryInfo_Oval,llvmo::TargetLibraryInfo_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::TargetLibraryInfo_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::TargetLibraryInfo_O>>::allocateClass();
        llvmo::TargetLibraryInfo_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::TargetLibraryInfo_O::static_className() % (void*)(llvmo::TargetLibraryInfo_O::static_allocator) );
    classllvmo__TargetLibraryInfo_Oval->setCreator(llvmo::TargetLibraryInfo_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::TargetLibraryInfo_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__TargetLibraryInfo_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__TerminatorInst_Oval]"));
    core::BuiltInClass_sp classllvmo__TerminatorInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__TerminatorInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__TerminatorInst_Oval,_lisp,llvmo::TerminatorInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::TerminatorInst_O>::id,llvmo::TerminatorInst_O::static_classSymbol());
    llvmo::TerminatorInst_O::___staticClass = classllvmo__TerminatorInst_Oval;
#ifdef USE_MPS
    llvmo::TerminatorInst_O::static_Kind = gctools::GCKind<llvmo::TerminatorInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__TerminatorInst_Oval,llvmo::TerminatorInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::TerminatorInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::TerminatorInst_O>>::allocateClass();
        llvmo::TerminatorInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::TerminatorInst_O::static_className() % (void*)(llvmo::TerminatorInst_O::static_allocator) );
    classllvmo__TerminatorInst_Oval->setCreator(llvmo::TerminatorInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::TerminatorInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__TerminatorInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__UnaryInstruction_Oval]"));
    core::BuiltInClass_sp classllvmo__UnaryInstruction_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__UnaryInstruction_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__UnaryInstruction_Oval,_lisp,llvmo::UnaryInstruction_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::UnaryInstruction_O>::id,llvmo::UnaryInstruction_O::static_classSymbol());
    llvmo::UnaryInstruction_O::___staticClass = classllvmo__UnaryInstruction_Oval;
#ifdef USE_MPS
    llvmo::UnaryInstruction_O::static_Kind = gctools::GCKind<llvmo::UnaryInstruction_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__UnaryInstruction_Oval,llvmo::UnaryInstruction_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::UnaryInstruction_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::UnaryInstruction_O>>::allocateClass();
        llvmo::UnaryInstruction_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::UnaryInstruction_O::static_className() % (void*)(llvmo::UnaryInstruction_O::static_allocator) );
    classllvmo__UnaryInstruction_Oval->setCreator(llvmo::UnaryInstruction_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::UnaryInstruction_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__UnaryInstruction_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__UndefValue_Oval]"));
    core::BuiltInClass_sp classllvmo__UndefValue_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__UndefValue_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__UndefValue_Oval,_lisp,llvmo::UndefValue_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::UndefValue_O>::id,llvmo::UndefValue_O::static_classSymbol());
    llvmo::UndefValue_O::___staticClass = classllvmo__UndefValue_Oval;
#ifdef USE_MPS
    llvmo::UndefValue_O::static_Kind = gctools::GCKind<llvmo::UndefValue_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__UndefValue_Oval,llvmo::UndefValue_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::UndefValue_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::UndefValue_O>>::allocateClass();
        llvmo::UndefValue_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::UndefValue_O::static_className() % (void*)(llvmo::UndefValue_O::static_allocator) );
    classllvmo__UndefValue_Oval->setCreator(llvmo::UndefValue_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::UndefValue_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__UndefValue_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__VectorType_Oval]"));
    core::BuiltInClass_sp classllvmo__VectorType_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__VectorType_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__VectorType_Oval,_lisp,llvmo::VectorType_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::VectorType_O>::id,llvmo::VectorType_O::static_classSymbol());
    llvmo::VectorType_O::___staticClass = classllvmo__VectorType_Oval;
#ifdef USE_MPS
    llvmo::VectorType_O::static_Kind = gctools::GCKind<llvmo::VectorType_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__VectorType_Oval,llvmo::VectorType_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::VectorType_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::VectorType_O>>::allocateClass();
        llvmo::VectorType_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::VectorType_O::static_className() % (void*)(llvmo::VectorType_O::static_allocator) );
    classllvmo__VectorType_Oval->setCreator(llvmo::VectorType_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::VectorType_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__VectorType_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_ClbindPkg
    LOG(BF("Creating class[classclbind__ClassRep_Oval]"));
    StandardClass_sp classclbind__ClassRep_Oval = StandardClass_O::createUncollectable();
    classclbind__ClassRep_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classclbind__ClassRep_Oval,_lisp,clbind::ClassRep_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<clbind::ClassRep_O>::id,clbind::ClassRep_O::static_classSymbol());
    clbind::ClassRep_O::___staticClass = classclbind__ClassRep_Oval;
#ifdef USE_MPS
    clbind::ClassRep_O::static_Kind = gctools::GCKind<clbind::ClassRep_O>::Kind;
#endif
    core::af_setf_findClass(classclbind__ClassRep_Oval,clbind::ClassRep_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<clbind::ClassRep_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<clbind::ClassRep_O>>::allocateClass();
        clbind::ClassRep_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% clbind::ClassRep_O::static_className() % (void*)(clbind::ClassRep_O::static_allocator) );
    classclbind__ClassRep_Oval->setCreator(clbind::ClassRep_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % clbind::ClassRep_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classclbind__ClassRep_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__FuncallableStandardClass_Oval]"));
    StandardClass_sp classcore__FuncallableStandardClass_Oval = StandardClass_O::createUncollectable();
    classcore__FuncallableStandardClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__FuncallableStandardClass_Oval,_lisp,core::FuncallableStandardClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::FuncallableStandardClass_O>::id,core::FuncallableStandardClass_O::static_classSymbol());
    core::FuncallableStandardClass_O::___staticClass = classcore__FuncallableStandardClass_Oval;
#ifdef USE_MPS
    core::FuncallableStandardClass_O::static_Kind = gctools::GCKind<core::FuncallableStandardClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__FuncallableStandardClass_Oval,core::FuncallableStandardClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::FuncallableStandardClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::FuncallableStandardClass_O>>::allocateClass();
        core::FuncallableStandardClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::FuncallableStandardClass_O::static_className() % (void*)(core::FuncallableStandardClass_O::static_allocator) );
    classcore__FuncallableStandardClass_Oval->setCreator(core::FuncallableStandardClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::FuncallableStandardClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__FuncallableStandardClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_CorePkg
    LOG(BF("Creating class[classcore__StandardClass_Oval]"));
    StandardClass_sp classcore__StandardClass_Oval = StandardClass_O::createUncollectable();
    classcore__StandardClass_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classcore__StandardClass_Oval,_lisp,core::StandardClass_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<core::StandardClass_O>::id,core::StandardClass_O::static_classSymbol());
    core::StandardClass_O::___staticClass = classcore__StandardClass_Oval;
#ifdef USE_MPS
    core::StandardClass_O::static_Kind = gctools::GCKind<core::StandardClass_O>::Kind;
#endif
    core::af_setf_findClass(classcore__StandardClass_Oval,core::StandardClass_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<core::StandardClass_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<core::StandardClass_O>>::allocateClass();
        core::StandardClass_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% core::StandardClass_O::static_className() % (void*)(core::StandardClass_O::static_allocator) );
    classcore__StandardClass_Oval->setCreator(core::StandardClass_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % core::StandardClass_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classcore__StandardClass_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__AllocaInst_Oval]"));
    core::BuiltInClass_sp classllvmo__AllocaInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__AllocaInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__AllocaInst_Oval,_lisp,llvmo::AllocaInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::AllocaInst_O>::id,llvmo::AllocaInst_O::static_classSymbol());
    llvmo::AllocaInst_O::___staticClass = classllvmo__AllocaInst_Oval;
#ifdef USE_MPS
    llvmo::AllocaInst_O::static_Kind = gctools::GCKind<llvmo::AllocaInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__AllocaInst_Oval,llvmo::AllocaInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::AllocaInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::AllocaInst_O>>::allocateClass();
        llvmo::AllocaInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::AllocaInst_O::static_className() % (void*)(llvmo::AllocaInst_O::static_allocator) );
    classllvmo__AllocaInst_Oval->setCreator(llvmo::AllocaInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::AllocaInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__AllocaInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__BranchInst_Oval]"));
    core::BuiltInClass_sp classllvmo__BranchInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__BranchInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__BranchInst_Oval,_lisp,llvmo::BranchInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::BranchInst_O>::id,llvmo::BranchInst_O::static_classSymbol());
    llvmo::BranchInst_O::___staticClass = classllvmo__BranchInst_Oval;
#ifdef USE_MPS
    llvmo::BranchInst_O::static_Kind = gctools::GCKind<llvmo::BranchInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__BranchInst_Oval,llvmo::BranchInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::BranchInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::BranchInst_O>>::allocateClass();
        llvmo::BranchInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::BranchInst_O::static_className() % (void*)(llvmo::BranchInst_O::static_allocator) );
    classllvmo__BranchInst_Oval->setCreator(llvmo::BranchInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::BranchInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__BranchInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ConstantDataArray_Oval]"));
    core::BuiltInClass_sp classllvmo__ConstantDataArray_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ConstantDataArray_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ConstantDataArray_Oval,_lisp,llvmo::ConstantDataArray_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ConstantDataArray_O>::id,llvmo::ConstantDataArray_O::static_classSymbol());
    llvmo::ConstantDataArray_O::___staticClass = classllvmo__ConstantDataArray_Oval;
#ifdef USE_MPS
    llvmo::ConstantDataArray_O::static_Kind = gctools::GCKind<llvmo::ConstantDataArray_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ConstantDataArray_Oval,llvmo::ConstantDataArray_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ConstantDataArray_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ConstantDataArray_O>>::allocateClass();
        llvmo::ConstantDataArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ConstantDataArray_O::static_className() % (void*)(llvmo::ConstantDataArray_O::static_allocator) );
    classllvmo__ConstantDataArray_Oval->setCreator(llvmo::ConstantDataArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ConstantDataArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ConstantDataArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__Function_Oval]"));
    core::BuiltInClass_sp classllvmo__Function_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__Function_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__Function_Oval,_lisp,llvmo::Function_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::Function_O>::id,llvmo::Function_O::static_classSymbol());
    llvmo::Function_O::___staticClass = classllvmo__Function_Oval;
#ifdef USE_MPS
    llvmo::Function_O::static_Kind = gctools::GCKind<llvmo::Function_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__Function_Oval,llvmo::Function_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::Function_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::Function_O>>::allocateClass();
        llvmo::Function_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::Function_O::static_className() % (void*)(llvmo::Function_O::static_allocator) );
    classllvmo__Function_Oval->setCreator(llvmo::Function_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::Function_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__Function_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__GlobalVariable_Oval]"));
    core::BuiltInClass_sp classllvmo__GlobalVariable_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__GlobalVariable_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__GlobalVariable_Oval,_lisp,llvmo::GlobalVariable_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::GlobalVariable_O>::id,llvmo::GlobalVariable_O::static_classSymbol());
    llvmo::GlobalVariable_O::___staticClass = classllvmo__GlobalVariable_Oval;
#ifdef USE_MPS
    llvmo::GlobalVariable_O::static_Kind = gctools::GCKind<llvmo::GlobalVariable_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__GlobalVariable_Oval,llvmo::GlobalVariable_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::GlobalVariable_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::GlobalVariable_O>>::allocateClass();
        llvmo::GlobalVariable_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::GlobalVariable_O::static_className() % (void*)(llvmo::GlobalVariable_O::static_allocator) );
    classllvmo__GlobalVariable_Oval->setCreator(llvmo::GlobalVariable_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::GlobalVariable_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__GlobalVariable_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__IndirectBrInst_Oval]"));
    core::BuiltInClass_sp classllvmo__IndirectBrInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__IndirectBrInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__IndirectBrInst_Oval,_lisp,llvmo::IndirectBrInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::IndirectBrInst_O>::id,llvmo::IndirectBrInst_O::static_classSymbol());
    llvmo::IndirectBrInst_O::___staticClass = classllvmo__IndirectBrInst_Oval;
#ifdef USE_MPS
    llvmo::IndirectBrInst_O::static_Kind = gctools::GCKind<llvmo::IndirectBrInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__IndirectBrInst_Oval,llvmo::IndirectBrInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::IndirectBrInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::IndirectBrInst_O>>::allocateClass();
        llvmo::IndirectBrInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::IndirectBrInst_O::static_className() % (void*)(llvmo::IndirectBrInst_O::static_allocator) );
    classllvmo__IndirectBrInst_Oval->setCreator(llvmo::IndirectBrInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::IndirectBrInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__IndirectBrInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__InvokeInst_Oval]"));
    core::BuiltInClass_sp classllvmo__InvokeInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__InvokeInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__InvokeInst_Oval,_lisp,llvmo::InvokeInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::InvokeInst_O>::id,llvmo::InvokeInst_O::static_classSymbol());
    llvmo::InvokeInst_O::___staticClass = classllvmo__InvokeInst_Oval;
#ifdef USE_MPS
    llvmo::InvokeInst_O::static_Kind = gctools::GCKind<llvmo::InvokeInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__InvokeInst_Oval,llvmo::InvokeInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::InvokeInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::InvokeInst_O>>::allocateClass();
        llvmo::InvokeInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::InvokeInst_O::static_className() % (void*)(llvmo::InvokeInst_O::static_allocator) );
    classllvmo__InvokeInst_Oval->setCreator(llvmo::InvokeInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::InvokeInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__InvokeInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__LoadInst_Oval]"));
    core::BuiltInClass_sp classllvmo__LoadInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__LoadInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__LoadInst_Oval,_lisp,llvmo::LoadInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::LoadInst_O>::id,llvmo::LoadInst_O::static_classSymbol());
    llvmo::LoadInst_O::___staticClass = classllvmo__LoadInst_Oval;
#ifdef USE_MPS
    llvmo::LoadInst_O::static_Kind = gctools::GCKind<llvmo::LoadInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__LoadInst_Oval,llvmo::LoadInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::LoadInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::LoadInst_O>>::allocateClass();
        llvmo::LoadInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::LoadInst_O::static_className() % (void*)(llvmo::LoadInst_O::static_allocator) );
    classllvmo__LoadInst_Oval->setCreator(llvmo::LoadInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::LoadInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__LoadInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ResumeInst_Oval]"));
    core::BuiltInClass_sp classllvmo__ResumeInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ResumeInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ResumeInst_Oval,_lisp,llvmo::ResumeInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ResumeInst_O>::id,llvmo::ResumeInst_O::static_classSymbol());
    llvmo::ResumeInst_O::___staticClass = classllvmo__ResumeInst_Oval;
#ifdef USE_MPS
    llvmo::ResumeInst_O::static_Kind = gctools::GCKind<llvmo::ResumeInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ResumeInst_Oval,llvmo::ResumeInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ResumeInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ResumeInst_O>>::allocateClass();
        llvmo::ResumeInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ResumeInst_O::static_className() % (void*)(llvmo::ResumeInst_O::static_allocator) );
    classllvmo__ResumeInst_Oval->setCreator(llvmo::ResumeInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ResumeInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ResumeInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__ReturnInst_Oval]"));
    core::BuiltInClass_sp classllvmo__ReturnInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__ReturnInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__ReturnInst_Oval,_lisp,llvmo::ReturnInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::ReturnInst_O>::id,llvmo::ReturnInst_O::static_classSymbol());
    llvmo::ReturnInst_O::___staticClass = classllvmo__ReturnInst_Oval;
#ifdef USE_MPS
    llvmo::ReturnInst_O::static_Kind = gctools::GCKind<llvmo::ReturnInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__ReturnInst_Oval,llvmo::ReturnInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::ReturnInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::ReturnInst_O>>::allocateClass();
        llvmo::ReturnInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::ReturnInst_O::static_className() % (void*)(llvmo::ReturnInst_O::static_allocator) );
    classllvmo__ReturnInst_Oval->setCreator(llvmo::ReturnInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::ReturnInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__ReturnInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__SwitchInst_Oval]"));
    core::BuiltInClass_sp classllvmo__SwitchInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__SwitchInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__SwitchInst_Oval,_lisp,llvmo::SwitchInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::SwitchInst_O>::id,llvmo::SwitchInst_O::static_classSymbol());
    llvmo::SwitchInst_O::___staticClass = classllvmo__SwitchInst_Oval;
#ifdef USE_MPS
    llvmo::SwitchInst_O::static_Kind = gctools::GCKind<llvmo::SwitchInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__SwitchInst_Oval,llvmo::SwitchInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::SwitchInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::SwitchInst_O>>::allocateClass();
        llvmo::SwitchInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::SwitchInst_O::static_className() % (void*)(llvmo::SwitchInst_O::static_allocator) );
    classllvmo__SwitchInst_Oval->setCreator(llvmo::SwitchInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::SwitchInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__SwitchInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__UnreachableInst_Oval]"));
    core::BuiltInClass_sp classllvmo__UnreachableInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__UnreachableInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__UnreachableInst_Oval,_lisp,llvmo::UnreachableInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::UnreachableInst_O>::id,llvmo::UnreachableInst_O::static_classSymbol());
    llvmo::UnreachableInst_O::___staticClass = classllvmo__UnreachableInst_Oval;
#ifdef USE_MPS
    llvmo::UnreachableInst_O::static_Kind = gctools::GCKind<llvmo::UnreachableInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__UnreachableInst_Oval,llvmo::UnreachableInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::UnreachableInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::UnreachableInst_O>>::allocateClass();
        llvmo::UnreachableInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::UnreachableInst_O::static_className() % (void*)(llvmo::UnreachableInst_O::static_allocator) );
    classllvmo__UnreachableInst_Oval->setCreator(llvmo::UnreachableInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::UnreachableInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__UnreachableInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#ifdef Use_LlvmoPkg
    LOG(BF("Creating class[classllvmo__VAArgInst_Oval]"));
    core::BuiltInClass_sp classllvmo__VAArgInst_Oval = core::BuiltInClass_O::createUncollectable();
    classllvmo__VAArgInst_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classllvmo__VAArgInst_Oval,_lisp,llvmo::VAArgInst_O::static_classSymbol());
    reg::lisp_associateClassIdWithClassSymbol(reg::registered_class<llvmo::VAArgInst_O>::id,llvmo::VAArgInst_O::static_classSymbol());
    llvmo::VAArgInst_O::___staticClass = classllvmo__VAArgInst_Oval;
#ifdef USE_MPS
    llvmo::VAArgInst_O::static_Kind = gctools::GCKind<llvmo::VAArgInst_O>::Kind;
#endif
    core::af_setf_findClass(classllvmo__VAArgInst_Oval,llvmo::VAArgInst_O::static_classSymbol(),true,_Nil<core::Environment_O>());
    {
        gctools::tagged_pointer<core::LispObjectCreator<llvmo::VAArgInst_O>> cb = gctools::ClassAllocator<core::LispObjectCreator<llvmo::VAArgInst_O>>::allocateClass();
        llvmo::VAArgInst_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% llvmo::VAArgInst_O::static_className() % (void*)(llvmo::VAArgInst_O::static_allocator) );
    classllvmo__VAArgInst_Oval->setCreator(llvmo::VAArgInst_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % llvmo::VAArgInst_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classllvmo__VAArgInst_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif
#endif // CREATE_CLASS
#undef CREATE_CLASS
#ifdef DUMP_INFO_CLASS // {
// Depends on nothing
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::ObjectSet_O @ %X") % classadapt__ObjectSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::ObjectSet_O::static_className() % adapt::ObjectSet_O::static_classSymbol() );
#endif
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::StringList_O @ %X") % classadapt__StringList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::StringList_O::static_className() % adapt::StringList_O::static_classSymbol() );
#endif
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::SymbolSet_O @ %X") % classadapt__SymbolSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::SymbolSet_O::static_className() % adapt::SymbolSet_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::T_O @ %X") % classcore__T_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::T_O::static_className() % core::T_O::static_classSymbol() );
#endif
#ifdef Use_AstToolingPkg

    LOG(BF("---    dump_info   --- className: asttooling::AstVisitor_O @ %X") % classasttooling__AstVisitor_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % asttooling::AstVisitor_O::static_className() % asttooling::AstVisitor_O::static_classSymbol() );
#endif
#ifdef Use_CffiPkg

    LOG(BF("---    dump_info   --- className: cffi::Pointer_O @ %X") % classcffi__Pointer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % cffi::Pointer_O::static_className() % cffi::Pointer_O::static_classSymbol() );
#endif
#ifdef Use_ClbindPkg

    LOG(BF("---    dump_info   --- className: clbind::ClassRegistry_O @ %X") % classclbind__ClassRegistry_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % clbind::ClassRegistry_O::static_className() % clbind::ClassRegistry_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Archive_O @ %X") % classcore__Archive_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Archive_O::static_className() % core::Archive_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Array_O @ %X") % classcore__Array_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Array_O::static_className() % core::Array_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Binder_O @ %X") % classcore__Binder_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Binder_O::static_className() % core::Binder_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::CandoException_O @ %X") % classcore__CandoException_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::CandoException_O::static_className() % core::CandoException_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Character_dummy_O @ %X") % classcore__Character_dummy_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Character_dummy_O::static_className() % core::Character_dummy_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Cons_O @ %X") % classcore__Cons_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Cons_O::static_className() % core::Cons_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::CxxObject_O @ %X") % classcore__CxxObject_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::CxxObject_O::static_className() % core::CxxObject_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::DirectoryEntry_O @ %X") % classcore__DirectoryEntry_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::DirectoryEntry_O::static_className() % core::DirectoryEntry_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Environment_O @ %X") % classcore__Environment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Environment_O::static_className() % core::Environment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ExternalObject_O @ %X") % classcore__ExternalObject_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ExternalObject_O::static_className() % core::ExternalObject_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::FileStatus_O @ %X") % classcore__FileStatus_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::FileStatus_O::static_className() % core::FileStatus_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Function_O @ %X") % classcore__Function_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Function_O::static_className() % core::Function_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::HashTable_O @ %X") % classcore__HashTable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::HashTable_O::static_className() % core::HashTable_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::IntArray_O @ %X") % classcore__IntArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::IntArray_O::static_className() % core::IntArray_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::InvocationHistoryFrameIterator_O @ %X") % classcore__InvocationHistoryFrameIterator_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::InvocationHistoryFrameIterator_O::static_className() % core::InvocationHistoryFrameIterator_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Iterator_O @ %X") % classcore__Iterator_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Iterator_O::static_className() % core::Iterator_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LambdaListHandler_O @ %X") % classcore__LambdaListHandler_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LambdaListHandler_O::static_className() % core::LambdaListHandler_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LightUserData_O @ %X") % classcore__LightUserData_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LightUserData_O::static_className() % core::LightUserData_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LoadTimeValues_O @ %X") % classcore__LoadTimeValues_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LoadTimeValues_O::static_className() % core::LoadTimeValues_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::MultiStringBuffer_O @ %X") % classcore__MultiStringBuffer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::MultiStringBuffer_O::static_className() % core::MultiStringBuffer_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Number_O @ %X") % classcore__Number_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Number_O::static_className() % core::Number_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Package_O @ %X") % classcore__Package_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Package_O::static_className() % core::Package_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Path_O @ %X") % classcore__Path_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Path_O::static_className() % core::Path_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Pathname_O @ %X") % classcore__Pathname_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Pathname_O::static_className() % core::Pathname_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Pointer_O @ %X") % classcore__Pointer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Pointer_O::static_className() % core::Pointer_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::PosixTimeDuration_O @ %X") % classcore__PosixTimeDuration_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::PosixTimeDuration_O::static_className() % core::PosixTimeDuration_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::PosixTime_O @ %X") % classcore__PosixTime_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::PosixTime_O::static_className() % core::PosixTime_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::RandomState_O @ %X") % classcore__RandomState_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::RandomState_O::static_className() % core::RandomState_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ReadTable_O @ %X") % classcore__ReadTable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ReadTable_O::static_className() % core::ReadTable_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Reader_O @ %X") % classcore__Reader_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Reader_O::static_className() % core::Reader_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Record_O @ %X") % classcore__Record_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Record_O::static_className() % core::Record_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::RegexMatch_O @ %X") % classcore__RegexMatch_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::RegexMatch_O::static_className() % core::RegexMatch_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Regex_O @ %X") % classcore__Regex_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Regex_O::static_className() % core::Regex_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SNode_O @ %X") % classcore__SNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SNode_O::static_className() % core::SNode_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SingleDispatchMethod_O @ %X") % classcore__SingleDispatchMethod_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SingleDispatchMethod_O::static_className() % core::SingleDispatchMethod_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SmallMap_O @ %X") % classcore__SmallMap_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SmallMap_O::static_className() % core::SmallMap_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SmallMultimap_O @ %X") % classcore__SmallMultimap_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SmallMultimap_O::static_className() % core::SmallMultimap_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SourceFileInfo_O @ %X") % classcore__SourceFileInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SourceFileInfo_O::static_className() % core::SourceFileInfo_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SourceManager_O @ %X") % classcore__SourceManager_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SourceManager_O::static_className() % core::SourceManager_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SourcePosInfo_O @ %X") % classcore__SourcePosInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SourcePosInfo_O::static_className() % core::SourcePosInfo_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StandardObject_O @ %X") % classcore__StandardObject_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StandardObject_O::static_className() % core::StandardObject_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Stream_O @ %X") % classcore__Stream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Stream_O::static_className() % core::Stream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StructureObject_O @ %X") % classcore__StructureObject_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StructureObject_O::static_className() % core::StructureObject_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SymbolToEnumConverter_O @ %X") % classcore__SymbolToEnumConverter_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SymbolToEnumConverter_O::static_className() % core::SymbolToEnumConverter_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Symbol_O @ %X") % classcore__Symbol_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Symbol_O::static_className() % core::Symbol_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::VaList_dummy_O @ %X") % classcore__VaList_dummy_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::VaList_dummy_O::static_className() % core::VaList_dummy_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::WeakHashTable_O @ %X") % classcore__WeakHashTable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::WeakHashTable_O::static_className() % core::WeakHashTable_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::WeakKeyMapping_O @ %X") % classcore__WeakKeyMapping_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::WeakKeyMapping_O::static_className() % core::WeakKeyMapping_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::WeakPointer_O @ %X") % classcore__WeakPointer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::WeakPointer_O::static_className() % core::WeakPointer_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::WrappedPointer_O @ %X") % classcore__WrappedPointer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::WrappedPointer_O::static_className() % core::WrappedPointer_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::AttributeSet_O @ %X") % classllvmo__AttributeSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::AttributeSet_O::static_className() % llvmo::AttributeSet_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Attribute_O @ %X") % classllvmo__Attribute_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Attribute_O::static_className() % llvmo::Attribute_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DebugInfo_O @ %X") % classllvmo__DebugInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DebugInfo_O::static_className() % llvmo::DebugInfo_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DebugLoc_O @ %X") % classllvmo__DebugLoc_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DebugLoc_O::static_className() % llvmo::DebugLoc_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::InsertPoint_O @ %X") % classllvmo__InsertPoint_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::InsertPoint_O::static_className() % llvmo::InsertPoint_O::static_classSymbol() );
#endif
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::IterateCons_O @ %X") % classadapt__IterateCons_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::IterateCons_O::static_className() % adapt::IterateCons_O::static_classSymbol() );
#endif
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::QDomNode_O @ %X") % classadapt__QDomNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::QDomNode_O::static_className() % adapt::QDomNode_O::static_classSymbol() );
#endif
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::StringSet_O @ %X") % classadapt__StringSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::StringSet_O::static_className() % adapt::StringSet_O::static_classSymbol() );
#endif
#ifdef Use_AdaptPkg

    LOG(BF("---    dump_info   --- className: adapt::SymbolList_O @ %X") % classadapt__SymbolList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::SymbolList_O::static_className() % adapt::SymbolList_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ActivationFrame_O @ %X") % classcore__ActivationFrame_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ActivationFrame_O::static_className() % core::ActivationFrame_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::AnsiStream_O @ %X") % classcore__AnsiStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::AnsiStream_O::static_className() % core::AnsiStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ArrayDisplaced_O @ %X") % classcore__ArrayDisplaced_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ArrayDisplaced_O::static_className() % core::ArrayDisplaced_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ArrayObjects_O @ %X") % classcore__ArrayObjects_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ArrayObjects_O::static_className() % core::ArrayObjects_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::BranchSNode_O @ %X") % classcore__BranchSNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::BranchSNode_O::static_className() % core::BranchSNode_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::CompiledFunction_O @ %X") % classcore__CompiledFunction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::CompiledFunction_O::static_className() % core::CompiledFunction_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Complex_O @ %X") % classcore__Complex_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Complex_O::static_className() % core::Complex_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::DirectoryIterator_O @ %X") % classcore__DirectoryIterator_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::DirectoryIterator_O::static_className() % core::DirectoryIterator_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ForeignData_O @ %X") % classcore__ForeignData_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ForeignData_O::static_className() % core::ForeignData_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::GlueEnvironment_O @ %X") % classcore__GlueEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::GlueEnvironment_O::static_className() % core::GlueEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::HashTableEq_O @ %X") % classcore__HashTableEq_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::HashTableEq_O::static_className() % core::HashTableEq_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::HashTableEql_O @ %X") % classcore__HashTableEql_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::HashTableEql_O::static_className() % core::HashTableEql_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::HashTableEqual_O @ %X") % classcore__HashTableEqual_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::HashTableEqual_O::static_className() % core::HashTableEqual_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::HashTableEqualp_O @ %X") % classcore__HashTableEqualp_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::HashTableEqualp_O::static_className() % core::HashTableEqualp_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Instance_O @ %X") % classcore__Instance_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Instance_O::static_className() % core::Instance_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LeafSNode_O @ %X") % classcore__LeafSNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LeafSNode_O::static_className() % core::LeafSNode_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LexicalEnvironment_O @ %X") % classcore__LexicalEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LexicalEnvironment_O::static_className() % core::LexicalEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LoadArchive_O @ %X") % classcore__LoadArchive_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LoadArchive_O::static_className() % core::LoadArchive_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LogicalPathname_O @ %X") % classcore__LogicalPathname_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LogicalPathname_O::static_className() % core::LogicalPathname_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Metaobject_O @ %X") % classcore__Metaobject_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Metaobject_O::static_className() % core::Metaobject_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Null_O @ %X") % classcore__Null_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Null_O::static_className() % core::Null_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Real_O @ %X") % classcore__Real_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Real_O::static_className() % core::Real_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::RecursiveDirectoryIterator_O @ %X") % classcore__RecursiveDirectoryIterator_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::RecursiveDirectoryIterator_O::static_className() % core::RecursiveDirectoryIterator_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SaveArchive_O @ %X") % classcore__SaveArchive_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SaveArchive_O::static_className() % core::SaveArchive_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SingleDispatchEffectiveMethodFunction_O @ %X") % classcore__SingleDispatchEffectiveMethodFunction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SingleDispatchEffectiveMethodFunction_O::static_className() % core::SingleDispatchEffectiveMethodFunction_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SingleDispatchGenericFunction_O @ %X") % classcore__SingleDispatchGenericFunction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SingleDispatchGenericFunction_O::static_className() % core::SingleDispatchGenericFunction_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SpecialForm_O @ %X") % classcore__SpecialForm_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SpecialForm_O::static_className() % core::SpecialForm_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::UserData_O @ %X") % classcore__UserData_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::UserData_O::static_className() % core::UserData_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Vector_O @ %X") % classcore__Vector_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Vector_O::static_className() % core::Vector_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::WeakKeyHashTable_O @ %X") % classcore__WeakKeyHashTable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::WeakKeyHashTable_O::static_className() % core::WeakKeyHashTable_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::BoundingBox_O @ %X") % classgeom__BoundingBox_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::BoundingBox_O::static_className() % geom::BoundingBox_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::Color_O @ %X") % classgeom__Color_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::Color_O::static_className() % geom::Color_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::OMatrix_O @ %X") % classgeom__OMatrix_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::OMatrix_O::static_className() % geom::OMatrix_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::OVector2_O @ %X") % classgeom__OVector2_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::OVector2_O::static_className() % geom::OVector2_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::OVector3_O @ %X") % classgeom__OVector3_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::OVector3_O::static_className() % geom::OVector3_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::APFloat_O @ %X") % classllvmo__APFloat_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::APFloat_O::static_className() % llvmo::APFloat_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::APInt_O @ %X") % classllvmo__APInt_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::APInt_O::static_className() % llvmo::APInt_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIArray_O @ %X") % classllvmo__DIArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIArray_O::static_className() % llvmo::DIArray_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIBasicType_O @ %X") % classllvmo__DIBasicType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIBasicType_O::static_className() % llvmo::DIBasicType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIBuilder_O @ %X") % classllvmo__DIBuilder_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIBuilder_O::static_className() % llvmo::DIBuilder_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DICompileUnit_O @ %X") % classllvmo__DICompileUnit_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DICompileUnit_O::static_className() % llvmo::DICompileUnit_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DICompositeType_O @ %X") % classllvmo__DICompositeType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DICompositeType_O::static_className() % llvmo::DICompositeType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIDerivedType_O @ %X") % classllvmo__DIDerivedType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIDerivedType_O::static_className() % llvmo::DIDerivedType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIDescriptor_O @ %X") % classllvmo__DIDescriptor_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIDescriptor_O::static_className() % llvmo::DIDescriptor_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIFile_O @ %X") % classllvmo__DIFile_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIFile_O::static_className() % llvmo::DIFile_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DILexicalBlock_O @ %X") % classllvmo__DILexicalBlock_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DILexicalBlock_O::static_className() % llvmo::DILexicalBlock_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIScope_O @ %X") % classllvmo__DIScope_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIScope_O::static_className() % llvmo::DIScope_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DISubprogram_O @ %X") % classllvmo__DISubprogram_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DISubprogram_O::static_className() % llvmo::DISubprogram_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DISubroutineType_O @ %X") % classllvmo__DISubroutineType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DISubroutineType_O::static_className() % llvmo::DISubroutineType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DITypeArray_O @ %X") % classllvmo__DITypeArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DITypeArray_O::static_className() % llvmo::DITypeArray_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DIType_O @ %X") % classllvmo__DIType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DIType_O::static_className() % llvmo::DIType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DataLayout_O @ %X") % classllvmo__DataLayout_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DataLayout_O::static_className() % llvmo::DataLayout_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::EngineBuilder_O @ %X") % classllvmo__EngineBuilder_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::EngineBuilder_O::static_className() % llvmo::EngineBuilder_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ExecutionEngine_O @ %X") % classllvmo__ExecutionEngine_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ExecutionEngine_O::static_className() % llvmo::ExecutionEngine_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::IRBuilderBase_O @ %X") % classllvmo__IRBuilderBase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::IRBuilderBase_O::static_className() % llvmo::IRBuilderBase_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::LLVMContext_O @ %X") % classllvmo__LLVMContext_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::LLVMContext_O::static_className() % llvmo::LLVMContext_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Linker_O @ %X") % classllvmo__Linker_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Linker_O::static_className() % llvmo::Linker_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::MCSubtargetInfo_O @ %X") % classllvmo__MCSubtargetInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::MCSubtargetInfo_O::static_className() % llvmo::MCSubtargetInfo_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Metadata_O @ %X") % classllvmo__Metadata_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Metadata_O::static_className() % llvmo::Metadata_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Module_O @ %X") % classllvmo__Module_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Module_O::static_className() % llvmo::Module_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::NamedMDNode_O @ %X") % classllvmo__NamedMDNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::NamedMDNode_O::static_className() % llvmo::NamedMDNode_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::PassManagerBase_O @ %X") % classllvmo__PassManagerBase_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::PassManagerBase_O::static_className() % llvmo::PassManagerBase_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::PassManagerBuilder_O @ %X") % classllvmo__PassManagerBuilder_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::PassManagerBuilder_O::static_className() % llvmo::PassManagerBuilder_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Pass_O @ %X") % classllvmo__Pass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Pass_O::static_className() % llvmo::Pass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::TargetMachine_O @ %X") % classllvmo__TargetMachine_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::TargetMachine_O::static_className() % llvmo::TargetMachine_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::TargetOptions_O @ %X") % classllvmo__TargetOptions_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::TargetOptions_O::static_className() % llvmo::TargetOptions_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Target_O @ %X") % classllvmo__Target_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Target_O::static_className() % llvmo::Target_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Triple_O @ %X") % classllvmo__Triple_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Triple_O::static_className() % llvmo::Triple_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Type_O @ %X") % classllvmo__Type_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Type_O::static_className() % llvmo::Type_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Value_O @ %X") % classllvmo__Value_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Value_O::static_className() % llvmo::Value_O::static_classSymbol() );
#endif
#ifdef Use_UnitsPkg

    LOG(BF("---    dump_info   --- className: units::Dimension_O @ %X") % classunits__Dimension_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::Dimension_O::static_className() % units::Dimension_O::static_classSymbol() );
#endif
#ifdef Use_UnitsPkg

    LOG(BF("---    dump_info   --- className: units::Quantity_O @ %X") % classunits__Quantity_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::Quantity_O::static_className() % units::Quantity_O::static_classSymbol() );
#endif
#ifdef Use_UnitsPkg

    LOG(BF("---    dump_info   --- className: units::Unit_O @ %X") % classunits__Unit_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::Unit_O::static_className() % units::Unit_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::BitVector_O @ %X") % classcore__BitVector_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::BitVector_O::static_className() % core::BitVector_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::BroadcastStream_O @ %X") % classcore__BroadcastStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::BroadcastStream_O::static_className() % core::BroadcastStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::CompileTimeEnvironment_O @ %X") % classcore__CompileTimeEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::CompileTimeEnvironment_O::static_className() % core::CompileTimeEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ConcatenatedStream_O @ %X") % classcore__ConcatenatedStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ConcatenatedStream_O::static_className() % core::ConcatenatedStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::EchoStream_O @ %X") % classcore__EchoStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::EchoStream_O::static_className() % core::EchoStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::FileStream_O @ %X") % classcore__FileStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::FileStream_O::static_className() % core::FileStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Float_O @ %X") % classcore__Float_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Float_O::static_className() % core::Float_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::FunctionFrame_O @ %X") % classcore__FunctionFrame_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::FunctionFrame_O::static_className() % core::FunctionFrame_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Rational_O @ %X") % classcore__Rational_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Rational_O::static_className() % core::Rational_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::RuntimeVisibleEnvironment_O @ %X") % classcore__RuntimeVisibleEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::RuntimeVisibleEnvironment_O::static_className() % core::RuntimeVisibleEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SexpLoadArchive_O @ %X") % classcore__SexpLoadArchive_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SexpLoadArchive_O::static_className() % core::SexpLoadArchive_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SexpSaveArchive_O @ %X") % classcore__SexpSaveArchive_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SexpSaveArchive_O::static_className() % core::SexpSaveArchive_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Specializer_O @ %X") % classcore__Specializer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Specializer_O::static_className() % core::Specializer_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StringStream_O @ %X") % classcore__StringStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StringStream_O::static_className() % core::StringStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::String_O @ %X") % classcore__String_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::String_O::static_className() % core::String_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SynonymStream_O @ %X") % classcore__SynonymStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SynonymStream_O::static_className() % core::SynonymStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::TagbodyFrame_O @ %X") % classcore__TagbodyFrame_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::TagbodyFrame_O::static_className() % core::TagbodyFrame_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::TwoWayStream_O @ %X") % classcore__TwoWayStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::TwoWayStream_O::static_className() % core::TwoWayStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ValueFrame_O @ %X") % classcore__ValueFrame_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ValueFrame_O::static_className() % core::ValueFrame_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::VectorObjects_O @ %X") % classcore__VectorObjects_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::VectorObjects_O::static_className() % core::VectorObjects_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::CoordinateArray_O @ %X") % classgeom__CoordinateArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::CoordinateArray_O::static_className() % geom::CoordinateArray_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Argument_O @ %X") % classllvmo__Argument_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Argument_O::static_className() % llvmo::Argument_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::BasicBlock_O @ %X") % classllvmo__BasicBlock_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::BasicBlock_O::static_className() % llvmo::BasicBlock_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::CompositeType_O @ %X") % classllvmo__CompositeType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::CompositeType_O::static_className() % llvmo::CompositeType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::FunctionPassManager_O @ %X") % classllvmo__FunctionPassManager_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::FunctionPassManager_O::static_className() % llvmo::FunctionPassManager_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::FunctionPass_O @ %X") % classllvmo__FunctionPass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::FunctionPass_O::static_className() % llvmo::FunctionPass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::FunctionType_O @ %X") % classllvmo__FunctionType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::FunctionType_O::static_className() % llvmo::FunctionType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::IRBuilder_O @ %X") % classllvmo__IRBuilder_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::IRBuilder_O::static_className() % llvmo::IRBuilder_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::IntegerType_O @ %X") % classllvmo__IntegerType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::IntegerType_O::static_className() % llvmo::IntegerType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::LLVMTargetMachine_O @ %X") % classllvmo__LLVMTargetMachine_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::LLVMTargetMachine_O::static_className() % llvmo::LLVMTargetMachine_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::MDNode_O @ %X") % classllvmo__MDNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::MDNode_O::static_className() % llvmo::MDNode_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::MDString_O @ %X") % classllvmo__MDString_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::MDString_O::static_className() % llvmo::MDString_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ModulePass_O @ %X") % classllvmo__ModulePass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ModulePass_O::static_className() % llvmo::ModulePass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::PassManager_O @ %X") % classllvmo__PassManager_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::PassManager_O::static_className() % llvmo::PassManager_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::TargetSubtargetInfo_O @ %X") % classllvmo__TargetSubtargetInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::TargetSubtargetInfo_O::static_className() % llvmo::TargetSubtargetInfo_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::User_O @ %X") % classllvmo__User_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::User_O::static_className() % llvmo::User_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ValueAsMetadata_O @ %X") % classllvmo__ValueAsMetadata_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ValueAsMetadata_O::static_className() % llvmo::ValueAsMetadata_O::static_classSymbol() );
#endif
#ifdef Use_UnitsPkg

    LOG(BF("---    dump_info   --- className: units::NamedUnit_O @ %X") % classunits__NamedUnit_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::NamedUnit_O::static_className() % units::NamedUnit_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::BlockEnvironment_O @ %X") % classcore__BlockEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::BlockEnvironment_O::static_className() % core::BlockEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::CatchEnvironment_O @ %X") % classcore__CatchEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::CatchEnvironment_O::static_className() % core::CatchEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Class_O @ %X") % classcore__Class_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Class_O::static_className() % core::Class_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::DoubleFloat_O @ %X") % classcore__DoubleFloat_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::DoubleFloat_O::static_className() % core::DoubleFloat_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::FunctionContainerEnvironment_O @ %X") % classcore__FunctionContainerEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::FunctionContainerEnvironment_O::static_className() % core::FunctionContainerEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::FunctionValueEnvironment_O @ %X") % classcore__FunctionValueEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::FunctionValueEnvironment_O::static_className() % core::FunctionValueEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::IOFileStream_O @ %X") % classcore__IOFileStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::IOFileStream_O::static_className() % core::IOFileStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::IOStreamStream_O @ %X") % classcore__IOStreamStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::IOStreamStream_O::static_className() % core::IOStreamStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Integer_O @ %X") % classcore__Integer_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Integer_O::static_className() % core::Integer_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::LongFloat_O @ %X") % classcore__LongFloat_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::LongFloat_O::static_className() % core::LongFloat_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::MacroletEnvironment_O @ %X") % classcore__MacroletEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::MacroletEnvironment_O::static_className() % core::MacroletEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Ratio_O @ %X") % classcore__Ratio_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Ratio_O::static_className() % core::Ratio_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ShortFloat_O @ %X") % classcore__ShortFloat_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ShortFloat_O::static_className() % core::ShortFloat_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SimpleBitVector_O @ %X") % classcore__SimpleBitVector_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SimpleBitVector_O::static_className() % core::SimpleBitVector_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SingleFloat_dummy_O @ %X") % classcore__SingleFloat_dummy_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SingleFloat_dummy_O::static_className() % core::SingleFloat_dummy_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StackValueEnvironment_O @ %X") % classcore__StackValueEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StackValueEnvironment_O::static_className() % core::StackValueEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Str_O @ %X") % classcore__Str_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Str_O::static_className() % core::Str_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StringInputStream_O @ %X") % classcore__StringInputStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StringInputStream_O::static_className() % core::StringInputStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StringOutputStream_O @ %X") % classcore__StringOutputStream_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StringOutputStream_O::static_className() % core::StringOutputStream_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::SymbolMacroletEnvironment_O @ %X") % classcore__SymbolMacroletEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::SymbolMacroletEnvironment_O::static_className() % core::SymbolMacroletEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::TagbodyEnvironment_O @ %X") % classcore__TagbodyEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::TagbodyEnvironment_O::static_className() % core::TagbodyEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::UnwindProtectEnvironment_O @ %X") % classcore__UnwindProtectEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::UnwindProtectEnvironment_O::static_className() % core::UnwindProtectEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ValueEnvironment_O @ %X") % classcore__ValueEnvironment_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ValueEnvironment_O::static_className() % core::ValueEnvironment_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::VectorObjectsWithFillPtr_O @ %X") % classcore__VectorObjectsWithFillPtr_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::VectorObjectsWithFillPtr_O::static_className() % core::VectorObjectsWithFillPtr_O::static_classSymbol() );
#endif
#ifdef Use_GeomPkg

    LOG(BF("---    dump_info   --- className: geom::CoordinateArrayWithHash_O @ %X") % classgeom__CoordinateArrayWithHash_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::CoordinateArrayWithHash_O::static_className() % geom::CoordinateArrayWithHash_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Constant_O @ %X") % classllvmo__Constant_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Constant_O::static_className() % llvmo::Constant_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ImmutablePass_O @ %X") % classllvmo__ImmutablePass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ImmutablePass_O::static_className() % llvmo::ImmutablePass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Instruction_O @ %X") % classllvmo__Instruction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Instruction_O::static_className() % llvmo::Instruction_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::SequentialType_O @ %X") % classllvmo__SequentialType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::SequentialType_O::static_className() % llvmo::SequentialType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::StructType_O @ %X") % classllvmo__StructType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::StructType_O::static_className() % llvmo::StructType_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Bignum_O @ %X") % classcore__Bignum_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Bignum_O::static_className() % core::Bignum_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::BuiltInClass_O @ %X") % classcore__BuiltInClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::BuiltInClass_O::static_className() % core::BuiltInClass_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::CxxClass_O @ %X") % classcore__CxxClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::CxxClass_O::static_className() % core::CxxClass_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::Fixnum_dummy_O @ %X") % classcore__Fixnum_dummy_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::Fixnum_dummy_O::static_className() % core::Fixnum_dummy_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::ForwardReferencedClass_O @ %X") % classcore__ForwardReferencedClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::ForwardReferencedClass_O::static_className() % core::ForwardReferencedClass_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StdClass_O @ %X") % classcore__StdClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StdClass_O::static_className() % core::StdClass_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StrWithFillPtr_O @ %X") % classcore__StrWithFillPtr_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StrWithFillPtr_O::static_className() % core::StrWithFillPtr_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StructureClass_O @ %X") % classcore__StructureClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StructureClass_O::static_className() % core::StructureClass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ArrayType_O @ %X") % classllvmo__ArrayType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ArrayType_O::static_className() % llvmo::ArrayType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::AtomicCmpXchgInst_O @ %X") % classllvmo__AtomicCmpXchgInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::AtomicCmpXchgInst_O::static_className() % llvmo::AtomicCmpXchgInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::AtomicRMWInst_O @ %X") % classllvmo__AtomicRMWInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::AtomicRMWInst_O::static_className() % llvmo::AtomicRMWInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::BlockAddress_O @ %X") % classllvmo__BlockAddress_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::BlockAddress_O::static_className() % llvmo::BlockAddress_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::CallInst_O @ %X") % classllvmo__CallInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::CallInst_O::static_className() % llvmo::CallInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantArray_O @ %X") % classllvmo__ConstantArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantArray_O::static_className() % llvmo::ConstantArray_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantDataSequential_O @ %X") % classllvmo__ConstantDataSequential_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantDataSequential_O::static_className() % llvmo::ConstantDataSequential_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantExpr_O @ %X") % classllvmo__ConstantExpr_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantExpr_O::static_className() % llvmo::ConstantExpr_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantFP_O @ %X") % classllvmo__ConstantFP_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantFP_O::static_className() % llvmo::ConstantFP_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantInt_O @ %X") % classllvmo__ConstantInt_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantInt_O::static_className() % llvmo::ConstantInt_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantPointerNull_O @ %X") % classllvmo__ConstantPointerNull_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantPointerNull_O::static_className() % llvmo::ConstantPointerNull_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantStruct_O @ %X") % classllvmo__ConstantStruct_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantStruct_O::static_className() % llvmo::ConstantStruct_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::DataLayoutPass_O @ %X") % classllvmo__DataLayoutPass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::DataLayoutPass_O::static_className() % llvmo::DataLayoutPass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::FenceInst_O @ %X") % classllvmo__FenceInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::FenceInst_O::static_className() % llvmo::FenceInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::GlobalValue_O @ %X") % classllvmo__GlobalValue_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::GlobalValue_O::static_className() % llvmo::GlobalValue_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::LandingPadInst_O @ %X") % classllvmo__LandingPadInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::LandingPadInst_O::static_className() % llvmo::LandingPadInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::PHINode_O @ %X") % classllvmo__PHINode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::PHINode_O::static_className() % llvmo::PHINode_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::PointerType_O @ %X") % classllvmo__PointerType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::PointerType_O::static_className() % llvmo::PointerType_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::StoreInst_O @ %X") % classllvmo__StoreInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::StoreInst_O::static_className() % llvmo::StoreInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::TargetLibraryInfo_O @ %X") % classllvmo__TargetLibraryInfo_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::TargetLibraryInfo_O::static_className() % llvmo::TargetLibraryInfo_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::TerminatorInst_O @ %X") % classllvmo__TerminatorInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::TerminatorInst_O::static_className() % llvmo::TerminatorInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::UnaryInstruction_O @ %X") % classllvmo__UnaryInstruction_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::UnaryInstruction_O::static_className() % llvmo::UnaryInstruction_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::UndefValue_O @ %X") % classllvmo__UndefValue_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::UndefValue_O::static_className() % llvmo::UndefValue_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::VectorType_O @ %X") % classllvmo__VectorType_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::VectorType_O::static_className() % llvmo::VectorType_O::static_classSymbol() );
#endif
#ifdef Use_ClbindPkg

    LOG(BF("---    dump_info   --- className: clbind::ClassRep_O @ %X") % classclbind__ClassRep_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % clbind::ClassRep_O::static_className() % clbind::ClassRep_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::FuncallableStandardClass_O @ %X") % classcore__FuncallableStandardClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::FuncallableStandardClass_O::static_className() % core::FuncallableStandardClass_O::static_classSymbol() );
#endif
#ifdef Use_CorePkg

    LOG(BF("---    dump_info   --- className: core::StandardClass_O @ %X") % classcore__StandardClass_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % core::StandardClass_O::static_className() % core::StandardClass_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::AllocaInst_O @ %X") % classllvmo__AllocaInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::AllocaInst_O::static_className() % llvmo::AllocaInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::BranchInst_O @ %X") % classllvmo__BranchInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::BranchInst_O::static_className() % llvmo::BranchInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ConstantDataArray_O @ %X") % classllvmo__ConstantDataArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ConstantDataArray_O::static_className() % llvmo::ConstantDataArray_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::Function_O @ %X") % classllvmo__Function_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::Function_O::static_className() % llvmo::Function_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::GlobalVariable_O @ %X") % classllvmo__GlobalVariable_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::GlobalVariable_O::static_className() % llvmo::GlobalVariable_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::IndirectBrInst_O @ %X") % classllvmo__IndirectBrInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::IndirectBrInst_O::static_className() % llvmo::IndirectBrInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::InvokeInst_O @ %X") % classllvmo__InvokeInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::InvokeInst_O::static_className() % llvmo::InvokeInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::LoadInst_O @ %X") % classllvmo__LoadInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::LoadInst_O::static_className() % llvmo::LoadInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ResumeInst_O @ %X") % classllvmo__ResumeInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ResumeInst_O::static_className() % llvmo::ResumeInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::ReturnInst_O @ %X") % classllvmo__ReturnInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::ReturnInst_O::static_className() % llvmo::ReturnInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::SwitchInst_O @ %X") % classllvmo__SwitchInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::SwitchInst_O::static_className() % llvmo::SwitchInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::UnreachableInst_O @ %X") % classllvmo__UnreachableInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::UnreachableInst_O::static_className() % llvmo::UnreachableInst_O::static_classSymbol() );
#endif
#ifdef Use_LlvmoPkg

    LOG(BF("---    dump_info   --- className: llvmo::VAArgInst_O @ %X") % classllvmo__VAArgInst_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % llvmo::VAArgInst_O::static_className() % llvmo::VAArgInst_O::static_classSymbol() );
#endif
#endif // } DUMP_INFO_CLASS
#undef DUMP_INFO_CLASS
#if defined(DEFINE_BASE_CLASSES) || defined(ALL_STAGES) // {
// Depends on nothing
#ifdef Use_AdaptPkg
classadapt__ObjectSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(adapt::T_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg
classadapt__StringList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(adapt::T_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg
classadapt__SymbolSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(adapt::T_O::static_classSymbol());
#endif
#ifdef Use_AstToolingPkg
classasttooling__AstVisitor_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CffiPkg
classcffi__Pointer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_ClbindPkg
classclbind__ClassRegistry_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Archive_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Array_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Binder_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__CandoException_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Character_dummy_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Cons_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__CxxObject_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__DirectoryEntry_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Environment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ExternalObject_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__FileStatus_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Function_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__HashTable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__IntArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__InvocationHistoryFrameIterator_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Iterator_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LambdaListHandler_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LightUserData_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LoadTimeValues_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__MultiStringBuffer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Number_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Package_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Path_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Pathname_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Pointer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__PosixTimeDuration_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__PosixTime_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__RandomState_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ReadTable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Reader_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Record_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__RegexMatch_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Regex_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SingleDispatchMethod_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SmallMap_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SmallMultimap_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SourceFileInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SourceManager_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SourcePosInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StandardObject_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Stream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StructureObject_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SymbolToEnumConverter_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Symbol_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__VaList_dummy_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__WeakHashTable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__WeakKeyMapping_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__WeakPointer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__WrappedPointer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__AttributeSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Attribute_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DebugInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DebugLoc_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__InsertPoint_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg
classadapt__IterateCons_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg
classadapt__QDomNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg
classadapt__StringSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg
classadapt__SymbolList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ActivationFrame_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Environment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__AnsiStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Stream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ArrayDisplaced_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Array_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ArrayObjects_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Array_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__BranchSNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::SNode_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__CompiledFunction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Function_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Complex_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Number_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__DirectoryIterator_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ForeignData_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__GlueEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Environment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__HashTableEq_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::HashTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__HashTableEql_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::HashTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__HashTableEqual_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::HashTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__HashTableEqualp_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::HashTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Instance_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Function_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LeafSNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::SNode_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LexicalEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Environment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LoadArchive_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Archive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LogicalPathname_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Pathname_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Metaobject_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::StandardObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Null_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Symbol_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Real_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Number_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__RecursiveDirectoryIterator_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SaveArchive_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Archive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SingleDispatchEffectiveMethodFunction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Function_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SingleDispatchGenericFunction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Function_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SpecialForm_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Function_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__UserData_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::LightUserData_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Vector_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Array_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__WeakKeyHashTable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::WeakHashTable_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__BoundingBox_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__Color_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__OMatrix_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__OVector2_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__OVector3_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__APFloat_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__APInt_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIBasicType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIBuilder_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DICompileUnit_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DICompositeType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIDerivedType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIDescriptor_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIFile_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DILexicalBlock_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIScope_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DISubprogram_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DISubroutineType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DITypeArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DIType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DataLayout_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__EngineBuilder_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ExecutionEngine_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__IRBuilderBase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__LLVMContext_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Linker_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__MCSubtargetInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Metadata_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Module_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__NamedMDNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__PassManagerBase_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__PassManagerBuilder_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Pass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__TargetMachine_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__TargetOptions_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Target_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Triple_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Type_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Value_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg
classunits__Dimension_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg
classunits__Quantity_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg
classunits__Unit_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__BitVector_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__BroadcastStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__CompileTimeEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::LexicalEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ConcatenatedStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__EchoStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__FileStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Float_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Real_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__FunctionFrame_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ActivationFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Rational_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Real_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__RuntimeVisibleEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::LexicalEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SexpLoadArchive_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::LoadArchive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SexpSaveArchive_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::SaveArchive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Specializer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Metaobject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StringStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__String_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SynonymStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__TagbodyFrame_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ActivationFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__TwoWayStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ValueFrame_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::ActivationFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__VectorObjects_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__CoordinateArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Argument_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Value_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__BasicBlock_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Value_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__CompositeType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Type_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__FunctionPassManager_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::PassManagerBase_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__FunctionPass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Pass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__FunctionType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Type_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__IRBuilder_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::IRBuilderBase_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__IntegerType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Type_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__LLVMTargetMachine_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TargetMachine_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__MDNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Metadata_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__MDString_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Metadata_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ModulePass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Pass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__PassManager_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::PassManagerBase_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__TargetSubtargetInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::MCSubtargetInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__User_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Value_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ValueAsMetadata_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Metadata_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg
classunits__NamedUnit_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(units::Unit_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__BlockEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__CatchEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Class_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Specializer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__DoubleFloat_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Float_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__FunctionContainerEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__FunctionValueEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::RuntimeVisibleEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__IOFileStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::FileStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__IOStreamStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::FileStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Integer_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Rational_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__LongFloat_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Float_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__MacroletEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Ratio_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Rational_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ShortFloat_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Float_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SimpleBitVector_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::BitVector_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SingleFloat_dummy_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Float_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StackValueEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Str_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::String_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StringInputStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::StringStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StringOutputStream_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::StringStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__SymbolMacroletEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__TagbodyEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::RuntimeVisibleEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__UnwindProtectEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ValueEnvironment_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::RuntimeVisibleEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__VectorObjectsWithFillPtr_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::VectorObjects_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg
classgeom__CoordinateArrayWithHash_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(geom::CoordinateArray_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Constant_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::User_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ImmutablePass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::ModulePass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Instruction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::User_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__SequentialType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::CompositeType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__StructType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::CompositeType_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Bignum_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Integer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__BuiltInClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Class_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__CxxClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Class_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__Fixnum_dummy_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Integer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__ForwardReferencedClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Class_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StdClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Class_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StrWithFillPtr_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Str_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StructureClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Class_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ArrayType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::SequentialType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__AtomicCmpXchgInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__AtomicRMWInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__BlockAddress_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__CallInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantDataSequential_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantExpr_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantFP_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantInt_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantPointerNull_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantStruct_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__DataLayoutPass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::ImmutablePass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__FenceInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__GlobalValue_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__LandingPadInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__PHINode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__PointerType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::SequentialType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__StoreInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__TargetLibraryInfo_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::ImmutablePass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__TerminatorInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__UnaryInstruction_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__UndefValue_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__VectorType_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::SequentialType_O::static_classSymbol());
#endif
#ifdef Use_ClbindPkg
classclbind__ClassRep_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::BuiltInClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__FuncallableStandardClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::StdClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg
classcore__StandardClass_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::StdClass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__AllocaInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::UnaryInstruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__BranchInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ConstantDataArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::ConstantDataSequential_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__Function_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::GlobalValue_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__GlobalVariable_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::GlobalValue_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__IndirectBrInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__InvokeInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__LoadInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::UnaryInstruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ResumeInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__ReturnInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__SwitchInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__UnreachableInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg
classllvmo__VAArgInst_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(llvmo::UnaryInstruction_O::static_classSymbol());
#endif
#endif // } DEFINE_BASE_CLASSES
#undef DEFINE_BASE_CLASSES
#if defined(DEFINE_CLASS_NAMES) || defined(ALL_STAGES) // {
// Depends on nothing
#ifdef Use_AdaptPkg

    classadapt__ObjectSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::ObjectSet_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg

    classadapt__StringList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::StringList_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg

    classadapt__SymbolSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::SymbolSet_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__T_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::T_O::static_classSymbol());
#endif
#ifdef Use_AstToolingPkg

    classasttooling__AstVisitor_Oval->__setupStage3NameAndCalculateClassPrecedenceList(asttooling::AstVisitor_O::static_classSymbol());
#endif
#ifdef Use_CffiPkg

    classcffi__Pointer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(cffi::Pointer_O::static_classSymbol());
#endif
#ifdef Use_ClbindPkg

    classclbind__ClassRegistry_Oval->__setupStage3NameAndCalculateClassPrecedenceList(clbind::ClassRegistry_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Archive_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Archive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Array_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Array_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Binder_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Binder_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__CandoException_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::CandoException_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Character_dummy_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Character_dummy_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Cons_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Cons_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__CxxObject_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::CxxObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__DirectoryEntry_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::DirectoryEntry_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Environment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Environment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ExternalObject_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ExternalObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__FileStatus_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::FileStatus_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Function_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Function_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__HashTable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::HashTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__IntArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::IntArray_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__InvocationHistoryFrameIterator_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::InvocationHistoryFrameIterator_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Iterator_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LambdaListHandler_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LambdaListHandler_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LightUserData_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LightUserData_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LoadTimeValues_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LoadTimeValues_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__MultiStringBuffer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::MultiStringBuffer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Number_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Number_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Package_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Package_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Path_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Path_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Pathname_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Pathname_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Pointer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Pointer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__PosixTimeDuration_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::PosixTimeDuration_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__PosixTime_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::PosixTime_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__RandomState_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::RandomState_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ReadTable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ReadTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Reader_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Reader_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Record_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Record_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__RegexMatch_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::RegexMatch_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Regex_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Regex_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SNode_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SingleDispatchMethod_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SingleDispatchMethod_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SmallMap_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SmallMap_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SmallMultimap_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SmallMultimap_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SourceFileInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SourceFileInfo_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SourceManager_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SourceManager_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SourcePosInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SourcePosInfo_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StandardObject_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StandardObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Stream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Stream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StructureObject_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StructureObject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SymbolToEnumConverter_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SymbolToEnumConverter_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Symbol_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Symbol_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__VaList_dummy_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::VaList_dummy_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__WeakHashTable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::WeakHashTable_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__WeakKeyMapping_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::WeakKeyMapping_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__WeakPointer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::WeakPointer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__WrappedPointer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::WrappedPointer_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__AttributeSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::AttributeSet_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Attribute_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Attribute_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DebugInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DebugInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DebugLoc_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DebugLoc_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__InsertPoint_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::InsertPoint_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg

    classadapt__IterateCons_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::IterateCons_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg

    classadapt__QDomNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::QDomNode_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg

    classadapt__StringSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::StringSet_O::static_classSymbol());
#endif
#ifdef Use_AdaptPkg

    classadapt__SymbolList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::SymbolList_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ActivationFrame_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ActivationFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__AnsiStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::AnsiStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ArrayDisplaced_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ArrayDisplaced_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ArrayObjects_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ArrayObjects_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__BranchSNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::BranchSNode_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__CompiledFunction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::CompiledFunction_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Complex_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Complex_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__DirectoryIterator_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::DirectoryIterator_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ForeignData_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ForeignData_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__GlueEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::GlueEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__HashTableEq_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::HashTableEq_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__HashTableEql_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::HashTableEql_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__HashTableEqual_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::HashTableEqual_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__HashTableEqualp_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::HashTableEqualp_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Instance_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Instance_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LeafSNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LeafSNode_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LexicalEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LexicalEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LoadArchive_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LoadArchive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LogicalPathname_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LogicalPathname_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Metaobject_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Metaobject_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Null_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Null_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Real_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Real_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__RecursiveDirectoryIterator_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::RecursiveDirectoryIterator_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SaveArchive_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SaveArchive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SingleDispatchEffectiveMethodFunction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SingleDispatchEffectiveMethodFunction_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SingleDispatchGenericFunction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SingleDispatchGenericFunction_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SpecialForm_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SpecialForm_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__UserData_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::UserData_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Vector_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__WeakKeyHashTable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::WeakKeyHashTable_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__BoundingBox_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::BoundingBox_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__Color_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::Color_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__OMatrix_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::OMatrix_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__OVector2_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::OVector2_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__OVector3_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::OVector3_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__APFloat_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::APFloat_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__APInt_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::APInt_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIArray_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIBasicType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIBasicType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIBuilder_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIBuilder_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DICompileUnit_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DICompileUnit_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DICompositeType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DICompositeType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIDerivedType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIDerivedType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIDescriptor_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIDescriptor_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIFile_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIFile_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DILexicalBlock_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DILexicalBlock_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIScope_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIScope_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DISubprogram_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DISubprogram_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DISubroutineType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DISubroutineType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DITypeArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DITypeArray_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DIType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DIType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DataLayout_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DataLayout_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__EngineBuilder_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::EngineBuilder_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ExecutionEngine_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ExecutionEngine_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__IRBuilderBase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::IRBuilderBase_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__LLVMContext_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::LLVMContext_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Linker_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Linker_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__MCSubtargetInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::MCSubtargetInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Metadata_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Metadata_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Module_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Module_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__NamedMDNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::NamedMDNode_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__PassManagerBase_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::PassManagerBase_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__PassManagerBuilder_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::PassManagerBuilder_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Pass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Pass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__TargetMachine_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::TargetMachine_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__TargetOptions_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::TargetOptions_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Target_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Target_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Triple_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Triple_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Type_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Type_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Value_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Value_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg

    classunits__Dimension_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::Dimension_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg

    classunits__Quantity_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::Quantity_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg

    classunits__Unit_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::Unit_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__BitVector_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::BitVector_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__BroadcastStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::BroadcastStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__CompileTimeEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::CompileTimeEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ConcatenatedStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ConcatenatedStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__EchoStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::EchoStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__FileStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::FileStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Float_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Float_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__FunctionFrame_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::FunctionFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Rational_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Rational_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__RuntimeVisibleEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::RuntimeVisibleEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SexpLoadArchive_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SexpLoadArchive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SexpSaveArchive_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SexpSaveArchive_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Specializer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Specializer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StringStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StringStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__String_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::String_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SynonymStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SynonymStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__TagbodyFrame_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::TagbodyFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__TwoWayStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::TwoWayStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ValueFrame_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ValueFrame_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__VectorObjects_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::VectorObjects_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__CoordinateArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::CoordinateArray_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Argument_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Argument_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__BasicBlock_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::BasicBlock_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__CompositeType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::CompositeType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__FunctionPassManager_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::FunctionPassManager_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__FunctionPass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::FunctionPass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__FunctionType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::FunctionType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__IRBuilder_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::IRBuilder_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__IntegerType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::IntegerType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__LLVMTargetMachine_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::LLVMTargetMachine_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__MDNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::MDNode_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__MDString_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::MDString_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ModulePass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ModulePass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__PassManager_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::PassManager_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__TargetSubtargetInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::TargetSubtargetInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__User_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::User_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ValueAsMetadata_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ValueAsMetadata_O::static_classSymbol());
#endif
#ifdef Use_UnitsPkg

    classunits__NamedUnit_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::NamedUnit_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__BlockEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::BlockEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__CatchEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::CatchEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Class_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Class_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__DoubleFloat_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::DoubleFloat_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__FunctionContainerEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::FunctionContainerEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__FunctionValueEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::FunctionValueEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__IOFileStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::IOFileStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__IOStreamStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::IOStreamStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Integer_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Integer_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__LongFloat_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::LongFloat_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__MacroletEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::MacroletEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Ratio_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Ratio_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ShortFloat_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ShortFloat_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SimpleBitVector_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SimpleBitVector_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SingleFloat_dummy_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SingleFloat_dummy_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StackValueEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StackValueEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Str_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Str_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StringInputStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StringInputStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StringOutputStream_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StringOutputStream_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__SymbolMacroletEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::SymbolMacroletEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__TagbodyEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::TagbodyEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__UnwindProtectEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::UnwindProtectEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ValueEnvironment_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ValueEnvironment_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__VectorObjectsWithFillPtr_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::VectorObjectsWithFillPtr_O::static_classSymbol());
#endif
#ifdef Use_GeomPkg

    classgeom__CoordinateArrayWithHash_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::CoordinateArrayWithHash_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Constant_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Constant_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ImmutablePass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ImmutablePass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Instruction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Instruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__SequentialType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::SequentialType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__StructType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::StructType_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Bignum_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Bignum_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__BuiltInClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::BuiltInClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__CxxClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::CxxClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__Fixnum_dummy_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::Fixnum_dummy_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__ForwardReferencedClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::ForwardReferencedClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StdClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StdClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StrWithFillPtr_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StrWithFillPtr_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StructureClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StructureClass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ArrayType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ArrayType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__AtomicCmpXchgInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::AtomicCmpXchgInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__AtomicRMWInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::AtomicRMWInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__BlockAddress_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::BlockAddress_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__CallInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::CallInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantArray_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantDataSequential_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantDataSequential_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantExpr_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantExpr_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantFP_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantFP_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantInt_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantInt_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantPointerNull_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantPointerNull_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantStruct_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantStruct_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__DataLayoutPass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::DataLayoutPass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__FenceInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::FenceInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__GlobalValue_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::GlobalValue_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__LandingPadInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::LandingPadInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__PHINode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::PHINode_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__PointerType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::PointerType_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__StoreInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::StoreInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__TargetLibraryInfo_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::TargetLibraryInfo_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__TerminatorInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::TerminatorInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__UnaryInstruction_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::UnaryInstruction_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__UndefValue_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::UndefValue_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__VectorType_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::VectorType_O::static_classSymbol());
#endif
#ifdef Use_ClbindPkg

    classclbind__ClassRep_Oval->__setupStage3NameAndCalculateClassPrecedenceList(clbind::ClassRep_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__FuncallableStandardClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::FuncallableStandardClass_O::static_classSymbol());
#endif
#ifdef Use_CorePkg

    classcore__StandardClass_Oval->__setupStage3NameAndCalculateClassPrecedenceList(core::StandardClass_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__AllocaInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::AllocaInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__BranchInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::BranchInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ConstantDataArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ConstantDataArray_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__Function_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::Function_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__GlobalVariable_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::GlobalVariable_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__IndirectBrInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::IndirectBrInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__InvokeInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::InvokeInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__LoadInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::LoadInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ResumeInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ResumeInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__ReturnInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::ReturnInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__SwitchInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::SwitchInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__UnreachableInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::UnreachableInst_O::static_classSymbol());
#endif
#ifdef Use_LlvmoPkg

    classllvmo__VAArgInst_Oval->__setupStage3NameAndCalculateClassPrecedenceList(llvmo::VAArgInst_O::static_classSymbol());
#endif
#endif // } DEFINE_CLASS_NAMES
#undef DEFINE_CLASS_NAMES
#if defined(EXPOSE_TO_CANDO) || defined(ALL_STAGES)
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__ObjectSet_O(core::Lisp_sp); // base(s): set(['adapt::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__ObjectSet_O");
    adapt::Register_adapt__ObjectSet_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__StringList_O(core::Lisp_sp); // base(s): set(['adapt::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__StringList_O");
    adapt::Register_adapt__StringList_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__SymbolSet_O(core::Lisp_sp); // base(s): set(['adapt::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__SymbolSet_O");
    adapt::Register_adapt__SymbolSet_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__T_O(core::Lisp_sp); // base(s): set(['core::_RootDummyClass'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__T_O");
    core::Register_core__T_O(_lisp); // base(s): set(['core::_RootDummyClass'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_AstToolingPkg
#ifdef EXTERN_REGISTER
extern void Register_asttooling__AstVisitor_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_asttooling__AstVisitor_O");
    asttooling::Register_asttooling__AstVisitor_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AstToolingPkg
#ifdef Use_CffiPkg
#ifdef EXTERN_REGISTER
extern void Register_cffi__Pointer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_cffi__Pointer_O");
    cffi::Register_cffi__Pointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CffiPkg
#ifdef Use_ClbindPkg
#ifdef EXTERN_REGISTER
extern void Register_clbind__ClassRegistry_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_clbind__ClassRegistry_O");
    clbind::Register_clbind__ClassRegistry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClbindPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Archive_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Archive_O");
    core::Register_core__Archive_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Array_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Array_O");
    core::Register_core__Array_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Binder_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Binder_O");
    core::Register_core__Binder_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__CandoException_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__CandoException_O");
    core::Register_core__CandoException_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Character_dummy_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Character_dummy_O");
    core::Register_core__Character_dummy_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Cons_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Cons_O");
    core::Register_core__Cons_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__CxxObject_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__CxxObject_O");
    core::Register_core__CxxObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__DirectoryEntry_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__DirectoryEntry_O");
    core::Register_core__DirectoryEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Environment_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Environment_O");
    core::Register_core__Environment_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ExternalObject_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ExternalObject_O");
    core::Register_core__ExternalObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__FileStatus_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__FileStatus_O");
    core::Register_core__FileStatus_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Function_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Function_O");
    core::Register_core__Function_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__HashTable_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__HashTable_O");
    core::Register_core__HashTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__IntArray_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__IntArray_O");
    core::Register_core__IntArray_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__InvocationHistoryFrameIterator_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__InvocationHistoryFrameIterator_O");
    core::Register_core__InvocationHistoryFrameIterator_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Iterator_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Iterator_O");
    core::Register_core__Iterator_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__LambdaListHandler_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LambdaListHandler_O");
    core::Register_core__LambdaListHandler_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__LightUserData_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LightUserData_O");
    core::Register_core__LightUserData_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__LoadTimeValues_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LoadTimeValues_O");
    core::Register_core__LoadTimeValues_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__MultiStringBuffer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__MultiStringBuffer_O");
    core::Register_core__MultiStringBuffer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Number_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Number_O");
    core::Register_core__Number_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Package_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Package_O");
    core::Register_core__Package_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Path_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Path_O");
    core::Register_core__Path_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Pathname_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Pathname_O");
    core::Register_core__Pathname_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Pointer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Pointer_O");
    core::Register_core__Pointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__PosixTimeDuration_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__PosixTimeDuration_O");
    core::Register_core__PosixTimeDuration_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__PosixTime_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__PosixTime_O");
    core::Register_core__PosixTime_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__RandomState_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__RandomState_O");
    core::Register_core__RandomState_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__ReadTable_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ReadTable_O");
    core::Register_core__ReadTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Reader_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Reader_O");
    core::Register_core__Reader_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Record_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Record_O");
    core::Register_core__Record_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__RegexMatch_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__RegexMatch_O");
    core::Register_core__RegexMatch_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Regex_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Regex_O");
    core::Register_core__Regex_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SNode_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SNode_O");
    core::Register_core__SNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SingleDispatchMethod_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SingleDispatchMethod_O");
    core::Register_core__SingleDispatchMethod_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SmallMap_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SmallMap_O");
    core::Register_core__SmallMap_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SmallMultimap_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SmallMultimap_O");
    core::Register_core__SmallMultimap_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SourceFileInfo_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SourceFileInfo_O");
    core::Register_core__SourceFileInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SourceManager_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SourceManager_O");
    core::Register_core__SourceManager_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SourcePosInfo_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SourcePosInfo_O");
    core::Register_core__SourcePosInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__StandardObject_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StandardObject_O");
    core::Register_core__StandardObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Stream_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Stream_O");
    core::Register_core__Stream_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__StructureObject_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StructureObject_O");
    core::Register_core__StructureObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SymbolToEnumConverter_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SymbolToEnumConverter_O");
    core::Register_core__SymbolToEnumConverter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Symbol_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Symbol_O");
    core::Register_core__Symbol_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__VaList_dummy_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__VaList_dummy_O");
    core::Register_core__VaList_dummy_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__WeakHashTable_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__WeakHashTable_O");
    core::Register_core__WeakHashTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__WeakKeyMapping_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__WeakKeyMapping_O");
    core::Register_core__WeakKeyMapping_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__WeakPointer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__WeakPointer_O");
    core::Register_core__WeakPointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__WrappedPointer_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__WrappedPointer_O");
    core::Register_core__WrappedPointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__AttributeSet_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__AttributeSet_O");
    llvmo::Register_llvmo__AttributeSet_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Attribute_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Attribute_O");
    llvmo::Register_llvmo__Attribute_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DebugInfo_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DebugInfo_O");
    llvmo::Register_llvmo__DebugInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DebugLoc_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DebugLoc_O");
    llvmo::Register_llvmo__DebugLoc_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__InsertPoint_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__InsertPoint_O");
    llvmo::Register_llvmo__InsertPoint_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__IterateCons_O(core::Lisp_sp); // base(s): set(['core::Iterator_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__IterateCons_O");
    adapt::Register_adapt__IterateCons_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__QDomNode_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__QDomNode_O");
    adapt::Register_adapt__QDomNode_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__StringSet_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__StringSet_O");
    adapt::Register_adapt__StringSet_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__SymbolList_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__SymbolList_O");
    adapt::Register_adapt__SymbolList_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ActivationFrame_O(core::Lisp_sp); // base(s): set(['core::Environment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ActivationFrame_O");
    core::Register_core__ActivationFrame_O(_lisp); // base(s): set(['core::Environment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ExtPkg
#ifdef EXTERN_REGISTER
extern void Register_core__AnsiStream_O(core::Lisp_sp); // base(s): set(['core::Stream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__AnsiStream_O");
    core::Register_core__AnsiStream_O(_lisp); // base(s): set(['core::Stream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ExtPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ArrayDisplaced_O(core::Lisp_sp); // base(s): set(['core::Array_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ArrayDisplaced_O");
    core::Register_core__ArrayDisplaced_O(_lisp); // base(s): set(['core::Array_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ArrayObjects_O(core::Lisp_sp); // base(s): set(['core::Array_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ArrayObjects_O");
    core::Register_core__ArrayObjects_O(_lisp); // base(s): set(['core::Array_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__BranchSNode_O(core::Lisp_sp); // base(s): set(['core::SNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__BranchSNode_O");
    core::Register_core__BranchSNode_O(_lisp); // base(s): set(['core::SNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__CompiledFunction_O(core::Lisp_sp); // base(s): set(['core::Function_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__CompiledFunction_O");
    core::Register_core__CompiledFunction_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Complex_O(core::Lisp_sp); // base(s): set(['core::Number_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Complex_O");
    core::Register_core__Complex_O(_lisp); // base(s): set(['core::Number_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__DirectoryIterator_O(core::Lisp_sp); // base(s): set(['core::Iterator_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__DirectoryIterator_O");
    core::Register_core__DirectoryIterator_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ForeignData_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ForeignData_O");
    core::Register_core__ForeignData_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__GlueEnvironment_O(core::Lisp_sp); // base(s): set(['core::Environment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__GlueEnvironment_O");
    core::Register_core__GlueEnvironment_O(_lisp); // base(s): set(['core::Environment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__HashTableEq_O(core::Lisp_sp); // base(s): set(['core::HashTable_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__HashTableEq_O");
    core::Register_core__HashTableEq_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__HashTableEql_O(core::Lisp_sp); // base(s): set(['core::HashTable_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__HashTableEql_O");
    core::Register_core__HashTableEql_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__HashTableEqual_O(core::Lisp_sp); // base(s): set(['core::HashTable_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__HashTableEqual_O");
    core::Register_core__HashTableEqual_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__HashTableEqualp_O(core::Lisp_sp); // base(s): set(['core::HashTable_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__HashTableEqualp_O");
    core::Register_core__HashTableEqualp_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Instance_O(core::Lisp_sp); // base(s): set(['core::Function_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Instance_O");
    core::Register_core__Instance_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__LeafSNode_O(core::Lisp_sp); // base(s): set(['core::SNode_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LeafSNode_O");
    core::Register_core__LeafSNode_O(_lisp); // base(s): set(['core::SNode_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__LexicalEnvironment_O(core::Lisp_sp); // base(s): set(['core::Environment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LexicalEnvironment_O");
    core::Register_core__LexicalEnvironment_O(_lisp); // base(s): set(['core::Environment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__LoadArchive_O(core::Lisp_sp); // base(s): set(['core::Archive_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LoadArchive_O");
    core::Register_core__LoadArchive_O(_lisp); // base(s): set(['core::Archive_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__LogicalPathname_O(core::Lisp_sp); // base(s): set(['core::Pathname_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LogicalPathname_O");
    core::Register_core__LogicalPathname_O(_lisp); // base(s): set(['core::Pathname_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Metaobject_O(core::Lisp_sp); // base(s): set(['core::StandardObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Metaobject_O");
    core::Register_core__Metaobject_O(_lisp); // base(s): set(['core::StandardObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Null_O(core::Lisp_sp); // base(s): set(['core::Symbol_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Null_O");
    core::Register_core__Null_O(_lisp); // base(s): set(['core::Symbol_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Real_O(core::Lisp_sp); // base(s): set(['core::Number_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Real_O");
    core::Register_core__Real_O(_lisp); // base(s): set(['core::Number_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__RecursiveDirectoryIterator_O(core::Lisp_sp); // base(s): set(['core::Iterator_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__RecursiveDirectoryIterator_O");
    core::Register_core__RecursiveDirectoryIterator_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SaveArchive_O(core::Lisp_sp); // base(s): set(['core::Archive_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SaveArchive_O");
    core::Register_core__SaveArchive_O(_lisp); // base(s): set(['core::Archive_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SingleDispatchEffectiveMethodFunction_O(core::Lisp_sp); // base(s): set(['core::Function_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SingleDispatchEffectiveMethodFunction_O");
    core::Register_core__SingleDispatchEffectiveMethodFunction_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SingleDispatchGenericFunction_O(core::Lisp_sp); // base(s): set(['core::Function_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SingleDispatchGenericFunction_O");
    core::Register_core__SingleDispatchGenericFunction_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SpecialForm_O(core::Lisp_sp); // base(s): set(['core::Function_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SpecialForm_O");
    core::Register_core__SpecialForm_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__UserData_O(core::Lisp_sp); // base(s): set(['core::LightUserData_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__UserData_O");
    core::Register_core__UserData_O(_lisp); // base(s): set(['core::LightUserData_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Vector_O(core::Lisp_sp); // base(s): set(['core::Array_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Vector_O");
    core::Register_core__Vector_O(_lisp); // base(s): set(['core::Array_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__WeakKeyHashTable_O(core::Lisp_sp); // base(s): set(['core::WeakHashTable_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__WeakKeyHashTable_O");
    core::Register_core__WeakKeyHashTable_O(_lisp); // base(s): set(['core::WeakHashTable_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__BoundingBox_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__BoundingBox_O");
    geom::Register_geom__BoundingBox_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__Color_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__Color_O");
    geom::Register_geom__Color_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__OMatrix_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__OMatrix_O");
    geom::Register_geom__OMatrix_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__OVector2_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__OVector2_O");
    geom::Register_geom__OVector2_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__OVector3_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__OVector3_O");
    geom::Register_geom__OVector3_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__APFloat_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__APFloat_O");
    llvmo::Register_llvmo__APFloat_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__APInt_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__APInt_O");
    llvmo::Register_llvmo__APInt_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIArray_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIArray_O");
    llvmo::Register_llvmo__DIArray_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIBasicType_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIBasicType_O");
    llvmo::Register_llvmo__DIBasicType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIBuilder_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIBuilder_O");
    llvmo::Register_llvmo__DIBuilder_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DICompileUnit_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DICompileUnit_O");
    llvmo::Register_llvmo__DICompileUnit_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DICompositeType_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DICompositeType_O");
    llvmo::Register_llvmo__DICompositeType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIDerivedType_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIDerivedType_O");
    llvmo::Register_llvmo__DIDerivedType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIDescriptor_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIDescriptor_O");
    llvmo::Register_llvmo__DIDescriptor_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIFile_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIFile_O");
    llvmo::Register_llvmo__DIFile_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DILexicalBlock_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DILexicalBlock_O");
    llvmo::Register_llvmo__DILexicalBlock_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIScope_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIScope_O");
    llvmo::Register_llvmo__DIScope_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DISubprogram_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DISubprogram_O");
    llvmo::Register_llvmo__DISubprogram_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DISubroutineType_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DISubroutineType_O");
    llvmo::Register_llvmo__DISubroutineType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DITypeArray_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DITypeArray_O");
    llvmo::Register_llvmo__DITypeArray_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DIType_O(core::Lisp_sp); // base(s): set(['llvmo::DebugInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DIType_O");
    llvmo::Register_llvmo__DIType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DataLayout_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DataLayout_O");
    llvmo::Register_llvmo__DataLayout_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__EngineBuilder_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__EngineBuilder_O");
    llvmo::Register_llvmo__EngineBuilder_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ExecutionEngine_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ExecutionEngine_O");
    llvmo::Register_llvmo__ExecutionEngine_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__IRBuilderBase_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__IRBuilderBase_O");
    llvmo::Register_llvmo__IRBuilderBase_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__LLVMContext_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__LLVMContext_O");
    llvmo::Register_llvmo__LLVMContext_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Linker_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Linker_O");
    llvmo::Register_llvmo__Linker_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__MCSubtargetInfo_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__MCSubtargetInfo_O");
    llvmo::Register_llvmo__MCSubtargetInfo_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Metadata_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Metadata_O");
    llvmo::Register_llvmo__Metadata_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Module_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Module_O");
    llvmo::Register_llvmo__Module_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__NamedMDNode_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__NamedMDNode_O");
    llvmo::Register_llvmo__NamedMDNode_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__PassManagerBase_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__PassManagerBase_O");
    llvmo::Register_llvmo__PassManagerBase_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__PassManagerBuilder_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__PassManagerBuilder_O");
    llvmo::Register_llvmo__PassManagerBuilder_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Pass_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Pass_O");
    llvmo::Register_llvmo__Pass_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__TargetMachine_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__TargetMachine_O");
    llvmo::Register_llvmo__TargetMachine_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__TargetOptions_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__TargetOptions_O");
    llvmo::Register_llvmo__TargetOptions_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Target_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Target_O");
    llvmo::Register_llvmo__Target_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Triple_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Triple_O");
    llvmo::Register_llvmo__Triple_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Type_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Type_O");
    llvmo::Register_llvmo__Type_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Value_O(core::Lisp_sp); // base(s): set(['core::ExternalObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Value_O");
    llvmo::Register_llvmo__Value_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__Dimension_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__Dimension_O");
    units::Register_units__Dimension_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__Quantity_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__Quantity_O");
    units::Register_units__Quantity_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__Unit_O(core::Lisp_sp); // base(s): set(['core::CxxObject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__Unit_O");
    units::Register_units__Unit_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__BitVector_O(core::Lisp_sp); // base(s): set(['core::Vector_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__BitVector_O");
    core::Register_core__BitVector_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__BroadcastStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__BroadcastStream_O");
    core::Register_core__BroadcastStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__CompileTimeEnvironment_O(core::Lisp_sp); // base(s): set(['core::LexicalEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__CompileTimeEnvironment_O");
    core::Register_core__CompileTimeEnvironment_O(_lisp); // base(s): set(['core::LexicalEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__ConcatenatedStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ConcatenatedStream_O");
    core::Register_core__ConcatenatedStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__EchoStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__EchoStream_O");
    core::Register_core__EchoStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__FileStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__FileStream_O");
    core::Register_core__FileStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Float_O(core::Lisp_sp); // base(s): set(['core::Real_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Float_O");
    core::Register_core__Float_O(_lisp); // base(s): set(['core::Real_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__FunctionFrame_O(core::Lisp_sp); // base(s): set(['core::ActivationFrame_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__FunctionFrame_O");
    core::Register_core__FunctionFrame_O(_lisp); // base(s): set(['core::ActivationFrame_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Rational_O(core::Lisp_sp); // base(s): set(['core::Real_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Rational_O");
    core::Register_core__Rational_O(_lisp); // base(s): set(['core::Real_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__RuntimeVisibleEnvironment_O(core::Lisp_sp); // base(s): set(['core::LexicalEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__RuntimeVisibleEnvironment_O");
    core::Register_core__RuntimeVisibleEnvironment_O(_lisp); // base(s): set(['core::LexicalEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SexpLoadArchive_O(core::Lisp_sp); // base(s): set(['core::LoadArchive_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SexpLoadArchive_O");
    core::Register_core__SexpLoadArchive_O(_lisp); // base(s): set(['core::LoadArchive_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SexpSaveArchive_O(core::Lisp_sp); // base(s): set(['core::SaveArchive_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SexpSaveArchive_O");
    core::Register_core__SexpSaveArchive_O(_lisp); // base(s): set(['core::SaveArchive_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__Specializer_O(core::Lisp_sp); // base(s): set(['core::Metaobject_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Specializer_O");
    core::Register_core__Specializer_O(_lisp); // base(s): set(['core::Metaobject_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__StringStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StringStream_O");
    core::Register_core__StringStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__String_O(core::Lisp_sp); // base(s): set(['core::Vector_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__String_O");
    core::Register_core__String_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__SynonymStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SynonymStream_O");
    core::Register_core__SynonymStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__TagbodyFrame_O(core::Lisp_sp); // base(s): set(['core::ActivationFrame_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__TagbodyFrame_O");
    core::Register_core__TagbodyFrame_O(_lisp); // base(s): set(['core::ActivationFrame_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__TwoWayStream_O(core::Lisp_sp); // base(s): set(['core::AnsiStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__TwoWayStream_O");
    core::Register_core__TwoWayStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ValueFrame_O(core::Lisp_sp); // base(s): set(['core::ActivationFrame_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ValueFrame_O");
    core::Register_core__ValueFrame_O(_lisp); // base(s): set(['core::ActivationFrame_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__VectorObjects_O(core::Lisp_sp); // base(s): set(['core::Vector_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__VectorObjects_O");
    core::Register_core__VectorObjects_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__CoordinateArray_O(core::Lisp_sp); // base(s): set(['core::Vector_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__CoordinateArray_O");
    geom::Register_geom__CoordinateArray_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Argument_O(core::Lisp_sp); // base(s): set(['llvmo::Value_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Argument_O");
    llvmo::Register_llvmo__Argument_O(_lisp); // base(s): set(['llvmo::Value_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__BasicBlock_O(core::Lisp_sp); // base(s): set(['llvmo::Value_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__BasicBlock_O");
    llvmo::Register_llvmo__BasicBlock_O(_lisp); // base(s): set(['llvmo::Value_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__CompositeType_O(core::Lisp_sp); // base(s): set(['llvmo::Type_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__CompositeType_O");
    llvmo::Register_llvmo__CompositeType_O(_lisp); // base(s): set(['llvmo::Type_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__FunctionPassManager_O(core::Lisp_sp); // base(s): set(['llvmo::PassManagerBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__FunctionPassManager_O");
    llvmo::Register_llvmo__FunctionPassManager_O(_lisp); // base(s): set(['llvmo::PassManagerBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__FunctionPass_O(core::Lisp_sp); // base(s): set(['llvmo::Pass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__FunctionPass_O");
    llvmo::Register_llvmo__FunctionPass_O(_lisp); // base(s): set(['llvmo::Pass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__FunctionType_O(core::Lisp_sp); // base(s): set(['llvmo::Type_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__FunctionType_O");
    llvmo::Register_llvmo__FunctionType_O(_lisp); // base(s): set(['llvmo::Type_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__IRBuilder_O(core::Lisp_sp); // base(s): set(['llvmo::IRBuilderBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__IRBuilder_O");
    llvmo::Register_llvmo__IRBuilder_O(_lisp); // base(s): set(['llvmo::IRBuilderBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__IntegerType_O(core::Lisp_sp); // base(s): set(['llvmo::Type_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__IntegerType_O");
    llvmo::Register_llvmo__IntegerType_O(_lisp); // base(s): set(['llvmo::Type_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__LLVMTargetMachine_O(core::Lisp_sp); // base(s): set(['llvmo::TargetMachine_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__LLVMTargetMachine_O");
    llvmo::Register_llvmo__LLVMTargetMachine_O(_lisp); // base(s): set(['llvmo::TargetMachine_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__MDNode_O(core::Lisp_sp); // base(s): set(['llvmo::Metadata_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__MDNode_O");
    llvmo::Register_llvmo__MDNode_O(_lisp); // base(s): set(['llvmo::Metadata_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__MDString_O(core::Lisp_sp); // base(s): set(['llvmo::Metadata_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__MDString_O");
    llvmo::Register_llvmo__MDString_O(_lisp); // base(s): set(['llvmo::Metadata_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ModulePass_O(core::Lisp_sp); // base(s): set(['llvmo::Pass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ModulePass_O");
    llvmo::Register_llvmo__ModulePass_O(_lisp); // base(s): set(['llvmo::Pass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__PassManager_O(core::Lisp_sp); // base(s): set(['llvmo::PassManagerBase_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__PassManager_O");
    llvmo::Register_llvmo__PassManager_O(_lisp); // base(s): set(['llvmo::PassManagerBase_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__TargetSubtargetInfo_O(core::Lisp_sp); // base(s): set(['llvmo::MCSubtargetInfo_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__TargetSubtargetInfo_O");
    llvmo::Register_llvmo__TargetSubtargetInfo_O(_lisp); // base(s): set(['llvmo::MCSubtargetInfo_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__User_O(core::Lisp_sp); // base(s): set(['llvmo::Value_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__User_O");
    llvmo::Register_llvmo__User_O(_lisp); // base(s): set(['llvmo::Value_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ValueAsMetadata_O(core::Lisp_sp); // base(s): set(['llvmo::Metadata_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ValueAsMetadata_O");
    llvmo::Register_llvmo__ValueAsMetadata_O(_lisp); // base(s): set(['llvmo::Metadata_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__NamedUnit_O(core::Lisp_sp); // base(s): set(['units::Unit_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__NamedUnit_O");
    units::Register_units__NamedUnit_O(_lisp); // base(s): set(['units::Unit_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__BlockEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__BlockEnvironment_O");
    core::Register_core__BlockEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__CatchEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__CatchEnvironment_O");
    core::Register_core__CatchEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Class_O(core::Lisp_sp); // base(s): set(['core::Specializer_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Class_O");
    core::Register_core__Class_O(_lisp); // base(s): set(['core::Specializer_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__DoubleFloat_O(core::Lisp_sp); // base(s): set(['core::Float_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__DoubleFloat_O");
    core::Register_core__DoubleFloat_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__FunctionContainerEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__FunctionContainerEnvironment_O");
    core::Register_core__FunctionContainerEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__FunctionValueEnvironment_O(core::Lisp_sp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__FunctionValueEnvironment_O");
    core::Register_core__FunctionValueEnvironment_O(_lisp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__IOFileStream_O(core::Lisp_sp); // base(s): set(['core::FileStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__IOFileStream_O");
    core::Register_core__IOFileStream_O(_lisp); // base(s): set(['core::FileStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__IOStreamStream_O(core::Lisp_sp); // base(s): set(['core::FileStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__IOStreamStream_O");
    core::Register_core__IOStreamStream_O(_lisp); // base(s): set(['core::FileStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Integer_O(core::Lisp_sp); // base(s): set(['core::Rational_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Integer_O");
    core::Register_core__Integer_O(_lisp); // base(s): set(['core::Rational_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__LongFloat_O(core::Lisp_sp); // base(s): set(['core::Float_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__LongFloat_O");
    core::Register_core__LongFloat_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__MacroletEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__MacroletEnvironment_O");
    core::Register_core__MacroletEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Ratio_O(core::Lisp_sp); // base(s): set(['core::Rational_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Ratio_O");
    core::Register_core__Ratio_O(_lisp); // base(s): set(['core::Rational_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__ShortFloat_O(core::Lisp_sp); // base(s): set(['core::Float_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ShortFloat_O");
    core::Register_core__ShortFloat_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__SimpleBitVector_O(core::Lisp_sp); // base(s): set(['core::BitVector_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SimpleBitVector_O");
    core::Register_core__SimpleBitVector_O(_lisp); // base(s): set(['core::BitVector_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__SingleFloat_dummy_O(core::Lisp_sp); // base(s): set(['core::Float_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SingleFloat_dummy_O");
    core::Register_core__SingleFloat_dummy_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__StackValueEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StackValueEnvironment_O");
    core::Register_core__StackValueEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Str_O(core::Lisp_sp); // base(s): set(['core::String_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Str_O");
    core::Register_core__Str_O(_lisp); // base(s): set(['core::String_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__StringInputStream_O(core::Lisp_sp); // base(s): set(['core::StringStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StringInputStream_O");
    core::Register_core__StringInputStream_O(_lisp); // base(s): set(['core::StringStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__StringOutputStream_O(core::Lisp_sp); // base(s): set(['core::StringStream_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StringOutputStream_O");
    core::Register_core__StringOutputStream_O(_lisp); // base(s): set(['core::StringStream_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__SymbolMacroletEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__SymbolMacroletEnvironment_O");
    core::Register_core__SymbolMacroletEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__TagbodyEnvironment_O(core::Lisp_sp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__TagbodyEnvironment_O");
    core::Register_core__TagbodyEnvironment_O(_lisp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__UnwindProtectEnvironment_O(core::Lisp_sp); // base(s): set(['core::CompileTimeEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__UnwindProtectEnvironment_O");
    core::Register_core__UnwindProtectEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ValueEnvironment_O(core::Lisp_sp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ValueEnvironment_O");
    core::Register_core__ValueEnvironment_O(_lisp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__VectorObjectsWithFillPtr_O(core::Lisp_sp); // base(s): set(['core::VectorObjects_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__VectorObjectsWithFillPtr_O");
    core::Register_core__VectorObjectsWithFillPtr_O(_lisp); // base(s): set(['core::VectorObjects_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__CoordinateArrayWithHash_O(core::Lisp_sp); // base(s): set(['geom::CoordinateArray_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__CoordinateArrayWithHash_O");
    geom::Register_geom__CoordinateArrayWithHash_O(_lisp); // base(s): set(['geom::CoordinateArray_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Constant_O(core::Lisp_sp); // base(s): set(['llvmo::User_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Constant_O");
    llvmo::Register_llvmo__Constant_O(_lisp); // base(s): set(['llvmo::User_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ImmutablePass_O(core::Lisp_sp); // base(s): set(['llvmo::ModulePass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ImmutablePass_O");
    llvmo::Register_llvmo__ImmutablePass_O(_lisp); // base(s): set(['llvmo::ModulePass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Instruction_O(core::Lisp_sp); // base(s): set(['llvmo::User_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Instruction_O");
    llvmo::Register_llvmo__Instruction_O(_lisp); // base(s): set(['llvmo::User_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__SequentialType_O(core::Lisp_sp); // base(s): set(['llvmo::CompositeType_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__SequentialType_O");
    llvmo::Register_llvmo__SequentialType_O(_lisp); // base(s): set(['llvmo::CompositeType_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__StructType_O(core::Lisp_sp); // base(s): set(['llvmo::CompositeType_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__StructType_O");
    llvmo::Register_llvmo__StructType_O(_lisp); // base(s): set(['llvmo::CompositeType_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Bignum_O(core::Lisp_sp); // base(s): set(['core::Integer_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Bignum_O");
    core::Register_core__Bignum_O(_lisp); // base(s): set(['core::Integer_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__BuiltInClass_O(core::Lisp_sp); // base(s): set(['core::Class_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__BuiltInClass_O");
    core::Register_core__BuiltInClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__CxxClass_O(core::Lisp_sp); // base(s): set(['core::Class_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__CxxClass_O");
    core::Register_core__CxxClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__Fixnum_dummy_O(core::Lisp_sp); // base(s): set(['core::Integer_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__Fixnum_dummy_O");
    core::Register_core__Fixnum_dummy_O(_lisp); // base(s): set(['core::Integer_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__ForwardReferencedClass_O(core::Lisp_sp); // base(s): set(['core::Class_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__ForwardReferencedClass_O");
    core::Register_core__ForwardReferencedClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClosPkg
#ifdef EXTERN_REGISTER
extern void Register_core__StdClass_O(core::Lisp_sp); // base(s): set(['core::Class_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StdClass_O");
    core::Register_core__StdClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClosPkg
#ifdef Use_CorePkg
#ifdef EXTERN_REGISTER
extern void Register_core__StrWithFillPtr_O(core::Lisp_sp); // base(s): set(['core::Str_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StrWithFillPtr_O");
    core::Register_core__StrWithFillPtr_O(_lisp); // base(s): set(['core::Str_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__StructureClass_O(core::Lisp_sp); // base(s): set(['core::Class_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StructureClass_O");
    core::Register_core__StructureClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ArrayType_O(core::Lisp_sp); // base(s): set(['llvmo::SequentialType_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ArrayType_O");
    llvmo::Register_llvmo__ArrayType_O(_lisp); // base(s): set(['llvmo::SequentialType_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__AtomicCmpXchgInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__AtomicCmpXchgInst_O");
    llvmo::Register_llvmo__AtomicCmpXchgInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__AtomicRMWInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__AtomicRMWInst_O");
    llvmo::Register_llvmo__AtomicRMWInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__BlockAddress_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__BlockAddress_O");
    llvmo::Register_llvmo__BlockAddress_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__CallInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__CallInst_O");
    llvmo::Register_llvmo__CallInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantArray_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantArray_O");
    llvmo::Register_llvmo__ConstantArray_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantDataSequential_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantDataSequential_O");
    llvmo::Register_llvmo__ConstantDataSequential_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantExpr_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantExpr_O");
    llvmo::Register_llvmo__ConstantExpr_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantFP_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantFP_O");
    llvmo::Register_llvmo__ConstantFP_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantInt_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantInt_O");
    llvmo::Register_llvmo__ConstantInt_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantPointerNull_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantPointerNull_O");
    llvmo::Register_llvmo__ConstantPointerNull_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantStruct_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantStruct_O");
    llvmo::Register_llvmo__ConstantStruct_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__DataLayoutPass_O(core::Lisp_sp); // base(s): set(['llvmo::ImmutablePass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__DataLayoutPass_O");
    llvmo::Register_llvmo__DataLayoutPass_O(_lisp); // base(s): set(['llvmo::ImmutablePass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__FenceInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__FenceInst_O");
    llvmo::Register_llvmo__FenceInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__GlobalValue_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__GlobalValue_O");
    llvmo::Register_llvmo__GlobalValue_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__LandingPadInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__LandingPadInst_O");
    llvmo::Register_llvmo__LandingPadInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__PHINode_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__PHINode_O");
    llvmo::Register_llvmo__PHINode_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__PointerType_O(core::Lisp_sp); // base(s): set(['llvmo::SequentialType_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__PointerType_O");
    llvmo::Register_llvmo__PointerType_O(_lisp); // base(s): set(['llvmo::SequentialType_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__StoreInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__StoreInst_O");
    llvmo::Register_llvmo__StoreInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__TargetLibraryInfo_O(core::Lisp_sp); // base(s): set(['llvmo::ImmutablePass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__TargetLibraryInfo_O");
    llvmo::Register_llvmo__TargetLibraryInfo_O(_lisp); // base(s): set(['llvmo::ImmutablePass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__TerminatorInst_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__TerminatorInst_O");
    llvmo::Register_llvmo__TerminatorInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__UnaryInstruction_O(core::Lisp_sp); // base(s): set(['llvmo::Instruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__UnaryInstruction_O");
    llvmo::Register_llvmo__UnaryInstruction_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__UndefValue_O(core::Lisp_sp); // base(s): set(['llvmo::Constant_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__UndefValue_O");
    llvmo::Register_llvmo__UndefValue_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__VectorType_O(core::Lisp_sp); // base(s): set(['llvmo::SequentialType_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__VectorType_O");
    llvmo::Register_llvmo__VectorType_O(_lisp); // base(s): set(['llvmo::SequentialType_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_ClbindPkg
#ifdef EXTERN_REGISTER
extern void Register_clbind__ClassRep_O(core::Lisp_sp); // base(s): set(['core::BuiltInClass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_clbind__ClassRep_O");
    clbind::Register_clbind__ClassRep_O(_lisp); // base(s): set(['core::BuiltInClass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClbindPkg
#ifdef Use_ClosPkg
#ifdef EXTERN_REGISTER
extern void Register_core__FuncallableStandardClass_O(core::Lisp_sp); // base(s): set(['core::StdClass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__FuncallableStandardClass_O");
    core::Register_core__FuncallableStandardClass_O(_lisp); // base(s): set(['core::StdClass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClosPkg
#ifdef Use_ClPkg
#ifdef EXTERN_REGISTER
extern void Register_core__StandardClass_O(core::Lisp_sp); // base(s): set(['core::StdClass_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_core__StandardClass_O");
    core::Register_core__StandardClass_O(_lisp); // base(s): set(['core::StdClass_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_ClPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__AllocaInst_O(core::Lisp_sp); // base(s): set(['llvmo::UnaryInstruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__AllocaInst_O");
    llvmo::Register_llvmo__AllocaInst_O(_lisp); // base(s): set(['llvmo::UnaryInstruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__BranchInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__BranchInst_O");
    llvmo::Register_llvmo__BranchInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ConstantDataArray_O(core::Lisp_sp); // base(s): set(['llvmo::ConstantDataSequential_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ConstantDataArray_O");
    llvmo::Register_llvmo__ConstantDataArray_O(_lisp); // base(s): set(['llvmo::ConstantDataSequential_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__Function_O(core::Lisp_sp); // base(s): set(['llvmo::GlobalValue_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__Function_O");
    llvmo::Register_llvmo__Function_O(_lisp); // base(s): set(['llvmo::GlobalValue_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__GlobalVariable_O(core::Lisp_sp); // base(s): set(['llvmo::GlobalValue_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__GlobalVariable_O");
    llvmo::Register_llvmo__GlobalVariable_O(_lisp); // base(s): set(['llvmo::GlobalValue_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__IndirectBrInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__IndirectBrInst_O");
    llvmo::Register_llvmo__IndirectBrInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__InvokeInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__InvokeInst_O");
    llvmo::Register_llvmo__InvokeInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__LoadInst_O(core::Lisp_sp); // base(s): set(['llvmo::UnaryInstruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__LoadInst_O");
    llvmo::Register_llvmo__LoadInst_O(_lisp); // base(s): set(['llvmo::UnaryInstruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ResumeInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ResumeInst_O");
    llvmo::Register_llvmo__ResumeInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__ReturnInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__ReturnInst_O");
    llvmo::Register_llvmo__ReturnInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__SwitchInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__SwitchInst_O");
    llvmo::Register_llvmo__SwitchInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__UnreachableInst_O(core::Lisp_sp); // base(s): set(['llvmo::TerminatorInst_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__UnreachableInst_O");
    llvmo::Register_llvmo__UnreachableInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
#ifdef EXTERN_REGISTER
extern void Register_llvmo__VAArgInst_O(core::Lisp_sp); // base(s): set(['llvmo::UnaryInstruction_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_llvmo__VAArgInst_O");
    llvmo::Register_llvmo__VAArgInst_O(_lisp); // base(s): set(['llvmo::UnaryInstruction_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_LlvmoPkg
#endif // EXPOSE_TO_CANDO
#undef EXPOSE_TO_CANDO
#ifdef EXPOSE_TO_PYTHON
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__ObjectSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__ObjectSet_O");
	Call_exposePython_adapt__ObjectSet_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__StringList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__StringList_O");
	Call_exposePython_adapt__StringList_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__SymbolSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__SymbolSet_O");
	Call_exposePython_adapt__SymbolSet_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__T_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__T_O");
	Call_exposePython_core__T_O(_lisp); // base(s): set(['core::_RootDummyClass'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_AstToolingPkg
extern void Call_exposePython_asttooling__AstVisitor_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: asttooling__AstVisitor_O");
	Call_exposePython_asttooling__AstVisitor_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_AstToolingPkg
#ifdef Use_CffiPkg
extern void Call_exposePython_cffi__Pointer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: cffi__Pointer_O");
	Call_exposePython_cffi__Pointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CffiPkg
#ifdef Use_ClbindPkg
extern void Call_exposePython_clbind__ClassRegistry_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: clbind__ClassRegistry_O");
	Call_exposePython_clbind__ClassRegistry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClbindPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Archive_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Archive_O");
	Call_exposePython_core__Archive_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Array_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Array_O");
	Call_exposePython_core__Array_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Binder_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Binder_O");
	Call_exposePython_core__Binder_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__CandoException_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__CandoException_O");
	Call_exposePython_core__CandoException_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Character_dummy_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Character_dummy_O");
	Call_exposePython_core__Character_dummy_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Cons_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Cons_O");
	Call_exposePython_core__Cons_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__CxxObject_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__CxxObject_O");
	Call_exposePython_core__CxxObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__DirectoryEntry_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__DirectoryEntry_O");
	Call_exposePython_core__DirectoryEntry_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Environment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Environment_O");
	Call_exposePython_core__Environment_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ExternalObject_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ExternalObject_O");
	Call_exposePython_core__ExternalObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__FileStatus_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__FileStatus_O");
	Call_exposePython_core__FileStatus_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Function_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Function_O");
	Call_exposePython_core__Function_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__HashTable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__HashTable_O");
	Call_exposePython_core__HashTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__IntArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__IntArray_O");
	Call_exposePython_core__IntArray_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__InvocationHistoryFrameIterator_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__InvocationHistoryFrameIterator_O");
	Call_exposePython_core__InvocationHistoryFrameIterator_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Iterator_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Iterator_O");
	Call_exposePython_core__Iterator_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__LambdaListHandler_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LambdaListHandler_O");
	Call_exposePython_core__LambdaListHandler_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__LightUserData_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LightUserData_O");
	Call_exposePython_core__LightUserData_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__LoadTimeValues_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LoadTimeValues_O");
	Call_exposePython_core__LoadTimeValues_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__MultiStringBuffer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__MultiStringBuffer_O");
	Call_exposePython_core__MultiStringBuffer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Number_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Number_O");
	Call_exposePython_core__Number_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Package_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Package_O");
	Call_exposePython_core__Package_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Path_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Path_O");
	Call_exposePython_core__Path_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Pathname_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Pathname_O");
	Call_exposePython_core__Pathname_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Pointer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Pointer_O");
	Call_exposePython_core__Pointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__PosixTimeDuration_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__PosixTimeDuration_O");
	Call_exposePython_core__PosixTimeDuration_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__PosixTime_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__PosixTime_O");
	Call_exposePython_core__PosixTime_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__RandomState_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__RandomState_O");
	Call_exposePython_core__RandomState_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__ReadTable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ReadTable_O");
	Call_exposePython_core__ReadTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Reader_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Reader_O");
	Call_exposePython_core__Reader_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Record_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Record_O");
	Call_exposePython_core__Record_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__RegexMatch_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__RegexMatch_O");
	Call_exposePython_core__RegexMatch_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Regex_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Regex_O");
	Call_exposePython_core__Regex_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SNode_O");
	Call_exposePython_core__SNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SingleDispatchMethod_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SingleDispatchMethod_O");
	Call_exposePython_core__SingleDispatchMethod_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SmallMap_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SmallMap_O");
	Call_exposePython_core__SmallMap_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SmallMultimap_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SmallMultimap_O");
	Call_exposePython_core__SmallMultimap_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SourceFileInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SourceFileInfo_O");
	Call_exposePython_core__SourceFileInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SourceManager_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SourceManager_O");
	Call_exposePython_core__SourceManager_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SourcePosInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SourcePosInfo_O");
	Call_exposePython_core__SourcePosInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__StandardObject_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StandardObject_O");
	Call_exposePython_core__StandardObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Stream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Stream_O");
	Call_exposePython_core__Stream_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__StructureObject_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StructureObject_O");
	Call_exposePython_core__StructureObject_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SymbolToEnumConverter_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SymbolToEnumConverter_O");
	Call_exposePython_core__SymbolToEnumConverter_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Symbol_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Symbol_O");
	Call_exposePython_core__Symbol_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__VaList_dummy_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__VaList_dummy_O");
	Call_exposePython_core__VaList_dummy_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__WeakHashTable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__WeakHashTable_O");
	Call_exposePython_core__WeakHashTable_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__WeakKeyMapping_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__WeakKeyMapping_O");
	Call_exposePython_core__WeakKeyMapping_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__WeakPointer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__WeakPointer_O");
	Call_exposePython_core__WeakPointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__WrappedPointer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__WrappedPointer_O");
	Call_exposePython_core__WrappedPointer_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__AttributeSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__AttributeSet_O");
	Call_exposePython_llvmo__AttributeSet_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Attribute_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Attribute_O");
	Call_exposePython_llvmo__Attribute_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DebugInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DebugInfo_O");
	Call_exposePython_llvmo__DebugInfo_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DebugLoc_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DebugLoc_O");
	Call_exposePython_llvmo__DebugLoc_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__InsertPoint_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__InsertPoint_O");
	Call_exposePython_llvmo__InsertPoint_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__IterateCons_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__IterateCons_O");
	Call_exposePython_adapt__IterateCons_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__QDomNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__QDomNode_O");
	Call_exposePython_adapt__QDomNode_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__StringSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__StringSet_O");
	Call_exposePython_adapt__StringSet_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__SymbolList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__SymbolList_O");
	Call_exposePython_adapt__SymbolList_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ActivationFrame_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ActivationFrame_O");
	Call_exposePython_core__ActivationFrame_O(_lisp); // base(s): set(['core::Environment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ExtPkg
extern void Call_exposePython_core__AnsiStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__AnsiStream_O");
	Call_exposePython_core__AnsiStream_O(_lisp); // base(s): set(['core::Stream_O'])
}
#endif // ifdef Use_ExtPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ArrayDisplaced_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ArrayDisplaced_O");
	Call_exposePython_core__ArrayDisplaced_O(_lisp); // base(s): set(['core::Array_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ArrayObjects_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ArrayObjects_O");
	Call_exposePython_core__ArrayObjects_O(_lisp); // base(s): set(['core::Array_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__BranchSNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__BranchSNode_O");
	Call_exposePython_core__BranchSNode_O(_lisp); // base(s): set(['core::SNode_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__CompiledFunction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__CompiledFunction_O");
	Call_exposePython_core__CompiledFunction_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Complex_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Complex_O");
	Call_exposePython_core__Complex_O(_lisp); // base(s): set(['core::Number_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__DirectoryIterator_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__DirectoryIterator_O");
	Call_exposePython_core__DirectoryIterator_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ForeignData_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ForeignData_O");
	Call_exposePython_core__ForeignData_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__GlueEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__GlueEnvironment_O");
	Call_exposePython_core__GlueEnvironment_O(_lisp); // base(s): set(['core::Environment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__HashTableEq_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__HashTableEq_O");
	Call_exposePython_core__HashTableEq_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__HashTableEql_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__HashTableEql_O");
	Call_exposePython_core__HashTableEql_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__HashTableEqual_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__HashTableEqual_O");
	Call_exposePython_core__HashTableEqual_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__HashTableEqualp_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__HashTableEqualp_O");
	Call_exposePython_core__HashTableEqualp_O(_lisp); // base(s): set(['core::HashTable_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Instance_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Instance_O");
	Call_exposePython_core__Instance_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__LeafSNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LeafSNode_O");
	Call_exposePython_core__LeafSNode_O(_lisp); // base(s): set(['core::SNode_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__LexicalEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LexicalEnvironment_O");
	Call_exposePython_core__LexicalEnvironment_O(_lisp); // base(s): set(['core::Environment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__LoadArchive_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LoadArchive_O");
	Call_exposePython_core__LoadArchive_O(_lisp); // base(s): set(['core::Archive_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__LogicalPathname_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LogicalPathname_O");
	Call_exposePython_core__LogicalPathname_O(_lisp); // base(s): set(['core::Pathname_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Metaobject_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Metaobject_O");
	Call_exposePython_core__Metaobject_O(_lisp); // base(s): set(['core::StandardObject_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Null_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Null_O");
	Call_exposePython_core__Null_O(_lisp); // base(s): set(['core::Symbol_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Real_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Real_O");
	Call_exposePython_core__Real_O(_lisp); // base(s): set(['core::Number_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__RecursiveDirectoryIterator_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__RecursiveDirectoryIterator_O");
	Call_exposePython_core__RecursiveDirectoryIterator_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SaveArchive_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SaveArchive_O");
	Call_exposePython_core__SaveArchive_O(_lisp); // base(s): set(['core::Archive_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SingleDispatchEffectiveMethodFunction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SingleDispatchEffectiveMethodFunction_O");
	Call_exposePython_core__SingleDispatchEffectiveMethodFunction_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SingleDispatchGenericFunction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SingleDispatchGenericFunction_O");
	Call_exposePython_core__SingleDispatchGenericFunction_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SpecialForm_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SpecialForm_O");
	Call_exposePython_core__SpecialForm_O(_lisp); // base(s): set(['core::Function_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__UserData_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__UserData_O");
	Call_exposePython_core__UserData_O(_lisp); // base(s): set(['core::LightUserData_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Vector_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Vector_O");
	Call_exposePython_core__Vector_O(_lisp); // base(s): set(['core::Array_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__WeakKeyHashTable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__WeakKeyHashTable_O");
	Call_exposePython_core__WeakKeyHashTable_O(_lisp); // base(s): set(['core::WeakHashTable_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__BoundingBox_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__BoundingBox_O");
	Call_exposePython_geom__BoundingBox_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__Color_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__Color_O");
	Call_exposePython_geom__Color_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__OMatrix_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__OMatrix_O");
	Call_exposePython_geom__OMatrix_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__OVector2_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__OVector2_O");
	Call_exposePython_geom__OVector2_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__OVector3_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__OVector3_O");
	Call_exposePython_geom__OVector3_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__APFloat_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__APFloat_O");
	Call_exposePython_llvmo__APFloat_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__APInt_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__APInt_O");
	Call_exposePython_llvmo__APInt_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIArray_O");
	Call_exposePython_llvmo__DIArray_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIBasicType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIBasicType_O");
	Call_exposePython_llvmo__DIBasicType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIBuilder_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIBuilder_O");
	Call_exposePython_llvmo__DIBuilder_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DICompileUnit_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DICompileUnit_O");
	Call_exposePython_llvmo__DICompileUnit_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DICompositeType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DICompositeType_O");
	Call_exposePython_llvmo__DICompositeType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIDerivedType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIDerivedType_O");
	Call_exposePython_llvmo__DIDerivedType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIDescriptor_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIDescriptor_O");
	Call_exposePython_llvmo__DIDescriptor_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIFile_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIFile_O");
	Call_exposePython_llvmo__DIFile_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DILexicalBlock_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DILexicalBlock_O");
	Call_exposePython_llvmo__DILexicalBlock_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIScope_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIScope_O");
	Call_exposePython_llvmo__DIScope_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DISubprogram_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DISubprogram_O");
	Call_exposePython_llvmo__DISubprogram_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DISubroutineType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DISubroutineType_O");
	Call_exposePython_llvmo__DISubroutineType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DITypeArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DITypeArray_O");
	Call_exposePython_llvmo__DITypeArray_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DIType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DIType_O");
	Call_exposePython_llvmo__DIType_O(_lisp); // base(s): set(['llvmo::DebugInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DataLayout_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DataLayout_O");
	Call_exposePython_llvmo__DataLayout_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__EngineBuilder_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__EngineBuilder_O");
	Call_exposePython_llvmo__EngineBuilder_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ExecutionEngine_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ExecutionEngine_O");
	Call_exposePython_llvmo__ExecutionEngine_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__IRBuilderBase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__IRBuilderBase_O");
	Call_exposePython_llvmo__IRBuilderBase_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__LLVMContext_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__LLVMContext_O");
	Call_exposePython_llvmo__LLVMContext_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Linker_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Linker_O");
	Call_exposePython_llvmo__Linker_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__MCSubtargetInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__MCSubtargetInfo_O");
	Call_exposePython_llvmo__MCSubtargetInfo_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Metadata_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Metadata_O");
	Call_exposePython_llvmo__Metadata_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Module_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Module_O");
	Call_exposePython_llvmo__Module_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__NamedMDNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__NamedMDNode_O");
	Call_exposePython_llvmo__NamedMDNode_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__PassManagerBase_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__PassManagerBase_O");
	Call_exposePython_llvmo__PassManagerBase_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__PassManagerBuilder_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__PassManagerBuilder_O");
	Call_exposePython_llvmo__PassManagerBuilder_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Pass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Pass_O");
	Call_exposePython_llvmo__Pass_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__TargetMachine_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__TargetMachine_O");
	Call_exposePython_llvmo__TargetMachine_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__TargetOptions_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__TargetOptions_O");
	Call_exposePython_llvmo__TargetOptions_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Target_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Target_O");
	Call_exposePython_llvmo__Target_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Triple_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Triple_O");
	Call_exposePython_llvmo__Triple_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Type_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Type_O");
	Call_exposePython_llvmo__Type_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Value_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Value_O");
	Call_exposePython_llvmo__Value_O(_lisp); // base(s): set(['core::ExternalObject_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__Dimension_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__Dimension_O");
	Call_exposePython_units__Dimension_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__Quantity_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__Quantity_O");
	Call_exposePython_units__Quantity_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__Unit_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__Unit_O");
	Call_exposePython_units__Unit_O(_lisp); // base(s): set(['core::CxxObject_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__BitVector_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__BitVector_O");
	Call_exposePython_core__BitVector_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__BroadcastStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__BroadcastStream_O");
	Call_exposePython_core__BroadcastStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__CompileTimeEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__CompileTimeEnvironment_O");
	Call_exposePython_core__CompileTimeEnvironment_O(_lisp); // base(s): set(['core::LexicalEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__ConcatenatedStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ConcatenatedStream_O");
	Call_exposePython_core__ConcatenatedStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__EchoStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__EchoStream_O");
	Call_exposePython_core__EchoStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__FileStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__FileStream_O");
	Call_exposePython_core__FileStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Float_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Float_O");
	Call_exposePython_core__Float_O(_lisp); // base(s): set(['core::Real_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__FunctionFrame_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__FunctionFrame_O");
	Call_exposePython_core__FunctionFrame_O(_lisp); // base(s): set(['core::ActivationFrame_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Rational_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Rational_O");
	Call_exposePython_core__Rational_O(_lisp); // base(s): set(['core::Real_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__RuntimeVisibleEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__RuntimeVisibleEnvironment_O");
	Call_exposePython_core__RuntimeVisibleEnvironment_O(_lisp); // base(s): set(['core::LexicalEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SexpLoadArchive_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SexpLoadArchive_O");
	Call_exposePython_core__SexpLoadArchive_O(_lisp); // base(s): set(['core::LoadArchive_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SexpSaveArchive_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SexpSaveArchive_O");
	Call_exposePython_core__SexpSaveArchive_O(_lisp); // base(s): set(['core::SaveArchive_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__Specializer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Specializer_O");
	Call_exposePython_core__Specializer_O(_lisp); // base(s): set(['core::Metaobject_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__StringStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StringStream_O");
	Call_exposePython_core__StringStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__String_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__String_O");
	Call_exposePython_core__String_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__SynonymStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SynonymStream_O");
	Call_exposePython_core__SynonymStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__TagbodyFrame_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__TagbodyFrame_O");
	Call_exposePython_core__TagbodyFrame_O(_lisp); // base(s): set(['core::ActivationFrame_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__TwoWayStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__TwoWayStream_O");
	Call_exposePython_core__TwoWayStream_O(_lisp); // base(s): set(['core::AnsiStream_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ValueFrame_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ValueFrame_O");
	Call_exposePython_core__ValueFrame_O(_lisp); // base(s): set(['core::ActivationFrame_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__VectorObjects_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__VectorObjects_O");
	Call_exposePython_core__VectorObjects_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__CoordinateArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__CoordinateArray_O");
	Call_exposePython_geom__CoordinateArray_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Argument_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Argument_O");
	Call_exposePython_llvmo__Argument_O(_lisp); // base(s): set(['llvmo::Value_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__BasicBlock_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__BasicBlock_O");
	Call_exposePython_llvmo__BasicBlock_O(_lisp); // base(s): set(['llvmo::Value_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__CompositeType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__CompositeType_O");
	Call_exposePython_llvmo__CompositeType_O(_lisp); // base(s): set(['llvmo::Type_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__FunctionPassManager_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__FunctionPassManager_O");
	Call_exposePython_llvmo__FunctionPassManager_O(_lisp); // base(s): set(['llvmo::PassManagerBase_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__FunctionPass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__FunctionPass_O");
	Call_exposePython_llvmo__FunctionPass_O(_lisp); // base(s): set(['llvmo::Pass_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__FunctionType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__FunctionType_O");
	Call_exposePython_llvmo__FunctionType_O(_lisp); // base(s): set(['llvmo::Type_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__IRBuilder_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__IRBuilder_O");
	Call_exposePython_llvmo__IRBuilder_O(_lisp); // base(s): set(['llvmo::IRBuilderBase_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__IntegerType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__IntegerType_O");
	Call_exposePython_llvmo__IntegerType_O(_lisp); // base(s): set(['llvmo::Type_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__LLVMTargetMachine_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__LLVMTargetMachine_O");
	Call_exposePython_llvmo__LLVMTargetMachine_O(_lisp); // base(s): set(['llvmo::TargetMachine_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__MDNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__MDNode_O");
	Call_exposePython_llvmo__MDNode_O(_lisp); // base(s): set(['llvmo::Metadata_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__MDString_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__MDString_O");
	Call_exposePython_llvmo__MDString_O(_lisp); // base(s): set(['llvmo::Metadata_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ModulePass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ModulePass_O");
	Call_exposePython_llvmo__ModulePass_O(_lisp); // base(s): set(['llvmo::Pass_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__PassManager_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__PassManager_O");
	Call_exposePython_llvmo__PassManager_O(_lisp); // base(s): set(['llvmo::PassManagerBase_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__TargetSubtargetInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__TargetSubtargetInfo_O");
	Call_exposePython_llvmo__TargetSubtargetInfo_O(_lisp); // base(s): set(['llvmo::MCSubtargetInfo_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__User_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__User_O");
	Call_exposePython_llvmo__User_O(_lisp); // base(s): set(['llvmo::Value_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ValueAsMetadata_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ValueAsMetadata_O");
	Call_exposePython_llvmo__ValueAsMetadata_O(_lisp); // base(s): set(['llvmo::Metadata_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__NamedUnit_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__NamedUnit_O");
	Call_exposePython_units__NamedUnit_O(_lisp); // base(s): set(['units::Unit_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__BlockEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__BlockEnvironment_O");
	Call_exposePython_core__BlockEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__CatchEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__CatchEnvironment_O");
	Call_exposePython_core__CatchEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Class_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Class_O");
	Call_exposePython_core__Class_O(_lisp); // base(s): set(['core::Specializer_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__DoubleFloat_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__DoubleFloat_O");
	Call_exposePython_core__DoubleFloat_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__FunctionContainerEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__FunctionContainerEnvironment_O");
	Call_exposePython_core__FunctionContainerEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__FunctionValueEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__FunctionValueEnvironment_O");
	Call_exposePython_core__FunctionValueEnvironment_O(_lisp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__IOFileStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__IOFileStream_O");
	Call_exposePython_core__IOFileStream_O(_lisp); // base(s): set(['core::FileStream_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__IOStreamStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__IOStreamStream_O");
	Call_exposePython_core__IOStreamStream_O(_lisp); // base(s): set(['core::FileStream_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Integer_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Integer_O");
	Call_exposePython_core__Integer_O(_lisp); // base(s): set(['core::Rational_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__LongFloat_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__LongFloat_O");
	Call_exposePython_core__LongFloat_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__MacroletEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__MacroletEnvironment_O");
	Call_exposePython_core__MacroletEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Ratio_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Ratio_O");
	Call_exposePython_core__Ratio_O(_lisp); // base(s): set(['core::Rational_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__ShortFloat_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ShortFloat_O");
	Call_exposePython_core__ShortFloat_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__SimpleBitVector_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SimpleBitVector_O");
	Call_exposePython_core__SimpleBitVector_O(_lisp); // base(s): set(['core::BitVector_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__SingleFloat_dummy_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SingleFloat_dummy_O");
	Call_exposePython_core__SingleFloat_dummy_O(_lisp); // base(s): set(['core::Float_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__StackValueEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StackValueEnvironment_O");
	Call_exposePython_core__StackValueEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Str_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Str_O");
	Call_exposePython_core__Str_O(_lisp); // base(s): set(['core::String_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__StringInputStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StringInputStream_O");
	Call_exposePython_core__StringInputStream_O(_lisp); // base(s): set(['core::StringStream_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__StringOutputStream_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StringOutputStream_O");
	Call_exposePython_core__StringOutputStream_O(_lisp); // base(s): set(['core::StringStream_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__SymbolMacroletEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__SymbolMacroletEnvironment_O");
	Call_exposePython_core__SymbolMacroletEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__TagbodyEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__TagbodyEnvironment_O");
	Call_exposePython_core__TagbodyEnvironment_O(_lisp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__UnwindProtectEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__UnwindProtectEnvironment_O");
	Call_exposePython_core__UnwindProtectEnvironment_O(_lisp); // base(s): set(['core::CompileTimeEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ValueEnvironment_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ValueEnvironment_O");
	Call_exposePython_core__ValueEnvironment_O(_lisp); // base(s): set(['core::RuntimeVisibleEnvironment_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__VectorObjectsWithFillPtr_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__VectorObjectsWithFillPtr_O");
	Call_exposePython_core__VectorObjectsWithFillPtr_O(_lisp); // base(s): set(['core::VectorObjects_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__CoordinateArrayWithHash_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__CoordinateArrayWithHash_O");
	Call_exposePython_geom__CoordinateArrayWithHash_O(_lisp); // base(s): set(['geom::CoordinateArray_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Constant_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Constant_O");
	Call_exposePython_llvmo__Constant_O(_lisp); // base(s): set(['llvmo::User_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ImmutablePass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ImmutablePass_O");
	Call_exposePython_llvmo__ImmutablePass_O(_lisp); // base(s): set(['llvmo::ModulePass_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Instruction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Instruction_O");
	Call_exposePython_llvmo__Instruction_O(_lisp); // base(s): set(['llvmo::User_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__SequentialType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__SequentialType_O");
	Call_exposePython_llvmo__SequentialType_O(_lisp); // base(s): set(['llvmo::CompositeType_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__StructType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__StructType_O");
	Call_exposePython_llvmo__StructType_O(_lisp); // base(s): set(['llvmo::CompositeType_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Bignum_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Bignum_O");
	Call_exposePython_core__Bignum_O(_lisp); // base(s): set(['core::Integer_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__BuiltInClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__BuiltInClass_O");
	Call_exposePython_core__BuiltInClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__CxxClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__CxxClass_O");
	Call_exposePython_core__CxxClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__Fixnum_dummy_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__Fixnum_dummy_O");
	Call_exposePython_core__Fixnum_dummy_O(_lisp); // base(s): set(['core::Integer_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__ForwardReferencedClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__ForwardReferencedClass_O");
	Call_exposePython_core__ForwardReferencedClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClosPkg
extern void Call_exposePython_core__StdClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StdClass_O");
	Call_exposePython_core__StdClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // ifdef Use_ClosPkg
#ifdef Use_CorePkg
extern void Call_exposePython_core__StrWithFillPtr_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StrWithFillPtr_O");
	Call_exposePython_core__StrWithFillPtr_O(_lisp); // base(s): set(['core::Str_O'])
}
#endif // ifdef Use_CorePkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__StructureClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StructureClass_O");
	Call_exposePython_core__StructureClass_O(_lisp); // base(s): set(['core::Class_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ArrayType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ArrayType_O");
	Call_exposePython_llvmo__ArrayType_O(_lisp); // base(s): set(['llvmo::SequentialType_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__AtomicCmpXchgInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__AtomicCmpXchgInst_O");
	Call_exposePython_llvmo__AtomicCmpXchgInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__AtomicRMWInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__AtomicRMWInst_O");
	Call_exposePython_llvmo__AtomicRMWInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__BlockAddress_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__BlockAddress_O");
	Call_exposePython_llvmo__BlockAddress_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__CallInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__CallInst_O");
	Call_exposePython_llvmo__CallInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantArray_O");
	Call_exposePython_llvmo__ConstantArray_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantDataSequential_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantDataSequential_O");
	Call_exposePython_llvmo__ConstantDataSequential_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantExpr_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantExpr_O");
	Call_exposePython_llvmo__ConstantExpr_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantFP_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantFP_O");
	Call_exposePython_llvmo__ConstantFP_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantInt_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantInt_O");
	Call_exposePython_llvmo__ConstantInt_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantPointerNull_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantPointerNull_O");
	Call_exposePython_llvmo__ConstantPointerNull_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantStruct_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantStruct_O");
	Call_exposePython_llvmo__ConstantStruct_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__DataLayoutPass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__DataLayoutPass_O");
	Call_exposePython_llvmo__DataLayoutPass_O(_lisp); // base(s): set(['llvmo::ImmutablePass_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__FenceInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__FenceInst_O");
	Call_exposePython_llvmo__FenceInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__GlobalValue_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__GlobalValue_O");
	Call_exposePython_llvmo__GlobalValue_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__LandingPadInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__LandingPadInst_O");
	Call_exposePython_llvmo__LandingPadInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__PHINode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__PHINode_O");
	Call_exposePython_llvmo__PHINode_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__PointerType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__PointerType_O");
	Call_exposePython_llvmo__PointerType_O(_lisp); // base(s): set(['llvmo::SequentialType_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__StoreInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__StoreInst_O");
	Call_exposePython_llvmo__StoreInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__TargetLibraryInfo_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__TargetLibraryInfo_O");
	Call_exposePython_llvmo__TargetLibraryInfo_O(_lisp); // base(s): set(['llvmo::ImmutablePass_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__TerminatorInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__TerminatorInst_O");
	Call_exposePython_llvmo__TerminatorInst_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__UnaryInstruction_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__UnaryInstruction_O");
	Call_exposePython_llvmo__UnaryInstruction_O(_lisp); // base(s): set(['llvmo::Instruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__UndefValue_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__UndefValue_O");
	Call_exposePython_llvmo__UndefValue_O(_lisp); // base(s): set(['llvmo::Constant_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__VectorType_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__VectorType_O");
	Call_exposePython_llvmo__VectorType_O(_lisp); // base(s): set(['llvmo::SequentialType_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_ClbindPkg
extern void Call_exposePython_clbind__ClassRep_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: clbind__ClassRep_O");
	Call_exposePython_clbind__ClassRep_O(_lisp); // base(s): set(['core::BuiltInClass_O'])
}
#endif // ifdef Use_ClbindPkg
#ifdef Use_ClosPkg
extern void Call_exposePython_core__FuncallableStandardClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__FuncallableStandardClass_O");
	Call_exposePython_core__FuncallableStandardClass_O(_lisp); // base(s): set(['core::StdClass_O'])
}
#endif // ifdef Use_ClosPkg
#ifdef Use_ClPkg
extern void Call_exposePython_core__StandardClass_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: core__StandardClass_O");
	Call_exposePython_core__StandardClass_O(_lisp); // base(s): set(['core::StdClass_O'])
}
#endif // ifdef Use_ClPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__AllocaInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__AllocaInst_O");
	Call_exposePython_llvmo__AllocaInst_O(_lisp); // base(s): set(['llvmo::UnaryInstruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__BranchInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__BranchInst_O");
	Call_exposePython_llvmo__BranchInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ConstantDataArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ConstantDataArray_O");
	Call_exposePython_llvmo__ConstantDataArray_O(_lisp); // base(s): set(['llvmo::ConstantDataSequential_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__Function_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__Function_O");
	Call_exposePython_llvmo__Function_O(_lisp); // base(s): set(['llvmo::GlobalValue_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__GlobalVariable_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__GlobalVariable_O");
	Call_exposePython_llvmo__GlobalVariable_O(_lisp); // base(s): set(['llvmo::GlobalValue_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__IndirectBrInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__IndirectBrInst_O");
	Call_exposePython_llvmo__IndirectBrInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__InvokeInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__InvokeInst_O");
	Call_exposePython_llvmo__InvokeInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__LoadInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__LoadInst_O");
	Call_exposePython_llvmo__LoadInst_O(_lisp); // base(s): set(['llvmo::UnaryInstruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ResumeInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ResumeInst_O");
	Call_exposePython_llvmo__ResumeInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__ReturnInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__ReturnInst_O");
	Call_exposePython_llvmo__ReturnInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__SwitchInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__SwitchInst_O");
	Call_exposePython_llvmo__SwitchInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__UnreachableInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__UnreachableInst_O");
	Call_exposePython_llvmo__UnreachableInst_O(_lisp); // base(s): set(['llvmo::TerminatorInst_O'])
}
#endif // ifdef Use_LlvmoPkg
#ifdef Use_LlvmoPkg
extern void Call_exposePython_llvmo__VAArgInst_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: llvmo__VAArgInst_O");
	Call_exposePython_llvmo__VAArgInst_O(_lisp); // base(s): set(['llvmo::UnaryInstruction_O'])
}
#endif // ifdef Use_LlvmoPkg
#endif // EXPOSE_TO_PYTHON
#undef EXPOSE_TO_PYTHON
#if defined(EXPAND_CLASS_MACROS)
_CLASS_MACRO(adapt::ObjectSet_O)
_CLASS_MACRO(adapt::StringList_O)
_CLASS_MACRO(adapt::SymbolSet_O)
_CLASS_MACRO(core::T_O)
_CLASS_MACRO(asttooling::AstVisitor_O)
_CLASS_MACRO(cffi::Pointer_O)
_CLASS_MACRO(clbind::ClassRegistry_O)
_CLASS_MACRO(core::Archive_O)
_CLASS_MACRO(core::Array_O)
_CLASS_MACRO(core::Binder_O)
_CLASS_MACRO(core::CandoException_O)
_CLASS_MACRO(core::Character_dummy_O)
_CLASS_MACRO(core::Cons_O)
_CLASS_MACRO(core::CxxObject_O)
_CLASS_MACRO(core::DirectoryEntry_O)
_CLASS_MACRO(core::Environment_O)
_CLASS_MACRO(core::ExternalObject_O)
_CLASS_MACRO(core::FileStatus_O)
_CLASS_MACRO(core::Function_O)
_CLASS_MACRO(core::HashTable_O)
_CLASS_MACRO(core::IntArray_O)
_CLASS_MACRO(core::InvocationHistoryFrameIterator_O)
_CLASS_MACRO(core::Iterator_O)
_CLASS_MACRO(core::LambdaListHandler_O)
_CLASS_MACRO(core::LightUserData_O)
_CLASS_MACRO(core::LoadTimeValues_O)
_CLASS_MACRO(core::MultiStringBuffer_O)
_CLASS_MACRO(core::Number_O)
_CLASS_MACRO(core::Package_O)
_CLASS_MACRO(core::Path_O)
_CLASS_MACRO(core::Pathname_O)
_CLASS_MACRO(core::Pointer_O)
_CLASS_MACRO(core::PosixTimeDuration_O)
_CLASS_MACRO(core::PosixTime_O)
_CLASS_MACRO(core::RandomState_O)
_CLASS_MACRO(core::ReadTable_O)
_CLASS_MACRO(core::Reader_O)
_CLASS_MACRO(core::Record_O)
_CLASS_MACRO(core::RegexMatch_O)
_CLASS_MACRO(core::Regex_O)
_CLASS_MACRO(core::SNode_O)
_CLASS_MACRO(core::SingleDispatchMethod_O)
_CLASS_MACRO(core::SmallMap_O)
_CLASS_MACRO(core::SmallMultimap_O)
_CLASS_MACRO(core::SourceFileInfo_O)
_CLASS_MACRO(core::SourceManager_O)
_CLASS_MACRO(core::SourcePosInfo_O)
_CLASS_MACRO(core::StandardObject_O)
_CLASS_MACRO(core::Stream_O)
_CLASS_MACRO(core::StructureObject_O)
_CLASS_MACRO(core::SymbolToEnumConverter_O)
_CLASS_MACRO(core::Symbol_O)
_CLASS_MACRO(core::VaList_dummy_O)
_CLASS_MACRO(core::WeakHashTable_O)
_CLASS_MACRO(core::WeakKeyMapping_O)
_CLASS_MACRO(core::WeakPointer_O)
_CLASS_MACRO(core::WrappedPointer_O)
_CLASS_MACRO(llvmo::AttributeSet_O)
_CLASS_MACRO(llvmo::Attribute_O)
_CLASS_MACRO(llvmo::DebugInfo_O)
_CLASS_MACRO(llvmo::DebugLoc_O)
_CLASS_MACRO(llvmo::InsertPoint_O)
_CLASS_MACRO(adapt::IterateCons_O)
_CLASS_MACRO(adapt::QDomNode_O)
_CLASS_MACRO(adapt::StringSet_O)
_CLASS_MACRO(adapt::SymbolList_O)
_CLASS_MACRO(core::ActivationFrame_O)
_CLASS_MACRO(core::AnsiStream_O)
_CLASS_MACRO(core::ArrayDisplaced_O)
_CLASS_MACRO(core::ArrayObjects_O)
_CLASS_MACRO(core::BranchSNode_O)
_CLASS_MACRO(core::CompiledFunction_O)
_CLASS_MACRO(core::Complex_O)
_CLASS_MACRO(core::DirectoryIterator_O)
_CLASS_MACRO(core::ForeignData_O)
_CLASS_MACRO(core::GlueEnvironment_O)
_CLASS_MACRO(core::HashTableEq_O)
_CLASS_MACRO(core::HashTableEql_O)
_CLASS_MACRO(core::HashTableEqual_O)
_CLASS_MACRO(core::HashTableEqualp_O)
_CLASS_MACRO(core::Instance_O)
_CLASS_MACRO(core::LeafSNode_O)
_CLASS_MACRO(core::LexicalEnvironment_O)
_CLASS_MACRO(core::LoadArchive_O)
_CLASS_MACRO(core::LogicalPathname_O)
_CLASS_MACRO(core::Metaobject_O)
_CLASS_MACRO(core::Null_O)
_CLASS_MACRO(core::Real_O)
_CLASS_MACRO(core::RecursiveDirectoryIterator_O)
_CLASS_MACRO(core::SaveArchive_O)
_CLASS_MACRO(core::SingleDispatchEffectiveMethodFunction_O)
_CLASS_MACRO(core::SingleDispatchGenericFunction_O)
_CLASS_MACRO(core::SpecialForm_O)
_CLASS_MACRO(core::UserData_O)
_CLASS_MACRO(core::Vector_O)
_CLASS_MACRO(core::WeakKeyHashTable_O)
_CLASS_MACRO(geom::BoundingBox_O)
_CLASS_MACRO(geom::Color_O)
_CLASS_MACRO(geom::OMatrix_O)
_CLASS_MACRO(geom::OVector2_O)
_CLASS_MACRO(geom::OVector3_O)
_CLASS_MACRO(llvmo::APFloat_O)
_CLASS_MACRO(llvmo::APInt_O)
_CLASS_MACRO(llvmo::DIArray_O)
_CLASS_MACRO(llvmo::DIBasicType_O)
_CLASS_MACRO(llvmo::DIBuilder_O)
_CLASS_MACRO(llvmo::DICompileUnit_O)
_CLASS_MACRO(llvmo::DICompositeType_O)
_CLASS_MACRO(llvmo::DIDerivedType_O)
_CLASS_MACRO(llvmo::DIDescriptor_O)
_CLASS_MACRO(llvmo::DIFile_O)
_CLASS_MACRO(llvmo::DILexicalBlock_O)
_CLASS_MACRO(llvmo::DIScope_O)
_CLASS_MACRO(llvmo::DISubprogram_O)
_CLASS_MACRO(llvmo::DISubroutineType_O)
_CLASS_MACRO(llvmo::DITypeArray_O)
_CLASS_MACRO(llvmo::DIType_O)
_CLASS_MACRO(llvmo::DataLayout_O)
_CLASS_MACRO(llvmo::EngineBuilder_O)
_CLASS_MACRO(llvmo::ExecutionEngine_O)
_CLASS_MACRO(llvmo::IRBuilderBase_O)
_CLASS_MACRO(llvmo::LLVMContext_O)
_CLASS_MACRO(llvmo::Linker_O)
_CLASS_MACRO(llvmo::MCSubtargetInfo_O)
_CLASS_MACRO(llvmo::Metadata_O)
_CLASS_MACRO(llvmo::Module_O)
_CLASS_MACRO(llvmo::NamedMDNode_O)
_CLASS_MACRO(llvmo::PassManagerBase_O)
_CLASS_MACRO(llvmo::PassManagerBuilder_O)
_CLASS_MACRO(llvmo::Pass_O)
_CLASS_MACRO(llvmo::TargetMachine_O)
_CLASS_MACRO(llvmo::TargetOptions_O)
_CLASS_MACRO(llvmo::Target_O)
_CLASS_MACRO(llvmo::Triple_O)
_CLASS_MACRO(llvmo::Type_O)
_CLASS_MACRO(llvmo::Value_O)
_CLASS_MACRO(units::Dimension_O)
_CLASS_MACRO(units::Quantity_O)
_CLASS_MACRO(units::Unit_O)
_CLASS_MACRO(core::BitVector_O)
_CLASS_MACRO(core::BroadcastStream_O)
_CLASS_MACRO(core::CompileTimeEnvironment_O)
_CLASS_MACRO(core::ConcatenatedStream_O)
_CLASS_MACRO(core::EchoStream_O)
_CLASS_MACRO(core::FileStream_O)
_CLASS_MACRO(core::Float_O)
_CLASS_MACRO(core::FunctionFrame_O)
_CLASS_MACRO(core::Rational_O)
_CLASS_MACRO(core::RuntimeVisibleEnvironment_O)
_CLASS_MACRO(core::SexpLoadArchive_O)
_CLASS_MACRO(core::SexpSaveArchive_O)
_CLASS_MACRO(core::Specializer_O)
_CLASS_MACRO(core::StringStream_O)
_CLASS_MACRO(core::String_O)
_CLASS_MACRO(core::SynonymStream_O)
_CLASS_MACRO(core::TagbodyFrame_O)
_CLASS_MACRO(core::TwoWayStream_O)
_CLASS_MACRO(core::ValueFrame_O)
_CLASS_MACRO(core::VectorObjects_O)
_CLASS_MACRO(geom::CoordinateArray_O)
_CLASS_MACRO(llvmo::Argument_O)
_CLASS_MACRO(llvmo::BasicBlock_O)
_CLASS_MACRO(llvmo::CompositeType_O)
_CLASS_MACRO(llvmo::FunctionPassManager_O)
_CLASS_MACRO(llvmo::FunctionPass_O)
_CLASS_MACRO(llvmo::FunctionType_O)
_CLASS_MACRO(llvmo::IRBuilder_O)
_CLASS_MACRO(llvmo::IntegerType_O)
_CLASS_MACRO(llvmo::LLVMTargetMachine_O)
_CLASS_MACRO(llvmo::MDNode_O)
_CLASS_MACRO(llvmo::MDString_O)
_CLASS_MACRO(llvmo::ModulePass_O)
_CLASS_MACRO(llvmo::PassManager_O)
_CLASS_MACRO(llvmo::TargetSubtargetInfo_O)
_CLASS_MACRO(llvmo::User_O)
_CLASS_MACRO(llvmo::ValueAsMetadata_O)
_CLASS_MACRO(units::NamedUnit_O)
_CLASS_MACRO(core::BlockEnvironment_O)
_CLASS_MACRO(core::CatchEnvironment_O)
_CLASS_MACRO(core::Class_O)
_CLASS_MACRO(core::DoubleFloat_O)
_CLASS_MACRO(core::FunctionContainerEnvironment_O)
_CLASS_MACRO(core::FunctionValueEnvironment_O)
_CLASS_MACRO(core::IOFileStream_O)
_CLASS_MACRO(core::IOStreamStream_O)
_CLASS_MACRO(core::Integer_O)
_CLASS_MACRO(core::LongFloat_O)
_CLASS_MACRO(core::MacroletEnvironment_O)
_CLASS_MACRO(core::Ratio_O)
_CLASS_MACRO(core::ShortFloat_O)
_CLASS_MACRO(core::SimpleBitVector_O)
_CLASS_MACRO(core::SingleFloat_dummy_O)
_CLASS_MACRO(core::StackValueEnvironment_O)
_CLASS_MACRO(core::Str_O)
_CLASS_MACRO(core::StringInputStream_O)
_CLASS_MACRO(core::StringOutputStream_O)
_CLASS_MACRO(core::SymbolMacroletEnvironment_O)
_CLASS_MACRO(core::TagbodyEnvironment_O)
_CLASS_MACRO(core::UnwindProtectEnvironment_O)
_CLASS_MACRO(core::ValueEnvironment_O)
_CLASS_MACRO(core::VectorObjectsWithFillPtr_O)
_CLASS_MACRO(geom::CoordinateArrayWithHash_O)
_CLASS_MACRO(llvmo::Constant_O)
_CLASS_MACRO(llvmo::ImmutablePass_O)
_CLASS_MACRO(llvmo::Instruction_O)
_CLASS_MACRO(llvmo::SequentialType_O)
_CLASS_MACRO(llvmo::StructType_O)
_CLASS_MACRO(core::Bignum_O)
_CLASS_MACRO(core::BuiltInClass_O)
_CLASS_MACRO(core::CxxClass_O)
_CLASS_MACRO(core::Fixnum_dummy_O)
_CLASS_MACRO(core::ForwardReferencedClass_O)
_CLASS_MACRO(core::StdClass_O)
_CLASS_MACRO(core::StrWithFillPtr_O)
_CLASS_MACRO(core::StructureClass_O)
_CLASS_MACRO(llvmo::ArrayType_O)
_CLASS_MACRO(llvmo::AtomicCmpXchgInst_O)
_CLASS_MACRO(llvmo::AtomicRMWInst_O)
_CLASS_MACRO(llvmo::BlockAddress_O)
_CLASS_MACRO(llvmo::CallInst_O)
_CLASS_MACRO(llvmo::ConstantArray_O)
_CLASS_MACRO(llvmo::ConstantDataSequential_O)
_CLASS_MACRO(llvmo::ConstantExpr_O)
_CLASS_MACRO(llvmo::ConstantFP_O)
_CLASS_MACRO(llvmo::ConstantInt_O)
_CLASS_MACRO(llvmo::ConstantPointerNull_O)
_CLASS_MACRO(llvmo::ConstantStruct_O)
_CLASS_MACRO(llvmo::DataLayoutPass_O)
_CLASS_MACRO(llvmo::FenceInst_O)
_CLASS_MACRO(llvmo::GlobalValue_O)
_CLASS_MACRO(llvmo::LandingPadInst_O)
_CLASS_MACRO(llvmo::PHINode_O)
_CLASS_MACRO(llvmo::PointerType_O)
_CLASS_MACRO(llvmo::StoreInst_O)
_CLASS_MACRO(llvmo::TargetLibraryInfo_O)
_CLASS_MACRO(llvmo::TerminatorInst_O)
_CLASS_MACRO(llvmo::UnaryInstruction_O)
_CLASS_MACRO(llvmo::UndefValue_O)
_CLASS_MACRO(llvmo::VectorType_O)
_CLASS_MACRO(clbind::ClassRep_O)
_CLASS_MACRO(core::FuncallableStandardClass_O)
_CLASS_MACRO(core::StandardClass_O)
_CLASS_MACRO(llvmo::AllocaInst_O)
_CLASS_MACRO(llvmo::BranchInst_O)
_CLASS_MACRO(llvmo::ConstantDataArray_O)
_CLASS_MACRO(llvmo::Function_O)
_CLASS_MACRO(llvmo::GlobalVariable_O)
_CLASS_MACRO(llvmo::IndirectBrInst_O)
_CLASS_MACRO(llvmo::InvokeInst_O)
_CLASS_MACRO(llvmo::LoadInst_O)
_CLASS_MACRO(llvmo::ResumeInst_O)
_CLASS_MACRO(llvmo::ReturnInst_O)
_CLASS_MACRO(llvmo::SwitchInst_O)
_CLASS_MACRO(llvmo::UnreachableInst_O)
_CLASS_MACRO(llvmo::VAArgInst_O)
#endif // EXPAND_CLASS_MACROS
#undef ALL_STAGES
// ---------------- after class initializers
#ifdef EXPOSE_TO_PYTHON
#endif // EXPOSE_TO_PYTHON
