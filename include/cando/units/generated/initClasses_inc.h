// start
// define cpp macros: SET_SYMBOL, CREATE_CLASS, SET_CLASS, DEFINE_BASE_CLASSES, DEFINE_CLASS_NAMES, EXPOSE_TO_CANDO 
// define cpp macro: ALL_STAGES to get the effect of defining all of the macros above
// define cpp macro: EXPOSE_PYTHON to expose python
// Associating namespace(units) with package(UnitsPkg)
 // class Dimension_O : public core::T_O
 // class NamedUnit_O : public Unit_O
 // class Quantity_O : public core::T_O
 // class Unit_O : public core::T_O
// Associating namespace(units) with package(UnitsPkg)
#ifdef HEADER_INCLUDES
#include "include/dimension.h"
#include "include/quantity.h"
#include "include/unit.h"
#include "include/namedUnit.h"
#endif // HEADER_INCLUDES
#undef HEADER_INCLUDES
#if defined(SET_SYMBOL) || defined(ALL_STAGES)
// requires LOOKUP_SYMBOL(pkg,symbolName) be defined
units::Dimension_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::Dimension_O::static_packageName(),units::Dimension_O::static_className()));
units::Quantity_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::Quantity_O::static_packageName(),units::Quantity_O::static_className()));
units::Unit_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::Unit_O::static_packageName(),units::Unit_O::static_className()));
units::NamedUnit_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(units::NamedUnit_O::static_packageName(),units::NamedUnit_O::static_className()));
#endif // SET_SYMBOL
#undef SET_SYMBOL
#if defined(CREATE_CLASS) || defined(ALL_STAGES)

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
        core::LispObjectCreator<units::Dimension_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<units::Dimension_O>>::allocateClass();
        units::Dimension_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::Dimension_O::static_className() % (void*)(units::Dimension_O::static_allocator) );
    classunits__Dimension_Oval->setCreator(units::Dimension_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::Dimension_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__Dimension_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<units::Quantity_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<units::Quantity_O>>::allocateClass();
        units::Quantity_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::Quantity_O::static_className() % (void*)(units::Quantity_O::static_allocator) );
    classunits__Quantity_Oval->setCreator(units::Quantity_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::Quantity_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__Quantity_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<units::Unit_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<units::Unit_O>>::allocateClass();
        units::Unit_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::Unit_O::static_className() % (void*)(units::Unit_O::static_allocator) );
    classunits__Unit_Oval->setCreator(units::Unit_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::Unit_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__Unit_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<units::NamedUnit_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<units::NamedUnit_O>>::allocateClass();
        units::NamedUnit_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% units::NamedUnit_O::static_className() % (void*)(units::NamedUnit_O::static_allocator) );
    classunits__NamedUnit_Oval->setCreator(units::NamedUnit_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % units::NamedUnit_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classunits__NamedUnit_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif // CREATE_CLASS
#undef CREATE_CLASS
#ifdef DUMP_INFO_CLASS // {
// Depends on nothing

    LOG(BF("---    dump_info   --- className: units::Dimension_O @ %X") % classunits__Dimension_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::Dimension_O::static_className() % units::Dimension_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: units::Quantity_O @ %X") % classunits__Quantity_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::Quantity_O::static_className() % units::Quantity_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: units::Unit_O @ %X") % classunits__Unit_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::Unit_O::static_className() % units::Unit_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: units::NamedUnit_O @ %X") % classunits__NamedUnit_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % units::NamedUnit_O::static_className() % units::NamedUnit_O::static_classSymbol() );
#endif // } DUMP_INFO_CLASS
#undef DUMP_INFO_CLASS
#if defined(DEFINE_BASE_CLASSES) || defined(ALL_STAGES) // {
// Depends on nothing
classunits__Dimension_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classunits__Quantity_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classunits__Unit_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classunits__NamedUnit_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(units::Unit_O::static_classSymbol());
#endif // } DEFINE_BASE_CLASSES
#undef DEFINE_BASE_CLASSES
#if defined(DEFINE_CLASS_NAMES) || defined(ALL_STAGES) // {
// Depends on nothing

    classunits__Dimension_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::Dimension_O::static_classSymbol());

    classunits__Quantity_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::Quantity_O::static_classSymbol());

    classunits__Unit_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::Unit_O::static_classSymbol());

    classunits__NamedUnit_Oval->__setupStage3NameAndCalculateClassPrecedenceList(units::NamedUnit_O::static_classSymbol());
#endif // } DEFINE_CLASS_NAMES
#undef DEFINE_CLASS_NAMES
#if defined(EXPOSE_TO_CANDO) || defined(ALL_STAGES)
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__Dimension_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__Dimension_O");
    units::Register_units__Dimension_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__Quantity_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__Quantity_O");
    units::Register_units__Quantity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
#ifdef EXTERN_REGISTER
extern void Register_units__Unit_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_units__Unit_O");
    units::Register_units__Unit_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_UnitsPkg
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
#endif // EXPOSE_TO_CANDO
#undef EXPOSE_TO_CANDO
#ifdef EXPOSE_TO_PYTHON
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__Dimension_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__Dimension_O");
	Call_exposePython_units__Dimension_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__Quantity_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__Quantity_O");
	Call_exposePython_units__Quantity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__Unit_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__Unit_O");
	Call_exposePython_units__Unit_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Call_exposePython_units__NamedUnit_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: units__NamedUnit_O");
	Call_exposePython_units__NamedUnit_O(_lisp); // base(s): set(['units::Unit_O'])
}
#endif // ifdef Use_UnitsPkg
#endif // EXPOSE_TO_PYTHON
#undef EXPOSE_TO_PYTHON
#if defined(EXPAND_CLASS_MACROS)
_CLASS_MACRO(units::Dimension_O)
_CLASS_MACRO(units::Quantity_O)
_CLASS_MACRO(units::Unit_O)
_CLASS_MACRO(units::NamedUnit_O)
#endif // EXPAND_CLASS_MACROS
#undef ALL_STAGES
// ---------------- after class initializers
#ifdef EXPOSE_TO_PYTHON
#endif // EXPOSE_TO_PYTHON
