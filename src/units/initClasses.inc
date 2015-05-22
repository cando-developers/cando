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
#include "dimension.h"
#include "quantity.h"
#include "unit.h"
#include "namedUnit.h"
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
// Depends on undefinedMetaClass
    core::MetaClass_sp undefinedMetaClass; undefinedMetaClass.reset();

    LOG(BF("Creating class[classunits__Dimension_Oval]"));
    core::BuiltInClass_sp classunits__Dimension_Oval = core::BuiltInClass_O::create(undefinedMetaClass);
    classunits__Dimension_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__Dimension_Oval,_lisp,units::Dimension_O::static_classSymbol());
    units::Dimension_O::___staticMetaClass = classunits__Dimension_Oval;
    _lisp->setf_findClass(units::Dimension_O::static_classSymbol(),classunits__Dimension_Oval);
    {
        AllocatorCallback cb = &new_Nil<units::Dimension_O>;
        units::Dimension_O::___set_static_newNil_callback(cb);
    }
    LOG(BF("Set _newNil_callback for class(%s) to %X")% units::Dimension_O::static_className() % (void*)(units::Dimension_O::static_newNil_callback()) );
    classunits__Dimension_Oval->setInstance_newNil_callback(units::Dimension_O::static_newNil_callback());
    {
        boost::shared_ptr<units::Dimension_O> nil_for_class(new units::Dimension_O(undefinedMetaClass));
        nil_for_class->__setWeakThis(nil_for_class);
        LOG(BF("Created nil for class[%s]") % units::Dimension_O::static_className() );
        units::Dimension_O::_nil = nil_for_class;
        classunits__Dimension_Oval->setInstanceNil(nil_for_class);
    }
    classunits__Dimension_Oval->setSupportsSlots(units::Dimension_O::static_supportsSlots());
    /* ----- the class and its nil are now defined but the _class of these objects is undefined - set in next stage  ----- */

    LOG(BF("Creating class[classunits__Quantity_Oval]"));
    core::BuiltInClass_sp classunits__Quantity_Oval = core::BuiltInClass_O::create(undefinedMetaClass);
    classunits__Quantity_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__Quantity_Oval,_lisp,units::Quantity_O::static_classSymbol());
    units::Quantity_O::___staticMetaClass = classunits__Quantity_Oval;
    _lisp->setf_findClass(units::Quantity_O::static_classSymbol(),classunits__Quantity_Oval);
    {
        AllocatorCallback cb = &new_Nil<units::Quantity_O>;
        units::Quantity_O::___set_static_newNil_callback(cb);
    }
    LOG(BF("Set _newNil_callback for class(%s) to %X")% units::Quantity_O::static_className() % (void*)(units::Quantity_O::static_newNil_callback()) );
    classunits__Quantity_Oval->setInstance_newNil_callback(units::Quantity_O::static_newNil_callback());
    {
        boost::shared_ptr<units::Quantity_O> nil_for_class(new units::Quantity_O(undefinedMetaClass));
        nil_for_class->__setWeakThis(nil_for_class);
        LOG(BF("Created nil for class[%s]") % units::Quantity_O::static_className() );
        units::Quantity_O::_nil = nil_for_class;
        classunits__Quantity_Oval->setInstanceNil(nil_for_class);
    }
    classunits__Quantity_Oval->setSupportsSlots(units::Quantity_O::static_supportsSlots());
    /* ----- the class and its nil are now defined but the _class of these objects is undefined - set in next stage  ----- */

    LOG(BF("Creating class[classunits__Unit_Oval]"));
    core::BuiltInClass_sp classunits__Unit_Oval = core::BuiltInClass_O::create(undefinedMetaClass);
    classunits__Unit_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__Unit_Oval,_lisp,units::Unit_O::static_classSymbol());
    units::Unit_O::___staticMetaClass = classunits__Unit_Oval;
    _lisp->setf_findClass(units::Unit_O::static_classSymbol(),classunits__Unit_Oval);
    {
        AllocatorCallback cb = &new_Nil<units::Unit_O>;
        units::Unit_O::___set_static_newNil_callback(cb);
    }
    LOG(BF("Set _newNil_callback for class(%s) to %X")% units::Unit_O::static_className() % (void*)(units::Unit_O::static_newNil_callback()) );
    classunits__Unit_Oval->setInstance_newNil_callback(units::Unit_O::static_newNil_callback());
    {
        boost::shared_ptr<units::Unit_O> nil_for_class(new units::Unit_O(undefinedMetaClass));
        nil_for_class->__setWeakThis(nil_for_class);
        LOG(BF("Created nil for class[%s]") % units::Unit_O::static_className() );
        units::Unit_O::_nil = nil_for_class;
        classunits__Unit_Oval->setInstanceNil(nil_for_class);
    }
    classunits__Unit_Oval->setSupportsSlots(units::Unit_O::static_supportsSlots());
    /* ----- the class and its nil are now defined but the _class of these objects is undefined - set in next stage  ----- */

    LOG(BF("Creating class[classunits__NamedUnit_Oval]"));
    core::BuiltInClass_sp classunits__NamedUnit_Oval = core::BuiltInClass_O::create(undefinedMetaClass);
    classunits__NamedUnit_Oval->__setup_stage1_with_sharedPtr_lisp_sid(classunits__NamedUnit_Oval,_lisp,units::NamedUnit_O::static_classSymbol());
    units::NamedUnit_O::___staticMetaClass = classunits__NamedUnit_Oval;
    _lisp->setf_findClass(units::NamedUnit_O::static_classSymbol(),classunits__NamedUnit_Oval);
    {
        AllocatorCallback cb = &new_Nil<units::NamedUnit_O>;
        units::NamedUnit_O::___set_static_newNil_callback(cb);
    }
    LOG(BF("Set _newNil_callback for class(%s) to %X")% units::NamedUnit_O::static_className() % (void*)(units::NamedUnit_O::static_newNil_callback()) );
    classunits__NamedUnit_Oval->setInstance_newNil_callback(units::NamedUnit_O::static_newNil_callback());
    {
        boost::shared_ptr<units::NamedUnit_O> nil_for_class(new units::NamedUnit_O(undefinedMetaClass));
        nil_for_class->__setWeakThis(nil_for_class);
        LOG(BF("Created nil for class[%s]") % units::NamedUnit_O::static_className() );
        units::NamedUnit_O::_nil = nil_for_class;
        classunits__NamedUnit_Oval->setInstanceNil(nil_for_class);
    }
    classunits__NamedUnit_Oval->setSupportsSlots(units::NamedUnit_O::static_supportsSlots());
    /* ----- the class and its nil are now defined but the _class of these objects is undefined - set in next stage  ----- */
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
classunits__Dimension_Oval->addInstanceBaseClass(core::T_O::static_classSymbol());
classunits__Quantity_Oval->addInstanceBaseClass(core::T_O::static_classSymbol());
classunits__Unit_Oval->addInstanceBaseClass(core::T_O::static_classSymbol());
classunits__NamedUnit_Oval->addInstanceBaseClass(units::Unit_O::static_classSymbol());
#endif // } DEFINE_BASE_CLASSES
#undef DEFINE_BASE_CLASSES
#if defined(DEFINE_CLASS_NAMES) || defined(ALL_STAGES) // {
 core::Package_sp _curPkg = _lisp->findPackage(CurrentPkg);
// Depends on nothing

    classunits__Dimension_Oval->__setup_stage3_name(units::Dimension_O::static_classSymbol());

    classunits__Quantity_Oval->__setup_stage3_name(units::Quantity_O::static_classSymbol());

    classunits__Unit_Oval->__setup_stage3_name(units::Unit_O::static_classSymbol());

    classunits__NamedUnit_Oval->__setup_stage3_name(units::NamedUnit_O::static_classSymbol());
#endif // } DEFINE_CLASS_NAMES
#undef DEFINE_CLASS_NAMES
#if defined(EXPOSE_TO_CANDO) || defined(ALL_STAGES)
#ifdef Use_UnitsPkg
extern void Register_units__Dimension_O(core::Lisp_sp); // base(s): set(['core::T_O'])
{_BLOCK_TRACE("initializing Register_units__Dimension_O");
    Register_units__Dimension_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Register_units__Quantity_O(core::Lisp_sp); // base(s): set(['core::T_O'])
{_BLOCK_TRACE("initializing Register_units__Quantity_O");
    Register_units__Quantity_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Register_units__Unit_O(core::Lisp_sp); // base(s): set(['core::T_O'])
{_BLOCK_TRACE("initializing Register_units__Unit_O");
    Register_units__Unit_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_UnitsPkg
#ifdef Use_UnitsPkg
extern void Register_units__NamedUnit_O(core::Lisp_sp); // base(s): set(['units::Unit_O'])
{_BLOCK_TRACE("initializing Register_units__NamedUnit_O");
    Register_units__NamedUnit_O(_lisp); // base(s): set(['units::Unit_O'])
}
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
#undef ALL_STAGES




