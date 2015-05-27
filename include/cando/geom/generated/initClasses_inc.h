// start
// define cpp macros: SET_SYMBOL, CREATE_CLASS, SET_CLASS, DEFINE_BASE_CLASSES, DEFINE_CLASS_NAMES, EXPOSE_TO_CANDO 
// define cpp macro: ALL_STAGES to get the effect of defining all of the macros above
// define cpp macro: EXPOSE_PYTHON to expose python
 // class BoundingBox_O : public core::T_O
 // class Color_O : public core::T_O
 // class CoordinateArray_O : public core::Vector_O
 // class CoordinateArrayWithHash_O : public CoordinateArray_O
// Associating namespace(geom) with package(GeomPkg)
 // class OMatrix_O : public core::T_O
 // class OVector2_O : public core::T_O
 // class OVector3_O : public core::T_O
#ifdef HEADER_INCLUDES
#include "include/boundingBox.h"
#include "include/color.h"
#include "include/coordinateArray.h"
#include "include/omatrix.h"
#include "include/ovector2.h"
#include "include/ovector3.h"
#endif // HEADER_INCLUDES
#undef HEADER_INCLUDES
#if defined(SET_SYMBOL) || defined(ALL_STAGES)
// requires LOOKUP_SYMBOL(pkg,symbolName) be defined
geom::BoundingBox_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::BoundingBox_O::static_packageName(),geom::BoundingBox_O::static_className()));
geom::Color_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::Color_O::static_packageName(),geom::Color_O::static_className()));
geom::CoordinateArray_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::CoordinateArray_O::static_packageName(),geom::CoordinateArray_O::static_className()));
geom::OMatrix_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::OMatrix_O::static_packageName(),geom::OMatrix_O::static_className()));
geom::OVector2_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::OVector2_O::static_packageName(),geom::OVector2_O::static_className()));
geom::OVector3_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::OVector3_O::static_packageName(),geom::OVector3_O::static_className()));
geom::CoordinateArrayWithHash_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(geom::CoordinateArrayWithHash_O::static_packageName(),geom::CoordinateArrayWithHash_O::static_className()));
#endif // SET_SYMBOL
#undef SET_SYMBOL
#if defined(CREATE_CLASS) || defined(ALL_STAGES)

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
        core::LispObjectCreator<geom::BoundingBox_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::BoundingBox_O>>::allocateClass();
        geom::BoundingBox_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::BoundingBox_O::static_className() % (void*)(geom::BoundingBox_O::static_allocator) );
    classgeom__BoundingBox_Oval->setCreator(geom::BoundingBox_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::BoundingBox_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__BoundingBox_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<geom::Color_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::Color_O>>::allocateClass();
        geom::Color_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::Color_O::static_className() % (void*)(geom::Color_O::static_allocator) );
    classgeom__Color_Oval->setCreator(geom::Color_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::Color_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__Color_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<geom::CoordinateArray_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::CoordinateArray_O>>::allocateClass();
        geom::CoordinateArray_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::CoordinateArray_O::static_className() % (void*)(geom::CoordinateArray_O::static_allocator) );
    classgeom__CoordinateArray_Oval->setCreator(geom::CoordinateArray_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::CoordinateArray_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__CoordinateArray_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<geom::OMatrix_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::OMatrix_O>>::allocateClass();
        geom::OMatrix_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::OMatrix_O::static_className() % (void*)(geom::OMatrix_O::static_allocator) );
    classgeom__OMatrix_Oval->setCreator(geom::OMatrix_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::OMatrix_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__OMatrix_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<geom::OVector2_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::OVector2_O>>::allocateClass();
        geom::OVector2_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::OVector2_O::static_className() % (void*)(geom::OVector2_O::static_allocator) );
    classgeom__OVector2_Oval->setCreator(geom::OVector2_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::OVector2_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__OVector2_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<geom::OVector3_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::OVector3_O>>::allocateClass();
        geom::OVector3_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::OVector3_O::static_className() % (void*)(geom::OVector3_O::static_allocator) );
    classgeom__OVector3_Oval->setCreator(geom::OVector3_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::OVector3_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__OVector3_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<geom::CoordinateArrayWithHash_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<geom::CoordinateArrayWithHash_O>>::allocateClass();
        geom::CoordinateArrayWithHash_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% geom::CoordinateArrayWithHash_O::static_className() % (void*)(geom::CoordinateArrayWithHash_O::static_allocator) );
    classgeom__CoordinateArrayWithHash_Oval->setCreator(geom::CoordinateArrayWithHash_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % geom::CoordinateArrayWithHash_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classgeom__CoordinateArrayWithHash_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif // CREATE_CLASS
#undef CREATE_CLASS
#ifdef DUMP_INFO_CLASS // {
// Depends on nothing

    LOG(BF("---    dump_info   --- className: geom::BoundingBox_O @ %X") % classgeom__BoundingBox_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::BoundingBox_O::static_className() % geom::BoundingBox_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: geom::Color_O @ %X") % classgeom__Color_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::Color_O::static_className() % geom::Color_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: geom::CoordinateArray_O @ %X") % classgeom__CoordinateArray_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::CoordinateArray_O::static_className() % geom::CoordinateArray_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: geom::OMatrix_O @ %X") % classgeom__OMatrix_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::OMatrix_O::static_className() % geom::OMatrix_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: geom::OVector2_O @ %X") % classgeom__OVector2_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::OVector2_O::static_className() % geom::OVector2_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: geom::OVector3_O @ %X") % classgeom__OVector3_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::OVector3_O::static_className() % geom::OVector3_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: geom::CoordinateArrayWithHash_O @ %X") % classgeom__CoordinateArrayWithHash_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % geom::CoordinateArrayWithHash_O::static_className() % geom::CoordinateArrayWithHash_O::static_classSymbol() );
#endif // } DUMP_INFO_CLASS
#undef DUMP_INFO_CLASS
#if defined(DEFINE_BASE_CLASSES) || defined(ALL_STAGES) // {
// Depends on nothing
classgeom__BoundingBox_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classgeom__Color_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classgeom__CoordinateArray_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Vector_O::static_classSymbol());
classgeom__OMatrix_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classgeom__OVector2_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classgeom__OVector3_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classgeom__CoordinateArrayWithHash_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(geom::CoordinateArray_O::static_classSymbol());
#endif // } DEFINE_BASE_CLASSES
#undef DEFINE_BASE_CLASSES
#if defined(DEFINE_CLASS_NAMES) || defined(ALL_STAGES) // {
// Depends on nothing

    classgeom__BoundingBox_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::BoundingBox_O::static_classSymbol());

    classgeom__Color_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::Color_O::static_classSymbol());

    classgeom__CoordinateArray_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::CoordinateArray_O::static_classSymbol());

    classgeom__OMatrix_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::OMatrix_O::static_classSymbol());

    classgeom__OVector2_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::OVector2_O::static_classSymbol());

    classgeom__OVector3_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::OVector3_O::static_classSymbol());

    classgeom__CoordinateArrayWithHash_Oval->__setupStage3NameAndCalculateClassPrecedenceList(geom::CoordinateArrayWithHash_O::static_classSymbol());
#endif // } DEFINE_CLASS_NAMES
#undef DEFINE_CLASS_NAMES
#if defined(EXPOSE_TO_CANDO) || defined(ALL_STAGES)
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__BoundingBox_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__BoundingBox_O");
    geom::Register_geom__BoundingBox_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__Color_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__Color_O");
    geom::Register_geom__Color_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
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
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__OMatrix_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__OMatrix_O");
    geom::Register_geom__OMatrix_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__OVector2_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__OVector2_O");
    geom::Register_geom__OVector2_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
#ifdef EXTERN_REGISTER
extern void Register_geom__OVector3_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_geom__OVector3_O");
    geom::Register_geom__OVector3_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_GeomPkg
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
#endif // EXPOSE_TO_CANDO
#undef EXPOSE_TO_CANDO
#ifdef EXPOSE_TO_PYTHON
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__BoundingBox_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__BoundingBox_O");
	Call_exposePython_geom__BoundingBox_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__Color_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__Color_O");
	Call_exposePython_geom__Color_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__CoordinateArray_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__CoordinateArray_O");
	Call_exposePython_geom__CoordinateArray_O(_lisp); // base(s): set(['core::Vector_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__OMatrix_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__OMatrix_O");
	Call_exposePython_geom__OMatrix_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__OVector2_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__OVector2_O");
	Call_exposePython_geom__OVector2_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__OVector3_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__OVector3_O");
	Call_exposePython_geom__OVector3_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_GeomPkg
#ifdef Use_GeomPkg
extern void Call_exposePython_geom__CoordinateArrayWithHash_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: geom__CoordinateArrayWithHash_O");
	Call_exposePython_geom__CoordinateArrayWithHash_O(_lisp); // base(s): set(['geom::CoordinateArray_O'])
}
#endif // ifdef Use_GeomPkg
#endif // EXPOSE_TO_PYTHON
#undef EXPOSE_TO_PYTHON
#if defined(EXPAND_CLASS_MACROS)
_CLASS_MACRO(geom::BoundingBox_O)
_CLASS_MACRO(geom::Color_O)
_CLASS_MACRO(geom::CoordinateArray_O)
_CLASS_MACRO(geom::OMatrix_O)
_CLASS_MACRO(geom::OVector2_O)
_CLASS_MACRO(geom::OVector3_O)
_CLASS_MACRO(geom::CoordinateArrayWithHash_O)
#endif // EXPAND_CLASS_MACROS
#undef ALL_STAGES
// ---------------- after class initializers
#ifdef EXPOSE_TO_PYTHON
#endif // EXPOSE_TO_PYTHON
