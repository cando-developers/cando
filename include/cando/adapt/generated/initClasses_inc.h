// start
// define cpp macros: SET_SYMBOL, CREATE_CLASS, SET_CLASS, DEFINE_BASE_CLASSES, DEFINE_CLASS_NAMES, EXPOSE_TO_CANDO 
// define cpp macro: ALL_STAGES to get the effect of defining all of the macros above
// define cpp macro: EXPOSE_PYTHON to expose python
// Associating namespace(adapt) with package(AdaptPkg)
 // class IterateCons_O : public core::Iterator_O
 // class ObjectSet_O : public T_O {
 // class QDomNode_O : public core::T_O
 // class StringList_O : public T_O {
 // class StringSet_O : public core::T_O {
 // class SymbolList_O : public core::T_O {
 // class SymbolSet_O : public T_O {
#ifdef HEADER_INCLUDES
#include "include/iterateCons.h"
#include "include/objectSet.h"
#include "include/quickDom.h"
#include "include/stringList.h"
#include "include/stringSet.h"
#include "include/symbolList.h"
#include "include/symbolSet.h"
#endif // HEADER_INCLUDES
#undef HEADER_INCLUDES
#if defined(SET_SYMBOL) || defined(ALL_STAGES)
// requires LOOKUP_SYMBOL(pkg,symbolName) be defined
adapt::IterateCons_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::IterateCons_O::static_packageName(),adapt::IterateCons_O::static_className()));
adapt::ObjectSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::ObjectSet_O::static_packageName(),adapt::ObjectSet_O::static_className()));
adapt::QDomNode_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::QDomNode_O::static_packageName(),adapt::QDomNode_O::static_className()));
adapt::StringList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::StringList_O::static_packageName(),adapt::StringList_O::static_className()));
adapt::StringSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::StringSet_O::static_packageName(),adapt::StringSet_O::static_className()));
adapt::SymbolList_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::SymbolList_O::static_packageName(),adapt::SymbolList_O::static_className()));
adapt::SymbolSet_O::___set_static_ClassSymbol(LOOKUP_SYMBOL(adapt::SymbolSet_O::static_packageName(),adapt::SymbolSet_O::static_className()));
#endif // SET_SYMBOL
#undef SET_SYMBOL
#if defined(CREATE_CLASS) || defined(ALL_STAGES)

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
        core::LispObjectCreator<adapt::IterateCons_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::IterateCons_O>>::allocateClass();
        adapt::IterateCons_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::IterateCons_O::static_className() % (void*)(adapt::IterateCons_O::static_allocator) );
    classadapt__IterateCons_Oval->setCreator(adapt::IterateCons_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::IterateCons_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__IterateCons_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<adapt::ObjectSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::ObjectSet_O>>::allocateClass();
        adapt::ObjectSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::ObjectSet_O::static_className() % (void*)(adapt::ObjectSet_O::static_allocator) );
    classadapt__ObjectSet_Oval->setCreator(adapt::ObjectSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::ObjectSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__ObjectSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<adapt::QDomNode_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::QDomNode_O>>::allocateClass();
        adapt::QDomNode_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::QDomNode_O::static_className() % (void*)(adapt::QDomNode_O::static_allocator) );
    classadapt__QDomNode_Oval->setCreator(adapt::QDomNode_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::QDomNode_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__QDomNode_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<adapt::StringList_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::StringList_O>>::allocateClass();
        adapt::StringList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::StringList_O::static_className() % (void*)(adapt::StringList_O::static_allocator) );
    classadapt__StringList_Oval->setCreator(adapt::StringList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::StringList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__StringList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<adapt::StringSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::StringSet_O>>::allocateClass();
        adapt::StringSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::StringSet_O::static_className() % (void*)(adapt::StringSet_O::static_allocator) );
    classadapt__StringSet_Oval->setCreator(adapt::StringSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::StringSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__StringSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<adapt::SymbolList_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::SymbolList_O>>::allocateClass();
        adapt::SymbolList_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::SymbolList_O::static_className() % (void*)(adapt::SymbolList_O::static_allocator) );
    classadapt__SymbolList_Oval->setCreator(adapt::SymbolList_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::SymbolList_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__SymbolList_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */

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
        core::LispObjectCreator<adapt::SymbolSet_O>* cb = gctools::ClassAllocator<core::LispObjectCreator<adapt::SymbolSet_O>>::allocateClass();
        adapt::SymbolSet_O::___set_static_creator(cb);
    }
    LOG(BF("Set static_allocator for class(%s) to %X")% adapt::SymbolSet_O::static_className() % (void*)(adapt::SymbolSet_O::static_allocator) );
    classadapt__SymbolSet_Oval->setCreator(adapt::SymbolSet_O::static_creator);
    {
        LOG(BF("Created nil for class[%s]") % adapt::SymbolSet_O::static_className() );
    }
    /* ----- the class and its nil are now defined and so is classadapt__SymbolSet_Oval::___staticClass but the class _Slots and _Signature_ClassSlots are undefined - set them both to _Nil<T_O>() in stage3   ----- */
#endif // CREATE_CLASS
#undef CREATE_CLASS
#ifdef DUMP_INFO_CLASS // {
// Depends on nothing

    LOG(BF("---    dump_info   --- className: adapt::IterateCons_O @ %X") % classadapt__IterateCons_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::IterateCons_O::static_className() % adapt::IterateCons_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: adapt::ObjectSet_O @ %X") % classadapt__ObjectSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::ObjectSet_O::static_className() % adapt::ObjectSet_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: adapt::QDomNode_O @ %X") % classadapt__QDomNode_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::QDomNode_O::static_className() % adapt::QDomNode_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: adapt::StringList_O @ %X") % classadapt__StringList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::StringList_O::static_className() % adapt::StringList_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: adapt::StringSet_O @ %X") % classadapt__StringSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::StringSet_O::static_className() % adapt::StringSet_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: adapt::SymbolList_O @ %X") % classadapt__SymbolList_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::SymbolList_O::static_className() % adapt::SymbolList_O::static_classSymbol() );

    LOG(BF("---    dump_info   --- className: adapt::SymbolSet_O @ %X") % classadapt__SymbolSet_Oval.get());
    LOG(BF("%s::static_classSymbol() = %d") % adapt::SymbolSet_O::static_className() % adapt::SymbolSet_O::static_classSymbol() );
#endif // } DUMP_INFO_CLASS
#undef DUMP_INFO_CLASS
#if defined(DEFINE_BASE_CLASSES) || defined(ALL_STAGES) // {
// Depends on nothing
classadapt__IterateCons_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::Iterator_O::static_classSymbol());
classadapt__ObjectSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(adapt::T_O::static_classSymbol());
classadapt__QDomNode_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classadapt__StringList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(adapt::T_O::static_classSymbol());
classadapt__StringSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classadapt__SymbolList_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(core::T_O::static_classSymbol());
classadapt__SymbolSet_Oval->addInstanceBaseClassDoNotCalculateClassPrecedenceList(adapt::T_O::static_classSymbol());
#endif // } DEFINE_BASE_CLASSES
#undef DEFINE_BASE_CLASSES
#if defined(DEFINE_CLASS_NAMES) || defined(ALL_STAGES) // {
// Depends on nothing

    classadapt__IterateCons_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::IterateCons_O::static_classSymbol());

    classadapt__ObjectSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::ObjectSet_O::static_classSymbol());

    classadapt__QDomNode_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::QDomNode_O::static_classSymbol());

    classadapt__StringList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::StringList_O::static_classSymbol());

    classadapt__StringSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::StringSet_O::static_classSymbol());

    classadapt__SymbolList_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::SymbolList_O::static_classSymbol());

    classadapt__SymbolSet_Oval->__setupStage3NameAndCalculateClassPrecedenceList(adapt::SymbolSet_O::static_classSymbol());
#endif // } DEFINE_CLASS_NAMES
#undef DEFINE_CLASS_NAMES
#if defined(EXPOSE_TO_CANDO) || defined(ALL_STAGES)
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
extern void Register_adapt__QDomNode_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__QDomNode_O");
    adapt::Register_adapt__QDomNode_O(_lisp); // base(s): set(['core::T_O'])
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
extern void Register_adapt__StringSet_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__StringSet_O");
    adapt::Register_adapt__StringSet_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // INVOKE_REGISTER
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
#ifdef EXTERN_REGISTER
extern void Register_adapt__SymbolList_O(core::Lisp_sp); // base(s): set(['core::T_O'])
#endif // EXTERN_REGISTER
#ifdef INVOKE_REGISTER
{_BLOCK_TRACE("initializing Register_adapt__SymbolList_O");
    adapt::Register_adapt__SymbolList_O(_lisp); // base(s): set(['core::T_O'])
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
#endif // EXPOSE_TO_CANDO
#undef EXPOSE_TO_CANDO
#ifdef EXPOSE_TO_PYTHON
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__IterateCons_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__IterateCons_O");
	Call_exposePython_adapt__IterateCons_O(_lisp); // base(s): set(['core::Iterator_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__ObjectSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__ObjectSet_O");
	Call_exposePython_adapt__ObjectSet_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__QDomNode_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__QDomNode_O");
	Call_exposePython_adapt__QDomNode_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__StringList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__StringList_O");
	Call_exposePython_adapt__StringList_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__StringSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__StringSet_O");
	Call_exposePython_adapt__StringSet_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__SymbolList_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__SymbolList_O");
	Call_exposePython_adapt__SymbolList_O(_lisp); // base(s): set(['core::T_O'])
}
#endif // ifdef Use_AdaptPkg
#ifdef Use_AdaptPkg
extern void Call_exposePython_adapt__SymbolSet_O(::core::Lisp_sp lisp);
{_DBG("exposing to python: adapt__SymbolSet_O");
	Call_exposePython_adapt__SymbolSet_O(_lisp); // base(s): set(['adapt::T_O'])
}
#endif // ifdef Use_AdaptPkg
#endif // EXPOSE_TO_PYTHON
#undef EXPOSE_TO_PYTHON
#if defined(EXPAND_CLASS_MACROS)
_CLASS_MACRO(adapt::IterateCons_O)
_CLASS_MACRO(adapt::ObjectSet_O)
_CLASS_MACRO(adapt::QDomNode_O)
_CLASS_MACRO(adapt::StringList_O)
_CLASS_MACRO(adapt::StringSet_O)
_CLASS_MACRO(adapt::SymbolList_O)
_CLASS_MACRO(adapt::SymbolSet_O)
#endif // EXPAND_CLASS_MACROS
#undef ALL_STAGES
// ---------------- after class initializers
#ifdef EXPOSE_TO_PYTHON
#endif // EXPOSE_TO_PYTHON
