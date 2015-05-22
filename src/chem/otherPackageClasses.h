//
//

LISP_OTHER_PACKAGE(CorePkg,addon::HighlightedObject_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::Iterator_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::ConditionBase_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::ConditionWithStackTrace_O,core::ConditionBase_O)
LISP_OTHER_PACKAGE(CorePkg,core::CandoException_O,core::ConditionWithStackTrace_O)
LISP_OTHER_PACKAGE(CorePkg,core::LispError_O,core::CandoException_O)
LISP_OTHER_PACKAGE(CorePkg,addon::GraphicsObject_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,addon::Render_O,addon::GraphicsObject_O)
LISP_OTHER_PACKAGE(CorePkg,core::SequenceBaseVirtual_O,core::T_O)
LISP_OTHER_PACKAGE(CorePkg,core::SequenceBase_O,core::SequenceBaseVirtual_O)
LISP_OTHER_PACKAGE(CorePkg,core::Sequence_O,core::SequenceBase_O)

