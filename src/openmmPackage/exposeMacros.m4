changecom(`//')

// ---------------------------------------------------------------------- 
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
//
// Expose openmm library
//
//
//
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// Define code 
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
// Divert to different streams for different aspects of the binding.

define(`fHeader',1)
define(`fExpose',2)
define(`fInit',3)
define(`fCode',4)
define(`fInitCode',5)
define(`fDeclarePredefinedSymbols',6)
define(`fStorePredefinedSymbols',7)
define(`fCreatePredefinedSymbols',8)
define(`fFunctions',9)
define(`fTranslate',10)

// ommSymbol(NoCutoff) --> declares _sym_ommPackage_NoCutoff

define(`ommSymbol',`
divert(fStorePredefinedSymbols)
SYMBOL_SC_($1);
divert
')


//
// omm_defun(_packageName_,_wrapperFunctionName_)
define(`omm_defun',`
divert(fFunctions)
core::af_def(OmmPkg,"$1",&$2);
divert
')


define(`ommGlobalInt',`
divert(fInit)
GLOBAL_INT($1,$2)
divert
')

define(`ommInitCode',`
divert(fInitCode)
    sInitCode << $1 << endl;
divert
')


define(`class_o',$1`_O')



define(`ommClass',`
define(`_M4_currentClass',$1)
define(`_M4_currentClassBase_O',$2)
divert // divert(fHeader)
SMART(_M4_currentClass);
class class_o(_M4_currentClass) : public _M4_currentClassBase_O
{
LISP_EXTERNAL_CLASS(OmmPkg,OpenMM::_M4_currentClass,class_o(_M4_currentClass),"_M4_currentClass",_M4_currentClassBase_O);
divert(fExpose)
EXPOSE_CLASS(omm,class_o(_M4_currentClass));
void class_o(_M4_currentClass)::exposePython(core::Lisp_sp lisp){}
void class_o(_M4_currentClass)::exposeCando(core::Lisp_sp lisp)
{
	core::externalClass_<class_o(_M4_currentClass)>(lisp->lisp())
divert(fTranslate)
divert
define(`_registerOmmClassMapping',`')
')









define(`ommClass_withTranslator',`
ommClass($1,$2)
divert(fExpose)
define(`_registerOmmClassMapping',`lisp->installGlobalInitializationCallback(&'class_o(_M4_currentClass)`::initializeGlobals);')
divert // divert(fHeader)
`public:
	static void initializeGlobals(core::Lisp_sp);
'
divert(fCode)
`void 'class_o(_M4_currentClass)`::initializeGlobals(core::Lisp_sp lisp)
{
	registerOmmClassMapping<'$1`,'class_o(_M4_currentClass)`>(lisp);
};'
divert
')


define(`lispInitialization',`
divert // divert(fHeader)
DECLARE_INIT()
divert(fCode)
`core::T_sp 'class_o(_M4_currentClass)`::__init__(core::Function_sp exec, core::Cons_sp args, core::Environment_sp bargs, core::Lisp_sp lisp)
{_G();
    try
    {
'
$2
`
    } catch (core::Condition& err )
    {
        err.conditionObject()->as<core::CandoException_O>()->appendMessage(BF("\n expected args: %s\n was passed: %s") % $1 % args->__repr__() );
	RETHROW(err.conditionObject());
    } catch (...)
    {
	throw;
    }
}
'
divert(fExpose)
`	.old_style_def_raw("core:__init__",&'class_o(_M4_currentClass)`::__init__,$1)'
divert
')


define(`ommMethod',`
divert(fExpose)
	.def("$1",&OpenMM::'_M4_currentClass`::$1)
divert
')

define(`oMethodDontExpose',`
divert // divert(fHeader)
public:
	$2;
divert
')

define(`oMethod',`
divert(fExpose)
	.def("$1",&class_o(_M4_currentClass)::$1)
divert // divert(fHeader)
	// This should be in the header
public:
	$2;
divert
')

define(`noWrapMethod',`
divert(fExpose)
	.old_style_def_raw("$1",&class_o(_M4_currentClass)::$1)
divert
')
define(`noWrapMethodWithArguments',`
divert(fExpose)
	.old_style_def_raw("$1",&class_o(_M4_currentClass)::$1,"$2")
divert
')



define(`end_ommClass',`
divert(fExpose)
;
_registerOmmClassMapping
divert // divert(fHeader)
DEFAULT_CTOR_DTOR(class_o(_M4_currentClass));
};
divert(fExpose)
}
define(`_M4_currentClass',`')
define(`_M4_currentClassBase_O',`')
divert
')



define(`end_ommClass_noDefaultCtorDtor',`
divert(fExpose)
;
_registerOmmClassMapping
divert // divert(fHeader)
};
divert(fExpose)
}
define(`_M4_currentClass',`')
define(`_M4_currentClassBase_O',`')
divert
')



define(`wrappedObject',`private: 
public:
       typedef $1 WrappedType;
       typedef $1* PointerToWrappedType;
       typedef PointerToWrappedType ExternalType;
       PointerToWrappedType $2;
public:
	virtual void* externalObject() const { return this->$2;};
	void set_externalObject(void* p) { this->$2 = static_cast<PointerToWrappedType>(p);};
	PointerToWrappedType wrappedPtr() { return this->$2;}
	')


define(`instanceVariable',`private:
'	$1;)

