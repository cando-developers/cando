#define DEBUG_LEVEL_FULL
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/package.h>
#include <clasp/core/environment.h>
#include <cando/candoBase/candoBasePackage.h>
#include <clasp/core/builtInClass.h>
#define HEADER_INCLUDES
#include <candoBase_initClasses_inc.h>


namespace kw {

#pragma GCC visibility push(default)
#define KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS
#pragma GCC visibility pop
};

namespace candoBase
{

#pragma GCC visibility push(default)
#define CandoBasePkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef CandoBasePkg_SYMBOLS
#pragma GCC visibility pop
};




namespace candoBase
{
#define EXPOSE_TO_CANDO
#define Use_CandoBasePkg
#define EXTERN_REGISTER
#include <candoBase_initClasses_inc.h>
#undef EXTERN_REGISTER
#undef Use_CandoBasePkg
#undef EXPOSE_TO_CANDO
};





namespace candoBase
{

    const char* CandoBase_nicknames[] = { "" };


#define SYMBOLS_STATIC
#include <symbols_scraped_inc.h>



    void CandoBaseExposer::expose(core::Lisp_sp lisp,WhatToExpose what) const
    {_G();
	switch (what)
	{
	case candoClasses:
	{
#define KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {kw::cname = _lisp->internWithPackageName(pkg,lispname); kw::cname->exportYourself(exportp);}
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS

#define CandoBasePkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {cname = _lisp->internWithPackageName(pkg,lispname); cname->exportYourself(exportp);}
#include <symbols_scraped_inc.h>
#undef DO_SYMBOL
#undef CandoBasePkg_SYMBOLS


#define ALL_STAGES
#define Use_CandoBasePkg
#define INVOKE_REGISTER
#define LOOKUP_SYMBOL(pkg,name) _lisp->internWithPackageName(pkg,name)
#include <candoBase_initClasses_inc.h>
#undef INVOKE_REGISTER
#undef Use_CandoBasePkg
#undef ALL_STAGES
	    /*! Use the core package */
	    this->package()->usePackage(lisp->findPackage(CorePkg));

	}
	break;
	case candoFunctions:
	{
	    // nothing
	}
	break;
	case candoGlobals:
	{
	}
	break;
	case pythonClasses:
	{
#define _DBG(x)
#define EXPOSE_TO_PYTHON
#define Use_CandoBasePkg
#include <candoBase_initClasses_inc.h>
#undef Use_CandoBasePkg
#undef EXPOSE_TO_PYTHON
	}
	break;
	case pythonGlobals:
	{
	    // nothing currently
	}
	break;
	case pythonFunctions:
	{
	    // nothing
	}
	break;
	};
    };
    
};





extern "C"
{
    bool ___user_libcandoBase()
    {
	static candoBase::candoBaseExposer* candoBasePkgP = NULL;
	if ( candoBasePkgP == NULL )
	{
	    printf("%s:%d - initializing candoBase library xxx\n", __FILE__, __LINE__ );
	    candoBasePkgP = new candoBase::candoBaseExposer(_lisp);
	    _lisp->installPackage(candoBasePkgP);
	    return true;
	} else return false;
    }

}




#if USE_INTRUSIVE_SMART_PTR==1
#define EXPAND_CLASS_MACROS
#define _CLASS_MACRO(_T_)				\
    STATIC_CLASS_INFO(_T_);			\
    INTRUSIVE_POINTER_REFERENCE_COUNT_ACCESSORS(_T_);
#include <candoBase_initClasses_inc.h>
#undef _CLASS_MACRO
#undef EXPAND_CLASS_MACROS
#endif
 
