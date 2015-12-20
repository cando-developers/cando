#define DEBUG_LEVEL_FULL
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/package.h>
#include <clasp/core/environment.h>
#include <cando/adapt/adaptPackage.h>
#include <clasp/core/builtInClass.h>
#define HEADER_INCLUDES
#include <generated/initClasses_inc.h>




namespace adapt
{

#pragma GCC visibility push(default)
#define AdaptPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif
#undef DO_SYMBOL
#undef AdaptPkg_SYMBOLS
#pragma GCC visibility pop
};




namespace adapt
{
#define EXPOSE_TO_CANDO
#define Use_AdaptPkg
#define EXTERN_REGISTER
#include <generated/initClasses_inc.h>
#undef EXTERN_REGISTER
#undef Use_AdaptPkg
#undef EXPOSE_TO_CANDO
};





namespace adapt
{

    const char* Adapt_nicknames[] = { "" };


#define SYMBOLS_STATIC
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif



    void AdaptExposer::expose(core::Lisp_sp lisp,WhatToExpose what) const
    {_G();
	switch (what)
	{
	case candoClasses:
	{
#if 0
#define KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {kw::cname = _lisp->internWithPackageName(pkg,core::lispify_symbol_name(lispname)); kw::cname->exportYourself(exportp);}
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS
#endif
          
#define AdaptPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {cname = _lisp->internWithPackageName(pkg,core::lispify_symbol_name(lispname)); cname->exportYourself(exportp);}
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif
#undef DO_SYMBOL
#undef AdaptPkg_SYMBOLS


#define ALL_STAGES
#define Use_AdaptPkg
#define INVOKE_REGISTER
#define LOOKUP_SYMBOL(pkg,name) _lisp->internWithPackageName(pkg,name)
#include <generated/initClasses_inc.h>
#undef INVOKE_REGISTER
#undef Use_AdaptPkg
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
#define Use_AdaptPkg
#include <generated/initClasses_inc.h>
#undef Use_AdaptPkg
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
    bool ___user_libadapt()
    {
	static adapt::AdaptExposer* adaptPkgP = NULL;
	if ( adaptPkgP == NULL )
	{
	    printf("%s:%d - initializing adapt library xxx\n", __FILE__, __LINE__ );
	    adaptPkgP = new adapt::AdaptExposer(_lisp);
	    _lisp->installPackage(adaptPkgP);
	    return true;
	} else return false;
    }

}




#if USE_INTRUSIVE_SMART_PTR==1
#define EXPAND_CLASS_MACROS
#define _CLASS_MACRO(_T_)				\
    STATIC_CLASS_INFO(_T_);
#define Use_AdaptPkg
#include <generated/initClasses_inc.h>
#undef Use_AdaptPkg
#undef _CLASS_MACRO
#undef EXPAND_CLASS_MACROS
#endif
 
