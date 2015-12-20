#define DEBUG_LEVEL_FULL
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/package.h>
#include <clasp/core/environment.h>
#include <cando/geom/geomPackage.h>
#include <clasp/core/builtInClass.h>
#define HEADER_INCLUDES
#include <generated/initClasses_inc.h>


namespace geom
{
#pragma GCC visibility push(default)
#define GeomPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif
#undef DO_SYMBOL
#undef GeomPkg_SYMBOLS
#pragma GCC visibility pop
};




namespace geom
{
#define EXPOSE_TO_CANDO
#define Use_GeomPkg
#define EXTERN_REGISTER
#include <generated/initClasses_inc.h>
#undef EXTERN_REGISTER
#undef Use_GeomPkg
#undef EXPOSE_TO_CANDO
};





namespace geom
{

    const char* Geom_nicknames[] = { "" };


#define SYMBOLS_STATIC
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
  #endif



    void GeomExposer::expose(core::Lisp_sp lisp,WhatToExpose what) const
    {_G();
	switch (what)
	{
	case candoClasses:
	{

#define GeomPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {cname = _lisp->internWithPackageName(pkg,core::lispify_symbol_name(lispname)); cname->exportYourself(exportp);}
  #ifndef SCRAPING
#include <generated/symbols_scraped_inc.h>
            #endif
#undef DO_SYMBOL
#undef GeomPkg_SYMBOLS


#define ALL_STAGES
#define Use_GeomPkg
#define INVOKE_REGISTER
#define LOOKUP_SYMBOL(pkg,name) _lisp->internWithPackageName(pkg,name)
#include <generated/initClasses_inc.h>
#undef INVOKE_REGISTER
#undef Use_GeomPkg
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
#define Use_GeomPkg
#include <generated/initClasses_inc.h>
#undef Use_GeomPkg
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
    bool ___user_libgeom()
    {
	static geom::GeomExposer* geomPkgP = NULL;
	if ( geomPkgP == NULL )
	{
	    printf("%s:%d - initializing geom library xxx\n", __FILE__, __LINE__ );
	    geomPkgP = new geom::GeomExposer(_lisp);
	    _lisp->installPackage(geomPkgP);
	    return true;
	} else return false;
    }

}




#if USE_INTRUSIVE_SMART_PTR==1
#define EXPAND_CLASS_MACROS
#define _CLASS_MACRO(_T_)				\
    STATIC_CLASS_INFO(_T_);
#define Use_GeomPkg
#include <generated/initClasses_inc.h>
#undef Use_GeomPkg
#undef _CLASS_MACRO
#undef EXPAND_CLASS_MACROS
#endif
 
