#include <boost/mpl/list.hpp>


#include "foundation.h"
#include "object.h"
#include "lisp.h"
#include "core/corePackage.h"
#include "core/vector3.h"
#include "core/builtInClass.h"
#include "openmmPackage.h"
#include "serializer.h"
#include "openmmExpose.h"
#include "core/wrappers.h"


#define HEADER_INCLUDES
#include "initClasses.inc"
#undef HEADER_INCLUDES


using namespace core;




namespace omm
{

#define SYMBOLS_STATIC
#include "symbols_scraped.inc"


    extern void initializeOpenmmConstants(Lisp_sp lisp);

    void OpenMMExposer::expose(core::Lisp_sp lisp,core::PackageExposer::WhatToExpose what) const
    {_G();
	switch (what)
	{
	case candoClasses:
	{
#define SYMBOLS_CREATE
#define SYMBOLS_EXPORT_LISP
#include "symbols_scraped.inc"

#define ALL_STAGES
#define Use_OmmPkg
#include "initClasses.inc"
#undef Use_OmmPkg
	    initializeSerializers(_lisp);
	}	
	break;
	case candoFunctions:
	{
	    //nothing
	};
	break;
	case candoGlobals:
	{
	initializeOpenmmConstants(_lisp);
	};
	break;
	case pythonClasses:
	case pythonFunctions:
	case pythonGlobals:
	{
	    IMPLEMENT_ME();
	}
	break;
	}
    }
	
};






