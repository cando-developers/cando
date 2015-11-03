#define DEBUG_LEVEL_FULL
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/builtInClass.h>
#include <clasp/core/lisp.h>
#include <cando/kinematics/atomTemplate.h>
#include <cando/kinematics/atomTree.h>
#include <cando/kinematics/kinematicsPackage.h>

#define HEADER_INCLUDES
#include INIT_CLASSES_INC_H
#undef HEADER_INCLUDES


namespace kinematics {
#pragma GCC visibility push(default)
#define KinPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkgName,lispName,export) core::Symbol_sp cname = UNDEFINED_SYMBOL;
#include SYMBOLS_SCRAPED_INC_H
#undef DO_SYMBOL
#undef KinPkg_SYMBOLS
#pragma GCC visibility pop
};


namespace kinematics
{
#define EXPOSE_TO_CANDO
#define Use_KinPkg
#define EXTERN_REGISTER
#include INIT_CLASSES_INC_H
#undef EXTERN_REGISTER
#undef Use_KinPkg
#undef EXPOSE_TO_CANDO
};






namespace kinematics
{


    SYMBOL_SC_(KinPkg,dihedrals);
    SYMBOL_SC_(KinPkg,sigmas);
    SYMBOL_SC_(KinPkg,indices);
    SYMBOL_SC_(KinPkg,count);
    SYMBOL_SC_(KinPkg,probability);
    SYMBOL_SC_(KinPkg,phi);
    SYMBOL_SC_(KinPkg,psi);
    SYMBOL_SC_(KinPkg,rotamer);
    SYMBOL_SC_(KinPkg,atom);
    SYMBOL_SC_(KinPkg,bonded);
    SYMBOL_SC_(KinPkg,delayed);
    SYMBOL_SC_(KinPkg,jump);
    SYMBOL_SC_(KinPkg,origin);
    SYMBOL_SC_(KinPkg,rootBonded);





    void KinematicsExposer::expose(core::Lisp_sp lisp,core::Exposer::WhatToExpose what) const
    {_G();
	switch (what)
	{
	case candoClasses:
	{_G();
#define KeywordPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {kw::cname = _lisp->internWithPackageName(pkg,lispname); kw::cname->exportYourself(exportp);}
#include SYMBOLS_SCRAPED_INC_H
#undef DO_SYMBOL
#undef KeywordPkg_SYMBOLS


#define KinPkg_SYMBOLS
#define DO_SYMBOL(cname,idx,pkg,lispname,exportp) {cname = _lisp->internWithPackageName(pkg,lispname); cname->exportYourself(exportp);}
#include SYMBOLS_SCRAPED_INC_H
#undef DO_SYMBOL
#undef KinPkg_SYMBOLS





#define ALL_STAGES
#define Use_KinPkg
#define INVOKE_REGISTER
#define LOOKUP_SYMBOL(pkg,name) _lisp->internWithPackageName(pkg,name)
#include INIT_CLASSES_INC_H
#undef INVOKE_REGISTER
#undef Use_KinPkg
#undef ALL_STAGES
	}
	break;
	case candoFunctions:
	{
	}
	break;
	case candoGlobals:
	{
	}
	break;
	case pythonClasses:
	{
	    IMPLEMENT_ME();
	}
	break;
	case pythonFunctions:
	{
	    IMPLEMENT_ME();
	}
	break;
	case pythonGlobals:
	{
	    IMPLEMENT_ME();
	}
	break;
	}
    }
	
};






extern "C"
{
    bool init_libkinematics()
    {
	static kinematics::KinematicsExposer* kinematicsPkgP = NULL;
	if ( kinematicsPkgP == NULL )
	{
	    kinematicsPkgP = new kinematics::KinematicsExposer(_lisp);
	    _lisp->installPackage(kinematicsPkgP);
	    return true;
	} else return false;
    }

}



#if USE_INTRUSIVE_SMART_PTR==1
#define EXPAND_CLASS_MACROS
#define _CLASS_MACRO(_T_)				\
    STATIC_CLASS_INFO(_T_);			\
    INTRUSIVE_POINTER_REFERENCE_COUNT_ACCESSORS(_T_);
#define Use_KinPkg
#include INIT_CLASSES_INC_H
#undef Use_KinPkg
#undef _CLASS_MACRO
#undef EXPAND_CLASS_MACROS
#endif
