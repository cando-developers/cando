#ifndef	candoBasePackage_H
#define candoBasePackage_H

#include "core/foundation.h"
#include "core/lisp.h"
#include "candoBasePackage.fwd.h"

namespace candoBase
{


#define SYMBOLS_EXTERN
#include "symbols_scraped_inc.h"

    extern const char* CandoBase_nicknames[];

    class CandoBaseExposer : public core::Exposer
    {
    public:
    CandoBaseExposer(core::Lisp_sp lisp) : Exposer(lisp,AddonPkg,Addon_nicknames) {};
	virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
