#ifndef	candoBasePackage_H
#define candoBasePackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <cando/candoBase/candoBasePackage.fwd.h>

namespace candoBase
{


#define SYMBOLS_EXTERN
#include <cando/candoBase/generated/symbols_scraped_inc.h>

    extern const char* CandoBase_nicknames[];

    class CandoBaseExposer : public core::Exposer
    {
    public:
    CandoBaseExposer(core::Lisp_sp lisp) : Exposer(lisp,CandoBasePkg,CandoBase_nicknames) {};
	virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
