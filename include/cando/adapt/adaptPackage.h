#ifndef	adaptPackage_H
#define adaptPackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <cando/adapt/adaptPackage.fwd.h>

namespace adapt
{


#define SYMBOLS_EXTERN
  #ifndef SCRAPING
    #include SYMBOLS_SCRAPED_INC_H
  #endif

    extern const char* Adapt_nicknames[];

    class AdaptExposer : public core::Exposer
    {
    public:
    AdaptExposer(core::Lisp_sp lisp) : Exposer(lisp,AdaptPkg,Adapt_nicknames) {};
	virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
