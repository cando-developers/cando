#ifndef	geomPackage_H
#define geomPackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <cando/geom/geomPackage.fwd.h>

namespace geom
{


#define SYMBOLS_EXTERN
#include SYMBOLS_SCRAPED_INC_H

    extern const char* Geom_nicknames[];

    class GeomExposer : public core::Exposer
    {
    public:
    GeomExposer(core::Lisp_sp lisp) : Exposer(lisp,GeomPkg,Geom_nicknames) {};
	virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
