#ifndef	unitsPackage_H
#define unitsPackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>




NAMESPACE_PACKAGE_ASSOCIATION(units,UnitsPkg,"UNITS");



namespace units
{


#define SYMBOLS_EXTERN
#include <cando/units/generated/symbols_scraped_inc.h>

    extern const char* Units_nicknames[];

    class UnitsExposer : public core::Exposer
    {
    public:
    UnitsExposer(core::Lisp_sp lisp) : Exposer(lisp,UnitsPkg,Units_nicknames) {};
	virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
