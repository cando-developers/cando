
#ifndef	openmmCando_H
#define openmmCando_H

#include "core/common.h"




NAMESPACE_PACKAGE_ASSOCIATION(omm,OmmPkg,"OMM");


namespace omm
{



#define SYMBOLS_EXTERN
#include "symbols_scraped.inc"



    class OpenMMExposer : public core::PackageExposer
    {
    private:
	int _stuff;
    public:
    OpenMMExposer(core::Lisp_sp lisp) : PackageExposer(lisp,OmmPkg) {};
        virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };



};
#endif
