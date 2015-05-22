#ifndef	kinematicsPackage_H
#define kinematicsPackage_H

#include "core/foundation.h"
#include "core/lisp.h"




NAMESPACE_PACKAGE_ASSOCIATION(kinematics,KinPkg,"KIN");



namespace kinematics
{


#define SYMBOLS_EXTERN
#include "symbols_scraped_inc.h"


    class KinematicsExposer : public core::PackageExposer
    {
    private:
    public:
    KinematicsExposer(core::Lisp_sp lisp) : PackageExposer(lisp,KinPkg) {};
        virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
