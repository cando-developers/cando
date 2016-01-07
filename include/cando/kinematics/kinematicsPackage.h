#ifndef	kinematicsPackage_H
#define kinematicsPackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>




NAMESPACE_PACKAGE_ASSOCIATION(kinematics,KinPkg,"KIN");



namespace kinematics
{


    class KinematicsExposer : public core::Exposer
    {
    private:
    public:
    KinematicsExposer(core::Lisp_sp lisp) : Exposer(lisp,KinPkg) {};
        virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
