#ifndef	chemPackage_H
#define chemPackage_H

#include "core/foundation.h"
#include "core/lisp.h"
#include "chemPackage.fwd.h"
#include "candoDatabase.fwd.h"

NAMESPACE_PACKAGE_ASSOCIATION(chemkw,ChemKwPkg,"CHEM-KEYWORD");
NAMESPACE_PACKAGE_ASSOCIATION(chem,ChemPkg,"CHEM");



namespace chem
{

    extern const char* Chem_nicknames[];

    FORWARD(Matter);
    typedef core::Symbol_sp MatterName;

    class ChemExposer : public core::Exposer
    {
    public:
    ChemExposer(core::Lisp_sp lisp) : Exposer(lisp,ChemPkg,Chem_nicknames) {};
        virtual void expose(core::Lisp_sp lisp, WhatToExpose what) const;
    };



    /*! Get the default CandoDatabase */
    CandoDatabase_sp getCandoDatabase();

};




#endif
