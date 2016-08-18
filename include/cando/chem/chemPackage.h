/*
    File: chemPackage.h
*/
/*
Open Source License
Copyright (c) 2016, Christian E. Schafmeister
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#ifndef	chemPackage_H
#define chemPackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>
#include <cando/chem/chemPackage.fwd.h>
#include <cando/chem/candoDatabase.fwd.h>

NAMESPACE_PACKAGE_ASSOCIATION(chemkw,ChemKwPkg,"KEYWORD");

// -------------------------------------
//
// Set up the CHEM package
//
PACKAGE_USE("COMMON-LISP");
NAMESPACE_PACKAGE_ASSOCIATION(chem,ChemPkg,"CHEM");



namespace chem
{

    extern const char* Chem_nicknames[];

    FORWARD(Matter);
    typedef core::Symbol_sp MatterName;

    class ChemExposer_O : public core::Exposer_O
    {
      LISP_CLASS(chem,ChemPkg,ChemExposer_O,"ChemExposer",core::Exposer_O);
    public:
    ChemExposer_O(core::Lisp_sp lisp) : Exposer_O(lisp,ChemPkg,Chem_nicknames) {};
        virtual void expose(core::Lisp_sp lisp, WhatToExpose what) const;
    };



    /*! Get the default CandoDatabase */
    CandoDatabase_sp getCandoDatabase();

};




#endif
