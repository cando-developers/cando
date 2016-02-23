/*
    File: openmmPackage.cc
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
#include <boost/mpl/list.hpp>


#include <foundation.h>
#include <object.h>
#include <lisp.h>
#include <clasp/core/corePackage.h>
#include <core/vector3.h>
#include <clasp/core/builtInClass.h>
#include <cando/openmmPackage/openmmPackage.h>
#include <cando/openmmPackage/serializer.h>
#include <cando/openmmPackage/openmmExpose.h>
#include <clasp/core/wrappers.h>


#define HEADER_INCLUDES
#include <initClasses.inc>
#undef HEADER_INCLUDES


using namespace core;




namespace omm
{

#define SYMBOLS_STATIC
#include <symbols_scraped.inc>


    extern void initializeOpenmmConstants(Lisp_sp lisp);

    void OpenMMExposer::expose(core::Lisp_sp lisp,core::PackageExposer::WhatToExpose what) const
    {
	switch (what)
	{
	case candoClasses:
	{
#define SYMBOLS_CREATE
#define SYMBOLS_EXPORT_LISP
#include <symbols_scraped.inc>

#define ALL_STAGES
#define Use_OmmPkg
#include <initClasses.inc>
#undef Use_OmmPkg
	    initializeSerializers(_lisp);
	}	
	break;
	case candoFunctions:
	{
	    //nothing
	};
	break;
	case candoGlobals:
	{
	initializeOpenmmConstants(_lisp);
	};
	break;
	case pythonClasses:
	case pythonFunctions:
	case pythonGlobals:
	{
	    IMPLEMENT_ME();
	}
	break;
	}
    }
	
};
