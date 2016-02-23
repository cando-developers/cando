/*
    File: gc_interface.h
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
#ifndef GC_INTERFACE_H
#define GC_INTERFACE_H
#include <clasp/core/foundation.h>

//
// All class forward declarations
//
namespace core {
class T_O;
class WrappedPointer_O;
class Functoid;
class Creator;
class Iterator_O;
    class SequenceStepper;
};
namespace asttooling {
    namespace internal {
        class MatcherDescriptor;
    };
};
namespace clbind {
class ConstructorCreator;
};


#ifdef USE_MPS
  #ifndef RUNNING_GC_BUILDER // when running the static analyzer - don't include the following
    #define DECLARE_FORWARDS
    #include GARBAGE_COLLECTION_INCLUDE
    #undef DECLARE_FORWARDS
  #endif // ifndef RUNNING_GC_BUILDER
#endif
#ifdef USE_BOEHM
  #ifdef USE_CXX_DYNAMIC_CAST
    // nothing
  #else
    #define DECLARE_FORWARDS
    #include GARBAGE_COLLECTION_INCLUDE
    #undef DECLARE_FORWARDS
  #endif
#endif

namespace gctools {

#ifdef USE_MPS
  #ifndef RUNNING_GC_BUILDER // when running the static analyzer - don't include the following
    #define GC_KIND_SELECTORS
    #include GARBAGE_COLLECTION_INCLUDE
    #undef GC_KIND_SELECTORS
  #endif // ifndef RUNNING_GC_BUILDER
#endif
#ifdef USE_BOEHM
  #ifdef USE_CXX_DYNAMIC_CAST
    // Nothing
  #else
    #define GC_KIND_SELECTORS
    #include GARBAGE_COLLECTION_INCLUDE
    #undef GC_KIND_SELECTORS
  #endif
#endif
};

#include <clasp/gctools/other_tagged_casts.h>


#endif
