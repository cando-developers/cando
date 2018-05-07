/*
    File: hold.h
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

#ifndef	Hold_H
#define	Hold_H

#include <cando/adapt/indexedObjectBag.h>

namespace core {
  typedef struct {
    int first_line;
    int first_column;
    int last_line;
    int last_column;
  } LispParserPos;


};


namespace chem {
  extern core::Symbol_sp& _sym_STARparserNodeHolderSTAR;

  template <class oClass>
    class Hold {
  public:
    core::Fixnum nodeHolderIndex;
    Hold( ) {
      gc::smart_ptr<oClass> obj = oClass::create();
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      core::Fixnum index = bag->set_next_available_entry(obj);
      this->nodeHolderIndex = index;
    };
    Hold( gctools::smart_ptr<oClass> obj ) {
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      core::Fixnum index = bag->set_next_available_entry(obj);
      this->nodeHolderIndex = index;
    };
    gc::smart_ptr<oClass> value() {
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      core::T_sp obj = bag->entry(this->nodeHolderIndex);
      return obj;
    }
    virtual ~Hold() {
      adapt::IndexedObjectBag_sp bag = gc::As<adapt::IndexedObjectBag_sp>(chem::_sym_STARparserNodeHolderSTAR->symbolValue());
      bag->remove_entry(this->nodeHolderIndex);
    };
  };

 /*!Use a GC compatible allocator here */
  template <class T>
    Hold<T>* createHolder(gc::smart_ptr<T> obj) {
    return new Hold<T>(obj);
  }

  template <class T>
    Hold<T>* createHolder() {
    return new Hold<T>();
  }

 /*! Use a corresponding GC compatible deallocator here */
  template <class T>
    void deleteHolder(Hold<T>* objP) {
    delete objP;
  };

};
#endif
