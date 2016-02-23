/*
    File: symbolList.h
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

//
// (C) 2004 Christian E. Schafmeister
//

#ifndef SymbolList_H
#define SymbolList_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <cando/adapt/adaptPackage.fwd.h>
namespace adapt {
  using namespace core;
SMART(Residue);
SMART(SymbolSet);
SMART(SymbolList);

/*!
	A class that stores a set of strings
*/
SMART(SymbolList);
 class SymbolList_O : public core::CxxObject_O {
   LISP_CLASS(adapt, AdaptPkg, SymbolList_O, "SymbolList",core::CxxObject_O);

public:
  bool fieldsp() const { return true;};
  void fields(core::Record_sp node);
GCPRIVATE:
  gctools::Vec0<core::Symbol_sp> _Contents;

public:
  typedef gctools::Vec0<core::Symbol_sp>::iterator iterator;
  typedef gctools::Vec0<Symbol_sp>::const_iterator const_iterator;

  iterator begin() { return this->_Contents.begin(); };
  iterator end() { return this->_Contents.end(); };
  const_iterator begin() const { return this->_Contents.begin(); };
  const_iterator end() const { return this->_Contents.end(); };

  int size() { return this->_Contents.size(); };
  void prepend(Symbol_sp s);
  void append(Symbol_sp s);
  void appendSymbolList(SymbolList_sp strings);
  void appendConsOfStrings(List_sp strings);

  bool contains(Symbol_sp nm);
  int indexOf(Symbol_sp nm);
  Symbol_sp get(uint i) {
    _G();
    ASSERT_lt(i, this->_Contents.size());
    return this->_Contents[i];
  };
  void clear();

  void setFromVectorStrings(VectorStrings vs);
  string asString();

  SymbolSet_sp asStringSet();

  void map(std::function<void(Symbol_sp)> const &fn);

  SymbolList_O(const SymbolList_O &ss); //!< Copy constructor

  DEFAULT_CTOR_DTOR(SymbolList_O);
};
};

TRANSLATE(adapt::SymbolList_O);
#endif
