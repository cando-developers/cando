/*
    File: symbolList.cc
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

#define DEBUG_LEVEL_NONE

//
// (C) 2004 Christian E. Schafmeister
//

#include <clasp/core/common.h>
#include <clasp/core/lisp.h>
#include <clasp/core/object.h>
#include <clasp/core/metaClass.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/array.h>
#include <cando/adapt/symbolList.h>
#include <cando/adapt/symbolSet.h>
#include <cando/adapt/adaptPackage.fwd.h>
#include <clasp/core/wrappers.h>


namespace adapt {





string SymbolList_O::asString() {
  
  stringstream ss;
  SymbolList_O::iterator si;
  string nm;
  ss.str("");
  for (si = this->_Contents.begin(); si != this->_Contents.end(); si++) {
    if (si != this->_Contents.begin()) {
      ss << " ";
    }
    ss << (*si)->formattedName(true);
  }
  return ss.str();
}

void SymbolList_O::fields(core::Record_sp node) {
  node->field_if_not_empty(INTERN_(kw,contents),this->_Contents);
}


void SymbolList_O::clear() {
  LOG(BF("SymbolList::clear size=%d") % (this->_Contents.size()));
  if (this->_Contents.size() == 0)
    return;
  this->_Contents.clear();
}

void SymbolList_O::prepend(Symbol_sp str) {
  
  this->_Contents.insert(this->_Contents.begin(), str);
}

void SymbolList_O::append(Symbol_sp str) {
  
  this->_Contents.push_back(str);
}

void SymbolList_O::appendSymbolList(SymbolList_sp s) {
  
  s->map([this](Symbol_sp s) {
                this->_Contents.push_back(s);
  });
}

void SymbolList_O::appendConsOfStrings(List_sp s) {
  

  for (auto cur : s) {
    Symbol_sp sym = _lisp->intern(gc::As<Str_sp>(oCar(cur))->get_std_string());
    this->_Contents.push_back(sym);
  }
}

void SymbolList_O::setFromVectorStrings(VectorStrings vs) {
  
  VectorStrings::iterator it;
  this->clear();
  for (it = vs.begin(); it != vs.end(); it++) {
    Symbol_sp sym = _lisp->intern(*it);
    this->append(sym);
  }
}

bool SymbolList_O::contains(Symbol_sp nm) {
  
  for (auto it = this->_Contents.begin(); it != this->_Contents.end(); it++) {
    if ((*it) == nm)
      return true;
  }
  return false;
}

int SymbolList_O::indexOf(Symbol_sp nm) {
  
  LOG(BF("Looking for string(%s)") % _rep_(nm));
  for (auto it = this->_Contents.begin(); it != this->_Contents.end(); it++) {
    LOG(BF("   Looking at SymbolList entry(%s)") % _rep_((*it)));
    if ((*it) == nm) {
      LOG(BF("Found Match!!!"));
      return it - this->_Contents.begin();
    }
  }
  SIMPLE_ERROR(BF("SymbolList does not contain: %s") % _rep_(nm));
}

void SymbolList_O::map(std::function<void(Symbol_sp)> const &fn) {
  for (auto it = this->_Contents.begin(); it != this->_Contents.end(); ++it) {
    Symbol_sp key = *it;
    fn(key);
  }
}
};
