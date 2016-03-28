/*
    File: stringList.h
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

#ifndef StringList_H
#define StringList_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/object.h>

namespace adapt {

SMART(StringList);

/*!
  A class that stores a set of strings
*/
SMART(StringList);
 class StringList_O : public core::General_O {
  LISP_CLASS(adapt, AdaptPkg, StringList_O, "StringList",core::General_O);

public:
#if defined(XML_ARCHIVE)
  void archiveBase(ArchiveP node);
#endif // defined(XML_ARCHIVE)
private:
  vector<string> _Contents;

public:
  typedef vector<string>::iterator iterator;

  iterator begin() { return this->_Contents.begin(); };
  iterator end() { return this->_Contents.end(); };

  string first() { return *(this->_Contents.begin()); };
  int size() { return this->_Contents.size(); };
  void prepend(const string &s);
  void append(const string &s);
  void appendStringList(StringList_sp strings);
  void appendConsOfStrings(List_sp strings);

  bool contains(const string &nm);
  int indexOf(const string &nm);
  string get(uint i) {
    _G();
    ASSERT_lt(i, this->_Contents.size());
    return this->_Contents[i];
  };
  void clear();

  void setFromVectorStrings(VectorStrings vs);
  void setFromString(const string &s);
  string asString();

  StringSet_sp asStringSet();

#ifdef USEBOOSTPYTHON
  void python_setFromList(boost::python::list res);
  boost::python::list python_asList();
#endif

  DEFAULT_CTOR_DTOR(StringList_O);
};
};

TRANSLATE(adapt::StringList_O);
#endif
