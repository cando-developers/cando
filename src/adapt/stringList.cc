/*
    File: stringList.cc
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
#include <cando/adapt/stringList.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/str.h>
#include <clasp/core/wrappers.h>

namespace adapt {





void StringList_O::setFromString(const string &ss) {
  
  VectorStrings parts;
  tokenize(ss, parts, " \n\t");
  this->setFromVectorStrings(parts);
}

CL_LISPIFY_NAME("asString");
CL_DEFMETHOD string StringList_O::asString() {
  
  stringstream ss;
  StringList_O::iterator si;
  string nm;
  ss.str("");
  for (si = this->_Contents.begin(); si != this->_Contents.end(); si++) {
    if (si != this->_Contents.begin()) {
      ss << " ";
    }
    ss << *si;
  }
  return ss.str();
}

#if defined(XML_ARCHIVE)
void StringList_O::archiveBase(ArchiveP node) {
  if (node->loading()) {
    VectorStrings v_Contents;
    node->getDataAsVectorOfStrings(v_Contents);
    this->_Contents.clear();
    this->setFromVectorStrings(v_Contents);
  } else {
    node->setCharacters(this->asString());
  }
}
#endif // defined(XML_ARCHIVE)

void StringList_O::clear() {
  LOG(BF("StringList::clear size=%d") % (this->_Contents.size()));
  if (this->_Contents.size() == 0)
    return;
  this->_Contents.clear();
}

void StringList_O::prepend(const string &str) {
  
  this->_Contents.insert(this->_Contents.begin(), str);
}

CL_LISPIFY_NAME("string-list-append");
CL_DEFMETHOD void StringList_O::append(const string &str) {
  
  this->_Contents.push_back(str);
}

void StringList_O::appendStringList(StringList_sp s) {
  
  StringList_O::iterator ic;
  for (ic = s->begin(); ic != s->end(); ic++) {
    this->_Contents.push_back(*ic);
  }
}

void StringList_O::appendConsOfStrings(List_sp s) {
  
  for (auto cur : s) {
    this->_Contents.push_back(gc::As<Str_sp>(oCar(cur))->get());
  }
}

void StringList_O::setFromVectorStrings(VectorStrings vs) {
  
  VectorStrings::iterator it;
  this->clear();
  for (it = vs.begin(); it != vs.end(); it++) {
    this->append(*it);
  }
}

bool StringList_O::contains(const string &nm) {
  
  StringList_O::iterator it;
  for (it = this->begin(); it != this->end(); it++) {
    if ((*it) == nm)
      return true;
  }
  return false;
}

int StringList_O::indexOf(const string &nm) {
  
  StringList_O::iterator it;
  int idx;
  LOG(BF("Looking for string(%s)") % nm.c_str());
  for (it = this->begin(), idx = 0; it != this->end(); it++, idx++) {
    LOG(BF("   Looking at StringList entry(%s)") % (*it).c_str());
    if ((*it) == nm) {
      LOG(BF("Found Match!!!"));
      return idx;
    }
  }
  SIMPLE_ERROR(BF("StringList does not contain: " + nm));
}
};
