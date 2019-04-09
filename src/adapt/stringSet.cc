/*
    File: stringSet.cc
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
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <cando/adapt/stringSet.h>
#include <clasp/core/symbolTable.h>
#include <clasp/core/multipleValues.h>
#include <cando/adapt/stringList.h>
#include <clasp/core/array.h>
#include <clasp/core/serialize.h>
#include <clasp/core/cons.h>
#include <clasp/core/numbers.h>
#include <clasp/core/wrappers.h>

namespace adapt {

/*
  __BEGIN_DOC(classes.StringSet.!class.StringSet)
  \requiredKeyed{entries:}{Cons::listOfStrings}

  Create a StringSet containing the strings in \sa{listOfStrings}.
  __END_DOC
*/

CL_LISPIFY_NAME("insertConsStrings");
CL_DEFMETHOD StringSet_sp StringSet_O::insertConsStrings(List_sp vals) {
  for (auto p : vals) {
    Str_sp t = gc::As<Str_sp>(oCar(p));
    this->insert(t->get());
  }
  return this->sharedThis<StringSet_O>();
}

void StringSet_O::insertVectorStrings(VectorStrings s) {
  VectorStrings::iterator i;
  for (i = s.begin(); i != s.end(); i++) {
    this->insert(*i);
  }
}

CL_LISPIFY_NAME("insertStringSet");
CL_DEFMETHOD void StringSet_O::insertStringSet(StringSet_sp s) {
  set<string>::iterator si;
  for (si = s->strs.begin(); si != s->strs.end(); si++) {
    this->insert(*si);
  }
}

void StringSet_O::insertStringList(StringList_sp s) {
  StringList_O::iterator si;
  for (si = s->begin(); si != s->end(); si++) {
    this->insert(*si);
  }
}

StringSet_sp StringSet_O::copy() {
  
  StringSet_O::iterator it;
  GC_COPY(StringSet_O, snew, *this);
  snew->clear();
  for (it = this->begin(); it != this->end(); it++)
    snew->insert(*it);
  return snew;
}

string StringSet_O::__repr__() const {
  stringstream ss;
  set<string>::iterator si;
  string nm;
  ss.str("");
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    if (si != this->strs.begin()) {
      ss << " ";
    }
    ss << *si;
  }
  return ss.str();
}
CL_LISPIFY_NAME("asString");
CL_DEFMETHOD string StringSet_O::asString() {
  stringstream ss;
  set<string>::iterator si;
  string nm;
  ss.str("");
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    if (si != this->strs.begin()) {
      ss << " ";
    }
    ss << *si;
  }
  return ss.str();
}

#if 0
void StringSet_O::archiveBase(ArchiveP node) {
  if (node->loading()) {
    Vector_sp vec = node->getVectorSNodes();
    this->strs.clear();
    for (int i(0), iEnd(vec->length()); i < iEnd; ++i) {
      LeafSNode_sp ln = gc::As<LeafSNode_sp>(vec->rowMajorAref(i));
      this->insert(gc::As<Str_sp>(ln->object())->get());
    }
  } else {
    SimpleVector_sp vec = SimpleVector_O::make(this->strs.size());
    int i(0);
    set<string>::iterator si;
    for (si = this->strs.begin(); si != this->strs.end(); si++) {
      (*vec)[i] = LeafSNode_O::create(Str_O::create(*si));
      ++i;
    }
    node->setVectorSNodesUnsafe(vec);
  }
}
#endif
CL_LISPIFY_NAME("contains");
CL_DEFMETHOD bool StringSet_O::contains(const string &s) {
  
  bool ye;
  ye = this->strs.count(s);
  return ye;
}

CL_LISPIFY_NAME("containsSubset");
CL_DEFMETHOD bool StringSet_O::containsSubset(StringSet_sp sub) {
  StringSet_O::iterator si;
  for (si = sub->begin(); si != sub->end(); si++) {
    if (!this->contains(*si)) {
      return false;
    }
  }
  return true;
}

CL_LISPIFY_NAME("clear");
CL_DEFMETHOD void StringSet_O::clear() {
  _OF();
  LOG(BF("StringSet::clear size=%d") % (this->strs.size()));
  if (this->strs.size() == 0)
    return;
  this->strs.clear();
}

/*! Check if the StringSet contains the same strings as this
 */
#if 0
bool StringSet_O::equal(T_sp obj) const {
  if (this->eq(obj))
    return true;
  if (gc::IsA<StringSet_sp>(obj)) {
    set<string>::iterator si;
    StringSet_sp ss = gc::As<StringSet_sp>(obj);
    if (this->strs.size() != ss->strs.size())
      return false;
    for (si = this->strs.begin(); si != this->strs.end(); si++) {
      if (!ss->contains(*si)) {
        return false;
      }
    }
    return true;
  }
  return false;
}
#endif

void StringSet_O::remove(const string &s) {
  _OF();
#ifdef DEBUG_ON
  if (this->strs.count(s) != 1) {
    SIMPLE_ERROR(BF("The string: %s was not found in StringSet") % s);
  }
#endif
  this->strs.erase(s);
}

void StringSet_O::initialize() {
  this->Base::initialize();
  this->rest = false;
}

StringSet_O::StringSet_O(const StringSet_O &ss) : CxxObject_O(ss) {
  this->rest = ss.rest;
  this->strs = ss.strs;
}

CL_LISPIFY_NAME("setFromString");
CL_DEFMETHOD void StringSet_O::setFromString(const string &s) {
  
  VectorStrings words;
  tokenize(s, words, "\t\n ");
  this->clear();
  this->insertVectorStrings(words);
}

CL_LISPIFY_NAME("StringSetUnion");
CL_DEFMETHOD StringSet_sp StringSet_O::setUnion(StringSet_sp b) {
  
  StringSet_sp nset;
  set<string>::iterator si;
  nset = StringSet_O::create();
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    nset->insert(*si);
  }
  for (si = b->strs.begin(); si != b->strs.end(); si++) {
    nset->insert(*si);
  }
  return nset;
}

CL_LISPIFY_NAME("StringSetIntersection");
CL_DEFMETHOD StringSet_sp StringSet_O::intersection(StringSet_sp b) {
  
  StringSet_sp nset;
  set<string>::iterator si;
  nset = StringSet_O::create();
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    LOG(BF("Looking for(%s)") % (*si).c_str());
    if (b->contains(*si)) {
      LOG(BF("Found it!!!"));
      nset->insert(*si);
    } else {
      LOG(BF("Not found"));
    }
  }
  return nset;
}

CL_LISPIFY_NAME("relativeComplement");
CL_DEFMETHOD StringSet_sp StringSet_O::relativeComplement(StringSet_sp b) {
  StringSet_sp nset;
  set<string>::iterator si;
  nset = StringSet_O::create();
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    if (!b->contains(*si)) {
      nset->insert(*si);
    }
  }
  return nset;
}

/*! Return a new set that takes every element of (this) in combination
  with every element in b separated by a comma
*/
CL_LISPIFY_NAME("cartesianProductInsert");
CL_DEFMETHOD StringSet_sp StringSet_O::cartesianProductInsert(string ins, StringSet_sp b) {
  
  StringSet_sp nset;
  set<string>::iterator si, bi;
  stringstream sstr;
  nset = StringSet_O::create();
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    for (bi = b->strs.begin(); bi != b->strs.end(); bi++) {
      sstr.str("");
      sstr << (*si) << ins << (*bi);
      nset->insert(sstr.str());
    }
  }
  return nset;
}

CL_LISPIFY_NAME("asCons");
CL_DEFMETHOD List_sp StringSet_O::asCons() const {
  
  List_sp cur = _Nil<T_O>();
  set<string>::iterator si;
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    Str_sp s = Str_O::create(*si);
    cur = Cons_O::create(s, cur);
  }
  return cur;
}

core::SimpleVector_sp StringSet_O::asVector() const {
  core::SimpleVector_sp vec = core::SimpleVector_O::make(this->strs.size());
  int i = 0;
  set<string>::iterator si;
  for (si = this->strs.begin(); si != this->strs.end(); si++) {
    core::SimpleBaseString_sp s = core::SimpleBaseString_O::make(*si);
    (*vec)[i] = s;
    ++i;
  }
  return vec;
}




#if 0
#ifdef USEBOOSTPYTHON
void	StringSet_O::python_setFromList(boost::python::list res) {
//string		i;
//boost::python::str		lval;
//VectorStrings	vs;
//VectorStrings::iterator vi;
  int			ll;
  string			x;

  ll = boost::python::extract<int>(res.attr("__len__")());
  for ( int i=0; i<ll; i++ ) {
    x = boost::python::extract<string>(res[i]);
//	printf( "Extracted element: %d = %s\n", i, x.c_str() );
    this->strs.insert(x);
  }
}



boost::python::list StringSet_O::python_asList() {
  string		i;
  int		id;
  boost::python::list	res;
  boost::python::str		lval;
  VectorStrings	vs;
  set<string>::iterator vi;

  for ( vi=this->strs.begin(); vi!=this->strs.end(); vi++ ) {
    i = *vi;
    LOG(BF("Appending %s to list") % (i.c_str() ) );
    lval = boost::python::str(i);
    res.append(lval);
  }
  return res;
}
#endif
#endif

};
