/*
    File: stringSet.h
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

#ifndef StringSet_H
#define StringSet_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/holder.h>
#include <clasp/core/lispVector.h>
#include <clasp/core/cons.h>
#include <cando/adapt/adaptPackage.fwd.h>
namespace adapt {
  using namespace core;
SMART(Residue);
SMART(StringSet);
SMART(StringList);

/*!
	A class that stores a set of strings
*/
SMART(StringSet);
 class StringSet_O : public core::CxxObject_O {
   LISP_CLASS(adapt, AdaptPkg, StringSet_O, "StringSet",core::CxxObject_O);
  DECLARE_INIT();

public:
  void initialize();
  void archiveBase(core::ArchiveP node);

private:
  set<string> strs;
  bool rest;

public:
public:
  typedef set<string>::iterator iterator;

  iterator begin() { return this->strs.begin(); };
  iterator end() { return this->strs.end(); };

  StringSet_sp copy();

  string first() { return *(this->strs.begin()); };
CL_NAME("size");
CL_DEFMETHOD   uint size() { return this->strs.size(); };
  void remove(const string &s);
  bool contains(const string &s);
  bool containsSubset(StringSet_sp sub);
CL_NAME("insert");
CL_DEFMETHOD   void insert(const string &s) { this->strs.insert(s); };
  void insertVectorStrings(core::VectorStrings s);
  void insertStringSet(StringSet_sp ss);
  void insertStringList(StringList_sp ss);
  StringSet_sp insertConsStrings(core::List_sp list);
  void clear();

  core::List_sp asCons() const;
  Vector_sp asVector() const;

  bool equal(T_sp ss) const;

  // Set theory operations

  //! A setUnion B = (x:x E A || x E B)
  StringSet_sp setUnion(StringSet_sp b);

  //! A intersection B = (x:x E A && x E B)
  StringSet_sp intersection(StringSet_sp b);

  //! A-B = (x: x E A && not x E B )
  StringSet_sp relativeComplement(StringSet_sp b);

CL_NAME("removeAll");
CL_DEFMETHOD   StringSet_sp removeAll(StringSet_sp b) { return this->relativeComplement(b); };

  //! AxB = ("x,y": x E A ; y E B )
CL_NAME("cartesianProduct");
CL_DEFMETHOD   StringSet_sp cartesianProduct(StringSet_sp b) {
    return this->cartesianProductInsert("->", b);
  };
  StringSet_sp cartesianProductInsert(string ins, StringSet_sp b);

  void setFromString(const string &s);

#if 0
	template <class oclass>
	void	setFromMapKeys(const StringMap<oclass>& map)
	{
	    this->clear;
	    StringMap<oclass>::iterator it;
	    for ( it=map.begin(); it!=map.end(); map++ )
	    {
		this->insert(it->first);
	    }
	}
#endif
  string asString();
  std::ostream &dumpToStream(std::ostream &o);

  string __repr__() const;

  StringSet_O(const StringSet_O &ss); //!< Copy constructor

  DEFAULT_CTOR(StringSet_O);
  virtual ~StringSet_O() {};
};
};

TRANSLATE(adapt::StringSet_O);
#endif
