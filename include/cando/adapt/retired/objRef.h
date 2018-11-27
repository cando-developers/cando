/*
    File: objRef.h
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
#ifndef ObjRef_H //[
#define ObjRef_H

#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>

namespace core {

SMART(ObjRef);
 class ObjRef_O : public core::CxxObject_O {
  LISP_CLASS(core, CorePkg, ObjRef_O, "ObjRef",core::CxxObject_O);

public:
  void initialize();
GCPRIVATE:
  string _Selector;
  string _Name;
  T_sp _SubRef;

public:
  static ObjRef_sp create(Lisp_sp e, const string &asString);
  static ObjRef_sp create2(Lisp_sp e, const string &asString) { return ObjRef_O::create(e, asString); };

public:
  string getSelector() { return this->_Selector; };
  void setSelector(const string &t) { this->_Selector = t; };
  string getName() { return this->_Name; };
  void setName(const string &t) { this->_Name = t; };
  ObjRef_sp getSubRef() { return this->_SubRef; };
  void setSubRef(ObjRef_sp o) { this->_SubRef = o; };

  string asString();

  /*! Follow the reference relative to the given object
	 */
  T_sp relativeTo(T_sp o);

  DEFAULT_CTOR_DTOR(ObjRef_O);
};
};
TRANSLATE(core::ObjRef_O);
#endif //]
