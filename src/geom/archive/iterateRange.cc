/*
    File: iterateRange.cc
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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE

#include <clasp/core/common.h>
#include <clasp/core/numbers.h>
#include <cando/geom/iterateRange.h>
#include <clasp/core/wrappers.h>

namespace geom {

IterateRange_sp IterateRange_O::create(core::LispPtr e, int ibegin, int iend) {
  auto r = gctools::GC<IterateRange_O>::allocate();
  r->_Begin = ibegin;
  r->_Step = 1;
  r->_End = iend;
  return r;
};

IterateRange_sp IterateRange_O::create(core::LispPtr e, int ibegin, int iend, int istep) {
  auto r = gctools::GC<IterateRange_O>::allocate();
  r->_Begin = ibegin;
  r->_Step = istep;
  r->_End = iend;
  return r;
};

void IterateRange_O::initialize() {
  this->Base::initialize();
  this->_Cur = core::clasp_make_fixnum(0);
}

void IterateRange_O::initRange(int b, int e) {
  this->_Begin = b;
  this->_End = e;
  this->_Step = 1;
}

void IterateRange_O::first() { this->_Cur->set(this->_Begin); };

void IterateRange_O::next() { this->_Cur = core::clasp_make_fixnum(this->_Cur->as_int() + this->_Step); };

bool IterateRange_O::isDone() {
  ASSERTNOTNULL(this->_Cur);
  return this->_Cur->get() >= this->_End;
};

core::T_sp IterateRange_O::currentObject() { return this->_Cur; };

}; // namespace geom
