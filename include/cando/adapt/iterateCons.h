/*
    File: iterateCons.h
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
#ifndef	IterateCons_H //[
#define IterateCons_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/iterator.h>
#include <cando/adapt/adaptPackage.fwd.h>

namespace adapt {


SMART(Cons);


SMART(IterateCons );
    class IterateCons_O : public core::Iterator_O
{
    LISP_CLASS(adapt,AdaptPkg,IterateCons_O,"IterateCons",core::Iterator_O);
public:
	void	initialize();
private:
	// instance variables
    core::List_sp	_Begin;
    core::List_sp	_Cur;
public:
    static IterateCons_sp create(core::List_sp beg)
    {
      auto  r  = gctools::GC<IterateCons_O>::allocate();
	r->init(beg);
	return r;
    };
public:

CL_NAME("INIT");
CL_DEFMETHOD     void init(core::List_sp c) { this->_Begin = c; };
	void first();
	void next();
	bool isDone();
    core::T_sp currentObject();



	DEFAULT_CTOR_DTOR(IterateCons_O);
};


};



#endif //]
