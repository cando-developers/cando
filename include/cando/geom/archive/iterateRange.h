/*
    File: iterateRange.h
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
#ifndef	IterateRange_H //[
#define IterateRange_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <clasp/core/numbers.fwd.h>
#include <clasp/core/iterator.h>
#include <cando/geom/geomPackage.h>

namespace geom {

SMART(IterateRange );
class IterateRange_O : public core::Iterator_O
{
    LISP_CLASS(geom,GeomPkg,IterateRange_O,"IterateRange",core::Iterator_O);
public:
	void	initialize();
private:
	// instance variables
	int	_Begin;
	int	_Step;
	int	_End;
	core::Fixnum_sp	_Cur;
public:
	static core::T_sp prim_range(core::Function_sp e,core::Cons_sp args, core::Environment_sp environ, core::Lisp_sp lisp);


	static IterateRange_sp create(core::Lisp_sp e,int ibegin, int iend);
	static IterateRange_sp create(core::Lisp_sp e,int ibegin, int iend, int istep);
public:
	void initRange(int b,int e);
	void first();
	void next();
	bool isDone();
	core::T_sp currentObject();

	IterateRange_O( const IterateRange_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(IterateRange_O);
};


};
TRANSLATE(geom::IterateRange_O);
#endif //]
