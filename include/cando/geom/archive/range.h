/*
    File: range.h
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
#ifndef	Range_H //[
#define Range_H



#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/numerics.h>
#include <cando/geom/geomPackage.fwd.h>

namespace geom {


    SMART(Range);
    class Range_O : public core::CxxObject_O
    {
	LISP_CLASS(geom,GeomPkg,Range_O,"Range",core::CxxObject_O);
	DECLARE_INIT();
    public:
	void	initialize();
    private:
	core::LongLongInt	_Begin;
	core::LongLongInt	_End;
	core::LongLongInt	_Step;
    public:
	static Range_sp create(core::LispPtr e,core::LongLongInt begin, core::LongLongInt end);
	static Range_sp create3(core::LispPtr e,core::LongLongInt begin, core::LongLongInt end, core::LongLongInt step);
    public:

	string __repr__() const;

	core::LongLongInt	begin() { return this->_Begin; };
	core::LongLongInt	end() { return this->_End; };
	core::LongLongInt	step() { return this->_Step; };
	core::Cons_sp asCons();

	Range_O( const Range_O& ss ); //!< Copy constructor

	DEFAULT_CTOR_DTOR(Range_O);
    };


};
TRANSLATE(geom::Range_O);
#endif //]
