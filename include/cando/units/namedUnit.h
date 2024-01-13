/*
    File: namedUnit.h
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
/* @(#)namedUnit.h
 */

#ifndef _NAMEDUNIT_H
#define _NAMEDUNIT_H 1


#include <clasp/core/common.h>
#include <cando/units/unit.h>
#include <cando/units/namedUnit.h>

namespace units
{

FORWARD(NamedUnit);
class NamedUnit_O : public Unit_O
{
    LISP_CLASS(units,UnitsPkg,NamedUnit_O,"NamedUnit",Unit_O);
//    DECLARE_ARCHIVE();
//    DECLARE_SERIALIZE();
    DEFAULT_CTOR_DTOR(NamedUnit_O);
public:
	void initialize();

public:
    static NamedUnit_sp create(core::Symbol_sp symbolId, string const& label, string const& code);
    static NamedUnit_sp create(core::Symbol_sp symbolId, double scale, int baseDimension, int power, string const& label, string const& code);
    static NamedUnit_sp create(core::Symbol_sp symbolId, double scale, Unit_sp unit, string const& label, string const& code);

private: // instance variables here
	/*! Store the symbol for the Unit */
	core::Symbol_sp		_Symbol;
	/*! Store the long label for the Unit (ie: meters) */
	string			_Label;
	/*! Store the code for this Unit (ie: m) */
	string			_Code;

public: // Functions here

};

}; /* units */



#endif /* _NAMEDUNIT_H */
