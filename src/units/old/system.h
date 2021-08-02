/*
    File: system.h
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
#ifndef	_units_System_H
#define _units_System_H

#include "core/common.h"
#include "unitsPackage.h"
#include "dimension.fwd.h"
#include "singleDimension.fwd.h"
#include "quantity.fwd.h"


namespace units
{

FORWARD(System);
class System_O : public Object_O
{
    LISP_CLASS(UnitsPkg,System_O,System,"Object_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(System_O);
public:
	void initialize();

public:
	static System_sp create(string const& systemName,
				uint distanceSymbolIndex,
				uint weightSymbolIndex,
				uint timeSymbolIndex,
				uint currentSymbolIndex,
				uint temperatureSymbolIndex,
				uint luminositySymbolIndex,
				uint amountSymbolIndex,
				core::LispPtr lisp);

private: // instance variables here
	/*! Name of this system */
	string			_Name;
	/*! What is each base dimension? */
	BaseUnit_sp	_BaseUnits[NumBaseDimensions];

public: // Functions here
	/*! Print the quantity in the standard units of this System*/
	string quantityAsString(double amount, Quantity_sp quantity) const;

	/*! Print the quantity in the standard units of this System*/
	string unitsAsString(double& amount, Quantity_sp quantity) const;

};

}; /* units */

TRANSLATE(units::System_O);

#endif /* _units_System_H */
