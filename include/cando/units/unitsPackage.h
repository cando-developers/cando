/*
    File: unitsPackage.h
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
#ifndef	unitsPackage_H
#define unitsPackage_H

#include <clasp/core/foundation.h>
#include <clasp/core/lisp.h>


/*
A unit:unit stores  a scale factor (double) and a vector of powers for the base SI units/dimensions.
A unit:quantity stores a value and a unit:unit.   The value can be a double-float, a vector of doubles, or other unitless values.
Then there is unit:named-unit - it inherits from unit:unit and it has some extra slots to store the long and short form of the name of the unit.
Cando's unit package works on the principle that every unit can be written in terms of some base units.
Base units are called unit:dimension - they have an implicit scale factor of 1.0 and an implicit power of 1 for their dimension.
The unit:dimension's correspond to the SI units for each kind of physical quantity.
Conversion of unit:quantity to different unit:units is done with respect to the base units/dimensions.
  */

PACKAGE_NICKNAME("U");
PACKAGE_USE("COMMON-LISP");
PACKAGE_SHADOW("*");
PACKAGE_SHADOW("/");
PACKAGE_SHADOW("+");
PACKAGE_SHADOW("-");
PACKAGE_SHADOW("SQRT");
PACKAGE_SHADOW(">");
PACKAGE_SHADOW(">");
PACKAGE_SHADOW(">=");
PACKAGE_SHADOW("<");
PACKAGE_SHADOW("<=");
NAMESPACE_PACKAGE_ASSOCIATION(units,UnitsPkg,"UNITS");



namespace units
{


    class UnitsExposer_O : public core::Exposer_O
    {
      LISP_CLASS(units,UnitsPkg,UnitsExposer_O,"UnitsExposer",core::Exposer_O);
    public:
    UnitsExposer_O(core::Lisp_sp lisp) : core::Exposer_O(lisp,UnitsPkg) {};
	virtual void expose(core::Lisp_sp lisp,WhatToExpose what) const;
    };




};


#endif
