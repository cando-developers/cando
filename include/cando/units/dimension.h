/*
    File: dimension.h
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
#ifndef	_units_Dimension_H
#define _units_Dimension_H

#include <clasp/core/common.h>
#include <cando/units/unitsPackage.h>
#include <cando/units/dimension.fwd.h>
namespace units
{

FORWARD(Dimension);
class Dimension_O : public core::CxxObject_O
{
    LISP_CLASS(units,UnitsPkg,Dimension_O,"Dimension",core::CxxObject_O);
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Dimension_O);
public:
	void initialize();

public:

private: // instance variables here
public: // Functions here

	static void zeroPowers(int powers[]);
	static void copyPowers(int powers[], int const other[]);
	virtual void adjustPowersAndAmountScale(int power, int powers[NumBaseDimensions],double& amountScale) const;

	static string baseDimensionUnitName(int i);

};

}; /* units */

TRANSLATE(units::Dimension_O);

#endif /* _units_Dimension_H */
