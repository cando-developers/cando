/*
    File: baseUnit.h
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
#ifndef	_units_BaseUnit_H
#define _units_BaseUnit_H

#include "core/common.h"
#include "unitsPackage.h"
#include "dimension.h"
#include "dimension.fwd.h"
namespace units
{

FORWARD(BaseUnit);

/*! A BaseUnit is a unit of a single base dimension
 */
class BaseUnit_O : public Unit_O
{
    LISP_CLASS(UnitsPkg,BaseUnit_O,BaseUnit,"Unit_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(BaseUnit_O);
public:
	void initialize();

public:
	static BaseUnit_sp create(uint symbolIndex, double scale, string const& label, BaseDimensionEnum baseDim, core::LispPtr lisp);

private: // instance variables here
	/*! Store the index into the SIUnits  */
	BaseDimensionEnum	_BaseDimension;
	/*! Store the scale factor for this BaseUnit */
//	double			_UnitToCanonicalScale;
	/*! Store the inverse scale factor for this BaseUnit */
	double			_CanonicalToUnitScale;
	/*! Store the label for this BaseUnit */
	string			_Label;
	/*! Store the symbol for this BaseUnit */
	string			_Symbol
public: // Functions here
//	double unitToCanonicalScale() const { return this->_UnitToCanonicalScale;};
	double canonicalToUnitScale() const { return this->_CanonicalToUnitScale;};
	string label() const { return this->_Label;};


	void adjustPowersAndAmountScale(int power, int powers[NumBaseDimensions],double& amountScale) const;


};

}; /* units */


#endif /* _units_BaseUnit_H */
