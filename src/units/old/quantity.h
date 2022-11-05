/*
    File: quantity.h
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
#ifndef	_units_Quantity_H
#define _units_Quantity_H

#include "core/common.h"
#include "unitsPackage.h"
#include "system.fwd.h"
#include "dimension.fwd.h"
#include "singleDimension.fwd.h"

namespace units
{

FORWARD(Quantity);
class Quantity_O : public Object_O
{
    friend class System_O;
    friend class Float_O;
    friend class CoordinateArray_O;
    friend class Xyz_O;
    LISP_CLASS(UnitsPkg,Quantity_O,Quantity,"Object_O");
//    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Quantity_O);
public:
	void initialize();
public:

private: // instance variables here

	/*! Store the system to print this quantity in */
	core::Symbol_sp	_System;
	/*! Store the power to which each BaseUnit is taken */
	int	_Powers[NumBaseDimensions];

private:
	void throwOnPowersMisMatch(int matchPowers[]) const;
	/*! Apply parseDimensionSymbols repeatedly to units/powers to fill powers[] and amountScale */
	void parseDimensionSymbols(int powers[], double& amountScale,
				      core::Symbol_sp unit1, int power1 );
	void parseDimensionSymbols(int powers[], double& amountScale,
				      core::Symbol_sp unit1, int power1,
				      core::Symbol_sp unit2, int power2 );

	core::Symbol_sp parseDimensionSymbolList(int powers[], double& amountScale,
					       core::Cons_sp list, bool allowLastSystemArgument) const;

	bool powersMatch(const int otherPowers[]) const;

public: // Functions here
	/*Set one dimension and return the multiplier for the amount */
	double setDimension(BaseUnit_sp dimension, int power);

	/*! Return the amount of this quantity - the quantities powers must all be zero
	  or an exception will be thrown. If they are all zero then we have a scalar
	  and return that amount */
	virtual Object_sp amount() const;

	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	virtual Object_sp value_in_units_lisp(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp);


	/*! Return true if the powers all match - then the quantities have the same dimension */
	bool dimensionsMatch(Quantity_sp other) const;

	/*! Return the reciprocal
	 Throw exception if divide by zero */
	virtual Quantity_sp negate() const;


	/*! Return the reciprocal
	 Throw exception if divide by zero */
	virtual Quantity_sp reciprocal() const;

	/*! Take the square root,
	 Throw exception if powers arent multiples of 2*/
	virtual Quantity_sp sqrt() const;


	/*! Multiply the quantity by a scalar */
	virtual Quantity_sp mulByScalar(double scalar) const;

	/*! Divide the quantity by a scalar */
	virtual Quantity_sp divByScalar(double scalar) const;


	/*! Raise the Quantity to an integer power - multiply the powers */
	virtual Quantity_sp pow(int power) const;

	/*! Return true if the amount is NaN */
	virtual bool isnan() const;

	/*! Return the amount and ignore the dimensions */
	virtual Object_sp unsafeAmount() const;


	/*! Dump the contents with dimensions */
	virtual string dump() const {_OF(); SUBCLASS_MUST_IMPLEMENT();};


	/*! Return the system for this Quantity */
	System_sp system() const;
	


	string __repr__() const;
};


	/*! Add the quantity to another quantity - this means powers must match */
Quantity_sp two_arg_add(T_sp first, T_sp second);

	/*! Subtract the quantity from another quantity - powers must match */
Quantity_sp two_arg_sub(T_sp first, T_sp second);

	/*! Multiply the quantity by another quantity - this means adding powers */
Quantity_sp two_arg_mul(T_sp first, T_sp second);

	/*! Divide the quantity by another quantity - this means adding powers */
Quantity_sp two_arg_div(T_sp first, T_sp second);



}; /* units */

TRANSLATE(units::Quantity_O);

#endif /* _units_Quantity_H */
