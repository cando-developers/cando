/*
    File: float.h
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
#ifndef	_units_Float_H
#define _units_Float_H

#include "core/common.h"
#include "unitsPackage.h"
#include "system.fwd.h"
#include "singleDimension.fwd.h"
#include "quantity.h"

namespace units
{

FORWARD(Float);
class Float_O : public Quantity_O
{
    friend class System_O;
    friend class Xyz_O;
    LISP_CLASS(UnitsPkg,Float_O,Float,"Quantity_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Float_O);
public:
	void initialize();
public:

private: // instance variables here
	/*! Store the amount of the quantity */
	double	_Amount;

private:
public:
	static Float_sp create(double val, core::Symbol_sp dimension, int power, core::Symbol_sp system,core::Lisp_sp lisp);

public: // Functions here
	void setDimension(BaseUnit_sp dimension, int power);

	/*! Return the amount of this quantity - the quantities powers must all be zero
	  or an exception will be thrown. If they are all zero then we have a scalar
	  and return that amount */
	Object_sp amount() const;

	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	Object_sp amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp);

	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	double amount_matchUnit(core::Symbol_sp unit, int power);


	/*! Add the quantity to another quantity - this means powers must match */
	Quantity_sp add(Quantity_sp other) const;

	/*! Subtract the quantity from another quantity - powers must match */
	Quantity_sp sub(Quantity_sp other) const;

	/*! Negate the quantity */
	Quantity_sp negate() const;


	/*! Multiply the quantity by another quantity - this means adding powers */
	Quantity_sp mul(Quantity_sp other) const;

	/*! Divide the quantity by another quantity - this means adding powers */
	Quantity_sp div(Quantity_sp other) const;


	/*! Return the reciprocal
	 Throw exception if divide by zero */
	Quantity_sp reciprocal() const;

	/*! Take the square root,
	 Throw exception if powers arent multiples of 2*/
	Quantity_sp sqrt() const;


	/*! Multiply the quantity by a scalar */
	Quantity_sp mulByScalar(double scalar) const;

	/*! Divide the quantity by a scalar */
	Quantity_sp divByScalar(double scalar) const;


	/*! Raise the Float to an integer power - multiply the powers */
	Quantity_sp pow(int power) const;

	/*! Return true if the amount is NaN */
	bool isnan() const;

	/*! Return the amount and ignore the dimensions */
	Object_sp unsafeAmount() const;



	/*! Dump the contents with dimensions */
	virtual string dump() const;


	


	string __repr__() const;
};

}; /* units */

TRANSLATE(units::Float_O);

#endif /* _units_Float_H */
