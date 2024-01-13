/*
    File: xyz.h
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
#ifndef	_units_Xyz_H
#define _units_Xyz_H

#include "core/common.h"
#include "core/vector3.h"
#include "unitsPackage.h"
#include "system.fwd.h"
#include "float.fwd.h"
#include "singleDimension.fwd.h"
#include "quantity.h"

namespace units
{

FORWARD(Xyz);
class Xyz_O : public Quantity_O
{
    friend class System_O;
    friend class CoordinateArray_O;
    LISP_CLASS(UnitsPkg,Xyz_O,Xyz,"Quantity_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Xyz_O);
public:
	void initialize();
public:

private: // instance variables here
	/*! Store the amount of the quantity */
	Vector3	_Pos;

private:
public:
	/*! Create an Xyz using the powers/system in template */
	static Xyz_sp create(Vector3 pos, Quantity_sp temp, core::LispPtr lisp);

public: // Functions here
	void setDimension(BaseUnit_sp dimension, int power);

	/*! Return the amount of this quantity - the quantities powers must all be zero
	  or an exception will be thrown. If they are all zero then we have a scalar
	  and return that amount */
	Object_sp amount() const;

	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	Object_sp amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp);


	/*! Add the quantity to another quantity - this means powers must match */
	Quantity_sp add(Quantity_sp other) const;

	/*! Subtract the quantity from another quantity - powers must match */
	Quantity_sp sub(Quantity_sp other) const;




	/*! Multiply the quantity by a scalar */
	Quantity_sp mulByFloat(Float_sp f) const;

	/*! Return true if the amount is NaN */
	bool isnan() const;

	/*! Return the amount and ignore the dimensions */
	Object_sp unsafeAmount() const;



	


	string __repr__() const;
};

}; /* units */


#endif /* _units_Xyz_H */
