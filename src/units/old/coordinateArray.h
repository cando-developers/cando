/*
    File: coordinateArray.h
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
#ifndef	_units_CoordinateArray_H
#define _units_CoordinateArray_H

#include "core/common.h"
#include "core/coordinateArray.h"
#include "unitsPackage.h"
#include "system.fwd.h"
#include "singleDimension.fwd.h"
#include "quantity.h"
#include "xyz.fwd.h"

namespace units
{

FORWARD(CoordinateArray);
class CoordinateArray_O : public Quantity_O
{
    friend class System_O;
    LISP_CLASS(UnitsPkg,CoordinateArray_O,CoordinateArray,"Quantity_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(CoordinateArray_O);
public:
	void initialize();
public:

private: // instance variables here
	/*! Store the amount of the quantity */
	core::CoordinateArray_sp	_Coordinates;

private:
public: // Functions here

	static CoordinateArray_sp create(core::CoordinateArray_sp coords, core::Symbol_sp dimension, int power, core::Symbol_sp system, core::LispPtr lisp);

	void setDimension(BaseUnit_sp dimension, int power);

	/*! Return the amount of this quantity - the quantities powers must all be zero
	  or an exception will be thrown. If they are all zero then we have a scalar
	  and return that amount */
	Object_sp amount() const;


	int size() const { return this->_Coordinates->size(); };
	/*! Append a vector with units */
	void append(Xyz_sp vec);



	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	Object_sp amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp);

	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	Object_sp getElement_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::LispPtr lisp);

	/*! Return one element with the appropriate units */
	Xyz_sp getElement(int i) const;

	/*! Return one element with the appropriate distance units */
	Vector3 getElement(int i, core::Symbol_sp unit1, int power1);

	/*! Return one element with the appropriate units */
	Vector3 getElement(int i, core::Symbol_sp unit1, int power1
			   ,core::Symbol_sp unit2, int power2 );


	CoordinateArray_sp copy() const;



	/*! Return the amount and ignore the dimensions */
	Object_sp unsafeAmount() const;




	/*! Dump the contents with dimensions */
	virtual string dump() const;


	


	string __repr__() const;
};

}; /* units */

TRANSLATE(units::CoordinateArray_O);

#endif /* _units_CoordinateArray_H */
