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

#include <clasp/core/common.h>
#include <cando/geom/vector3.h>
#include <cando/units/unit.h>
#include <cando/units/unitsPackage.h>

namespace units
{

FORWARD(Quantity);
class Quantity_O : public core::CxxObject_O
{
    LISP_CLASS(units,UnitsPkg,Quantity_O,"Quantity",core::CxxObject_O);
//    DECLARE_ARCHIVE();
//    DECLARE_SERIALIZE();
    DEFAULT_CTOR_DTOR(Quantity_O);
public:
	void initialize();

private: // instance variables here
	core::T_sp	_Value;
	Unit_sp		_Unit;

public:

    static Quantity_sp make(core::T_sp val, Unit_sp unit);

	/*! Create the quantity from a value and unit */
	static Quantity_sp create( core::T_sp val, Unit_sp unit);

	/*! Create the quantity from a value and unit */
	static Quantity_sp create( double val, Unit_sp unit);

public: // Handle different kinds of values
	static void throwOnInvalidValue(core::T_sp obj);
	/*! Create a deep copy of the value and scale it by the conversion factor */
	static core::T_sp copyAndScaleValue(core::T_sp val, double conversion);

	/*! Create a deep copy of the element and scale it by the conversion factor */
	static core::T_sp copyAndScaleValueElement(core::T_sp val, uint index, double conversion);

	/*! Create a deep copy of the element */
	static core::T_sp copyValueElement(core::T_sp val, uint index);

	/*! Set the element of a value */
	static void setValueElement(core::T_sp obj, uint index,	
				    core::T_sp newVal);

	/*! Return the size of the value */
	static uint sizeOfValue(core::T_sp val);

	/*! Test if value is nan */
	static bool isnanValue(core::T_sp val);

public: // Functions here

        /*! For serialization */
        bool fieldsp() const override { return true;}
        void fields(core::Record_sp node) override;
  
	/*! Return true if the Unit or Quantity other is compatible with my unit */
	bool is_compatible(core::T_sp other) const;

	/*Return the value in the desired Unit. Throw an exception if there
	  is a Unit mismatch */
	core::T_sp value_in_unit(Unit_sp other, int power=1) const;

	/*Return the value as a double. Throw an exception if there
	  is a Unit mismatch */
	double value_in_unit_asReal(Unit_sp other, int power=1) const;


	/*! Return the size of the value if it is an array */
	int size() const;

	/*! Get an element if the Value is a SimpleVectorCoordinate otherwise throw an exception */
	core::T_sp getElement_in_unit(uint index, Unit_sp other) const;

	/*! Get an element if the Value is a SimpleVectorCoordinate otherwise throw an exception */
	Quantity_sp getElement(uint index) const;

	/*! Get an element if the Value is a SimpleVectorCoordinate otherwise throw an exception */
	Vector3 getElement_in_unit_asVector3(uint index, Unit_sp other) const;

	/*! Set an element of the value with a Quantity */
	core::T_sp setElement(uint index, Quantity_sp other);

	string __repr__() const;

	/*! Return the raw value */
CL_LISPIFY_NAME("rawValue");
CL_DEFMETHOD 	core::T_sp rawValue() const { return this->_Value;};
	
	/*! Return the raw unit */
CL_LISPIFY_NAME("rawUnit");
CL_DEFMETHOD 	Unit_sp rawUnit() const { return this->_Unit;};

	/*! Return a string representation of the value:unit */
	string rawAsString() const;

	/*! Multiply Quantity by something else */
	core::T_sp operator*(core::T_sp other) const;

	/*! Divide Quantity by other value */
	core::T_sp operator/(core::T_sp other) const;

	/*! Add Quantity */
	core::T_sp operator+(core::T_sp other) const;

	/*! Subtract Quantity */
	core::T_sp operator-(core::T_sp other) const;

	/*! Return my Quantity to power (pow) */
	core::T_sp power(int pow) const;


	/*! Return the square root Quantity - throw an exception if there are powers that are not evenly divisible by 2 */
	core::T_sp sqrt() const;

	/*! Deep copy the Quantity */
	core::T_sp deepCopy() const;

	/*! Return true if the Value is nan */
	bool isnan() const; // This should be const


	/*! Comparison operators */
	bool operator<(core::T_sp obj) const;
	bool operator<=(core::T_sp obj) const;
	bool operator>(core::T_sp obj) const;
	bool operator>=(core::T_sp obj) const;
	
};

}; /* units */


#endif /* _units_Quantity_H */
