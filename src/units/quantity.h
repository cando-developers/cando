#ifndef	_units_Quantity_H
#define _units_Quantity_H

#include "core/common.h"
#include "addon/vector3.h"
#include "unit.h"
#include "unitsPackage.h"

namespace units
{

FORWARD(Quantity);
class Quantity_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(units,UnitsPkg,Quantity_O,"Quantity");
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
	static Quantity_sp create( core::T_sp val, Unit_sp unit, core::Lisp_sp lisp );

	/*! Create the quantity from a value and unit */
	static Quantity_sp create( double val, Unit_sp unit, core::Lisp_sp lisp );

public: // Handle different kinds of values
	static void throwOnInvalidValue(core::T_sp obj,core::Lisp_sp lisp);
	/*! Create a deep copy of the value and scale it by the conversion factor */
	static core::T_sp copyAndScaleValue(core::T_sp val, double conversion, core::Lisp_sp lisp );

	/*! Create a deep copy of the element and scale it by the conversion factor */
	static core::T_sp copyAndScaleValueElement(core::T_sp val, uint index, double conversion, core::Lisp_sp lisp );

	/*! Create a deep copy of the element */
	static core::T_sp copyValueElement(core::T_sp val, uint index, core::Lisp_sp lisp );

	/*! Set the element of a value */
	static void setValueElement(core::T_sp obj, uint index,	
				    core::T_sp newVal, core::Lisp_sp lisp);

	/*! Return the size of the value */
	static uint sizeOfValue(core::T_sp val, core::Lisp_sp lisp);

	/*! Test if value is nan */
	static bool isnanValue(core::T_sp val, core::Lisp_sp lisp);

public: // Functions here
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

	/*! Get an element if the Value is a CoordinateArray otherwise throw an exception */
	core::T_sp getElement_in_unit(uint index, Unit_sp other) const;

	/*! Get an element if the Value is a CoordinateArray otherwise throw an exception */
	Quantity_sp getElement(uint index) const;

	/*! Get an element if the Value is a CoordinateArray otherwise throw an exception */
	Vector3 getElement_in_unit_asVector3(uint index, Unit_sp other) const;

	/*! Set an element of the value with a Quantity */
	core::T_sp setElement(uint index, Quantity_sp other);

	string __repr__() const;

	/*! Return the raw value */
	core::T_sp rawValue() const { return this->_Value;};
	
	/*! Return the raw unit */
	Unit_sp rawUnit() const { return this->_Unit;};

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

TRANSLATE(units::Quantity_O);

#endif /* _units_Quantity_H */


