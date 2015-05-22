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
	static Xyz_sp create(Vector3 pos, Quantity_sp temp, core::Lisp_sp lisp);

public: // Functions here
	void setDimension(BaseUnit_sp dimension, int power);

	/*! Return the amount of this quantity - the quantities powers must all be zero
	  or an exception will be thrown. If they are all zero then we have a scalar
	  and return that amount */
	Object_sp amount() const;

	/*! Return the amount of this quantity - if the base dimensions powers
	  that are provided as  arguments match my powers */
	Object_sp amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp);


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

TRANSLATE(units::Xyz_O);

#endif /* _units_Xyz_H */


