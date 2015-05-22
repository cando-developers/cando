#ifndef	_units_System_H
#define _units_System_H

#include "core/common.h"
#include "unitsPackage.h"
#include "dimension.fwd.h"
#include "singleDimension.fwd.h"
#include "quantity.fwd.h"


namespace units
{

FORWARD(System);
class System_O : public Object_O
{
    LISP_CLASS(UnitsPkg,System_O,System,"Object_O");
    DECLARE_INIT();
    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(System_O);
public:
	void initialize();

public:
	static System_sp create(string const& systemName,
				uint distanceSymbolIndex,
				uint weightSymbolIndex,
				uint timeSymbolIndex,
				uint currentSymbolIndex,
				uint temperatureSymbolIndex,
				uint luminositySymbolIndex,
				uint amountSymbolIndex,
				core::Lisp_sp lisp);

private: // instance variables here
	/*! Name of this system */
	string			_Name;
	/*! What is each base dimension? */
	BaseUnit_sp	_BaseUnits[NumBaseDimensions];

public: // Functions here
	/*! Print the quantity in the standard units of this System*/
	string quantityAsString(double amount, Quantity_sp quantity) const;

	/*! Print the quantity in the standard units of this System*/
	string unitsAsString(double& amount, Quantity_sp quantity) const;

};

}; /* units */

TRANSLATE(units::System_O);

#endif /* _units_System_H */


