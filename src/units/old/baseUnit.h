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
	static BaseUnit_sp create(uint symbolIndex, double scale, string const& label, BaseDimensionEnum baseDim, core::Lisp_sp lisp);

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

TRANSLATE(units::BaseUnit_O);

#endif /* _units_BaseUnit_H */


