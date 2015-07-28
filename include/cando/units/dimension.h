#ifndef	_units_Dimension_H
#define _units_Dimension_H

#include <clasp/core/common.h>
#include <cando/units/unitsPackage.h>
#include <cando/units/dimension.fwd.h>
namespace units
{

FORWARD(Dimension);
class Dimension_O : public core::CxxObject_O
{
    LISP_BASE1(core::CxxObject_O);
    LISP_CLASS(units,UnitsPkg,Dimension_O,"Dimension");
//    DECLARE_ARCHIVE();
    DEFAULT_CTOR_DTOR(Dimension_O);
public:
	void initialize();

public:

private: // instance variables here
public: // Functions here

	static void zeroPowers(int powers[]);
	static void copyPowers(int powers[], int const other[]);
	virtual void adjustPowersAndAmountScale(int power, int powers[NumBaseDimensions],double& amountScale) const;

	static string baseDimensionUnitName(int i);

};

}; /* units */

TRANSLATE(units::Dimension_O);

#endif /* _units_Dimension_H */


