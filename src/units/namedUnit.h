/* @(#)namedUnit.h
 */

#ifndef _NAMEDUNIT_H
#define _NAMEDUNIT_H 1


#include <clasp/core/common.h>
#include <cando/units/unit.h>
#include <cando/units/namedUnit.h>

namespace units
{

FORWARD(NamedUnit);
class NamedUnit_O : public Unit_O
{
    LISP_BASE1(Unit_O);
    LISP_CLASS(units,UnitsPkg,NamedUnit_O,"NamedUnit");
//    DECLARE_ARCHIVE();
//    DECLARE_SERIALIZE();
    DEFAULT_CTOR_DTOR(NamedUnit_O);
public:
	void initialize();

public:
    static NamedUnit_sp create(core::Symbol_sp symbolId, string const& label, string const& code, core::Lisp_sp lisp);
    static NamedUnit_sp create(core::Symbol_sp symbolId, double scale, int baseDimension, int power, string const& label, string const& code, core::Lisp_sp lisp);
    static NamedUnit_sp create(core::Symbol_sp symbolId, double scale, Unit_sp unit, string const& label, string const& code, core::Lisp_sp lisp);

private: // instance variables here
	/*! Store the symbol for the Unit */
	core::Symbol_sp		_Symbol;
	/*! Store the long label for the Unit (ie: meters) */
	string			_Label;
	/*! Store the code for this Unit (ie: m) */
	string			_Code;

public: // Functions here

};

}; /* units */

TRANSLATE(units::NamedUnit_O);


#endif /* _NAMEDUNIT_H */

