/*
    File: unitsPackage.cc
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

This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University
at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define DEBUG_LEVEL_NONE
#include <clasp/core/foundation.h>
#include <clasp/core/object.h>
#include <clasp/core/lisp.h>
#include <clasp/core/corePackage.fwd.h>
#include <clasp/core/package.h>
#include <clasp/core/symbolTable.h>
#include <cando/units/unit.h>
#include <cando/units/quantity.h>
#include <cando/units/namedUnit.h>
#include <cando/units/unitsPackage.h>

namespace units {
#if 0
#define EXPOSE_TO_CANDO
#define Use_UnitsPkg
#define EXTERN_REGISTER
#include <clasp/core/initClasses.h>
#undef EXTERN_REGISTER
#undef Use_UnitsPkg
#undef EXPOSE_TO_CANDO
#endif
}; // namespace units

namespace units {

const char *Units_nicknames[] = {""};

SYMBOL_EXPORT_SC_(UnitsPkg, meters);
SYMBOL_EXPORT_SC_(UnitsPkg, kilometers);
SYMBOL_EXPORT_SC_(UnitsPkg, millimeters);
SYMBOL_EXPORT_SC_(UnitsPkg, centimeters);
SYMBOL_EXPORT_SC_(UnitsPkg, micrometers);
SYMBOL_EXPORT_SC_(UnitsPkg, nanometers);
SYMBOL_EXPORT_SC_(UnitsPkg, picometers);
SYMBOL_EXPORT_SC_(UnitsPkg, angstroms);

SYMBOL_EXPORT_SC_(UnitsPkg, liters);
SYMBOL_EXPORT_SC_(UnitsPkg, milliliters);
SYMBOL_EXPORT_SC_(UnitsPkg, microliters);
SYMBOL_EXPORT_SC_(UnitsPkg, nanoliters);
SYMBOL_EXPORT_SC_(UnitsPkg, picoliters);
SYMBOL_EXPORT_SC_(UnitsPkg, femtoliters);

SYMBOL_EXPORT_SC_(UnitsPkg, molar);
SYMBOL_EXPORT_SC_(UnitsPkg, millimolar);
SYMBOL_EXPORT_SC_(UnitsPkg, micromolar);
SYMBOL_EXPORT_SC_(UnitsPkg, nanomolar);
SYMBOL_EXPORT_SC_(UnitsPkg, picomolar);
SYMBOL_EXPORT_SC_(UnitsPkg, femtomolar);

SYMBOL_EXPORT_SC_(UnitsPkg, kilograms);
SYMBOL_EXPORT_SC_(UnitsPkg, grams);
SYMBOL_EXPORT_SC_(UnitsPkg, milligrams);
SYMBOL_EXPORT_SC_(UnitsPkg, micrograms);
SYMBOL_EXPORT_SC_(UnitsPkg, nanograms);
SYMBOL_EXPORT_SC_(UnitsPkg, daltons);

SYMBOL_EXPORT_SC_(UnitsPkg, seconds);
SYMBOL_EXPORT_SC_(UnitsPkg, milliseconds);
SYMBOL_EXPORT_SC_(UnitsPkg, microseconds);
SYMBOL_EXPORT_SC_(UnitsPkg, nanoseconds);
SYMBOL_EXPORT_SC_(UnitsPkg, picoseconds);
SYMBOL_EXPORT_SC_(UnitsPkg, femtoseconds);

SYMBOL_EXPORT_SC_(UnitsPkg, per_picosecond);

SYMBOL_EXPORT_SC_(UnitsPkg, degrees);
SYMBOL_EXPORT_SC_(UnitsPkg, radians);

SYMBOL_EXPORT_SC_(UnitsPkg, amperes);
SYMBOL_EXPORT_SC_(UnitsPkg, coulomb);
SYMBOL_EXPORT_SC_(UnitsPkg, elementaryCharge);

SYMBOL_EXPORT_SC_(UnitsPkg, kelvin);

SYMBOL_EXPORT_SC_(UnitsPkg, candelas);

SYMBOL_EXPORT_SC_(UnitsPkg, moles);
SYMBOL_EXPORT_SC_(UnitsPkg, millimoles);
SYMBOL_EXPORT_SC_(UnitsPkg, micromoles);
SYMBOL_EXPORT_SC_(UnitsPkg, nanomoles);
SYMBOL_EXPORT_SC_(UnitsPkg, picomoles);
SYMBOL_EXPORT_SC_(UnitsPkg, femtomoles);

// MultiDimensions
SYMBOL_EXPORT_SC_(UnitsPkg, Joule);
SYMBOL_EXPORT_SC_(UnitsPkg, joule);
SYMBOL_EXPORT_SC_(UnitsPkg, Joules);
SYMBOL_EXPORT_SC_(UnitsPkg, joules);
SYMBOL_EXPORT_SC_(UnitsPkg, kiloJoule);
SYMBOL_EXPORT_SC_(UnitsPkg, kilojoule);
SYMBOL_EXPORT_SC_(UnitsPkg, kiloJoules);
SYMBOL_EXPORT_SC_(UnitsPkg, kilojoules);
SYMBOL_EXPORT_SC_(UnitsPkg, kiloJoules_per_mole);
SYMBOL_EXPORT_SC_(UnitsPkg, kilojoules_per_mole);
SYMBOL_EXPORT_SC_(UnitsPkg, kiloCalories_per_mole);
SYMBOL_EXPORT_SC_(UnitsPkg, kilocalorie);
SYMBOL_EXPORT_SC_(UnitsPkg, kilocalories_per_mole);
SYMBOL_EXPORT_SC_(UnitsPkg, nm_per_ps);

void UnitsExposer_O::expose(core::LispPtr lisp, WhatToExpose what) const {
  switch (what) {
  case candoClasses: {
    /*! Use the core package */
    //          this->package()->usePackage(lisp->findPackage(CorePkg));
  } break;
  case candoFunctions: {
    // nothing
  } break;
  case candoGlobals: {
#define DEFDIM(sym, scale, code, dim, pwr)                                                                                         \
  Unit_sp unit_##sym =                                                                                                             \
      (_sym_##sym)->exportYourself()->defparameter(NamedUnit_O::create(_sym_##sym, scale, dim, pwr, code, code)).as<Unit_O>();

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    DEFDIM(meters, 1.0, "m", distanceMeter, 1);
    DEFDIM(kilometers, 1000.0, "km", distanceMeter, 1);
    DEFDIM(millimeters, 1.0e-3, "mm", distanceMeter, 1);
    DEFDIM(centimeters, 1.0e-2, "cm", distanceMeter, 1);
    DEFDIM(micrometers, 1.0e-6, "um", distanceMeter, 1);
    DEFDIM(nanometers, 1.0e-9, "nm", distanceMeter, 1);
    DEFDIM(picometers, 1.0e-12, "pm", distanceMeter, 1);
    DEFDIM(angstroms, 1.0e-10, "Angstroms", distanceMeter, 1);

    DEFDIM(liters, 0.001, "l", distanceMeter, 3);
    DEFDIM(milliliters, 0.001 * 1.0e-3, "ml", distanceMeter, 3);
    DEFDIM(microliters, 0.001 * 1.0e-6, "ul", distanceMeter, 3);
    DEFDIM(nanoliters, 0.001 * 1.0e-9, "nl", distanceMeter, 3);
    DEFDIM(picoliters, 0.001 * 1.0e-12, "pl", distanceMeter, 3);
    DEFDIM(femtoliters, 0.001 * 1.0e-15, "fl", distanceMeter, 3);

    DEFDIM(kilograms, 1.0, "kg", weightKilogram, 1);
    DEFDIM(grams, 1.0e-3, "g", weightKilogram, 1);
    DEFDIM(milligrams, 1.0e-6, "mg", weightKilogram, 1);
    DEFDIM(micrograms, 1.0e-9, "ug", weightKilogram, 1);
    DEFDIM(nanograms, 1.0e-12, "ng", weightKilogram, 1);

    DEFDIM(seconds, 1.0, "s", timeSecond, 1);
    DEFDIM(milliseconds, 1.0e-3, "ms", timeSecond, 1);
    DEFDIM(microseconds, 1.0e-6, "us", timeSecond, 1);
    DEFDIM(nanoseconds, 1.0e-9, "ns", timeSecond, 1);
    DEFDIM(picoseconds, 1.0e-12, "ps", timeSecond, 1);
    DEFDIM(femtoseconds, 1.0e-15, "fs", timeSecond, 1);

    DEFDIM(per_picosecond, 1.0e+12, "per_ps", timeSecond, -1);

    DEFDIM(radians, 1.0, "deg", radianAngle, 1);
    DEFDIM(degrees, 0.017453292519943295, "rad", radianAngle, 1); // 0.017... = pi/180

    DEFDIM(amperes, 1.0, "A", currentAmpere, 1);
    NamedUnit_sp unit_coulomb = NamedUnit_O::create(_sym_coulomb, "coulomb", "C");
    unit_coulomb->incorporateUnit(unit_amperes, 1.0, 1);
    unit_coulomb->incorporateUnit(unit_seconds, 1.0, -1);
    _sym_coulomb->exportYourself()->defparameter(unit_coulomb);

    NamedUnit_sp unit_elementaryCharge =
        NamedUnit_O::create(_sym_elementaryCharge, 1.60217656535e-19, unit_coulomb, "elementaryCharge", "e");
    _sym_elementaryCharge->exportYourself()->defparameter(unit_elementaryCharge);

    DEFDIM(kelvin, 1.0, "K", temperatureKelvin, 1);

    DEFDIM(candelas, 1.0, "cd", luminosityCandela, 1);

    DEFDIM(moles, 1.0, "mol", amountMole, 1);
    DEFDIM(millimoles, 1.0e-3, "mmol", amountMole, 1);
    DEFDIM(micromoles, 1.0e-6, "umol", amountMole, 1);
    DEFDIM(nanomoles, 1.0e-9, "nmol", amountMole, 1);
    DEFDIM(picomoles, 1.0e-12, "pmol", amountMole, 1);
    DEFDIM(femtomoles, 1.0e-15, "fmol", amountMole, 1);

    NamedUnit_sp unit_molar = NamedUnit_O::create(_sym_molar, "molar", "M");
    unit_molar->incorporateUnit(unit_moles, 1.0, 1);
    unit_molar->incorporateUnit(unit_liters, 1.0, -1);
    _sym_molar->exportYourself()->defparameter(unit_molar);
    NamedUnit_sp unit_millimolar = NamedUnit_O::create(_sym_molar, "millimolar", "mM");
    unit_millimolar->incorporateUnit(unit_moles, 1.0e-3, 1);
    unit_millimolar->incorporateUnit(unit_liters, 1.0, -1);
    _sym_millimolar->exportYourself()->defparameter(unit_millimolar);
    NamedUnit_sp unit_micromolar = NamedUnit_O::create(_sym_molar, "micromolar", "uM");
    unit_micromolar->incorporateUnit(unit_moles, 1.0e-6, 1);
    unit_micromolar->incorporateUnit(unit_liters, 1.0, -1);
    _sym_micromolar->exportYourself()->defparameter(unit_micromolar);
    NamedUnit_sp unit_nanomolar = NamedUnit_O::create(_sym_molar, "nanomolar", "nM");
    unit_nanomolar->incorporateUnit(unit_moles, 1.0e-9, 1);
    unit_nanomolar->incorporateUnit(unit_liters, 1.0, -1);
    _sym_nanomolar->exportYourself()->defparameter(unit_nanomolar);
    NamedUnit_sp unit_picomolar = NamedUnit_O::create(_sym_molar, "picomolar", "pM");
    unit_picomolar->incorporateUnit(unit_moles, 1.0e-12, 1);
    unit_picomolar->incorporateUnit(unit_liters, 1.0, -1);
    _sym_picomolar->exportYourself()->defparameter(unit_picomolar);
    NamedUnit_sp unit_femtomolar = NamedUnit_O::create(_sym_molar, "femtomolar", "fM");
    unit_femtomolar->incorporateUnit(unit_moles, 1.0e-15, 1);
    unit_femtomolar->incorporateUnit(unit_liters, 1.0, -1);
    _sym_femtomolar->exportYourself()->defparameter(unit_femtomolar);

    NamedUnit_sp unit_joules = NamedUnit_O::create(_sym_joules, "joules", "J");
    unit_joules->incorporateUnit(unit_kilograms, 1.0, 1);
    unit_joules->incorporateUnit(unit_meters, 1.0, 2);
    unit_joules->incorporateUnit(unit_seconds, 1.0, -2);
    _sym_joules->exportYourself()->defparameter(unit_joules);
    NamedUnit_sp unit_kilojoules = NamedUnit_O::create(_sym_kilojoules, 1000.0, unit_joules, "kilojoules", "kJ");
    _sym_kilojoules->exportYourself()->defparameter(unit_kilojoules);

    NamedUnit_sp unit_kilojoules_per_mole = NamedUnit_O::create(_sym_kilojoules_per_mole,
                                                                1.0, unit_kilojoules,
                                                                "kiloJoules_per_mole","kj/mole");
    unit_kilojoules_per_mole->incorporateUnit(unit_moles,1.0,-1);
    _sym_kilojoules_per_mole->exportYourself()->defparameter(unit_kilojoules_per_mole);


    NamedUnit_sp unit_kilocalorie = NamedUnit_O::create(_sym_kilocalorie, 4.184, unit_kilojoules, "kilocalorie", "kcal");
    _sym_kilocalorie->defparameter(unit_kilocalorie);

    NamedUnit_sp unit_kilocalories_per_mole =
        NamedUnit_O::create(_sym_kilocalories_per_mole, 4.184, unit_kilojoules_per_mole, "kilocalories_per_mole", "kcal/mole");
    _sym_kilocalories_per_mole->exportYourself()->defparameter(unit_kilocalories_per_mole);

    NamedUnit_sp unit_nm_per_ps = NamedUnit_O::create(_sym_nm_per_ps, "nanometers_per_picosecond", "nm/ps");
    unit_nm_per_ps->incorporateUnit(unit_nanometers, 1.0, 1);
    unit_nm_per_ps->incorporateUnit(unit_picoseconds, 1.0, -1);
    _sym_nm_per_ps->exportYourself()->defparameter(unit_nm_per_ps);

    NamedUnit_sp unit_daltons = NamedUnit_O::create(_sym_daltons, "Daltons", "Da");
    unit_daltons->incorporateUnit(unit_grams, 1.0, 1);
    unit_daltons->incorporateUnit(unit_moles, 1.0, -1);
    _sym_daltons->exportYourself()->defparameter(unit_daltons);

#pragma clang diagnostic pop
  } break;
  case pythonClasses: {
  } break;
  case pythonGlobals: {
    // nothing currently
  } break;
  case pythonFunctions: {
    // nothing
  } break;
  };
};

}; // namespace units

extern "C" {
bool init_libunits() {
  static units::UnitsExposer_O *unitsPkgP = NULL;
  if (unitsPkgP == NULL) {
    unitsPkgP = new units::UnitsExposer_O(_lisp);
    _lisp->installPackage(unitsPkgP);
    return true;
  } else
    return false;
}
}
