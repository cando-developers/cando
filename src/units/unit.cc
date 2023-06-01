/*
    File: unit.cc
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

#include <clasp/core/common.h>
#include <cando/units/dimension.fwd.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/array.h>
#include <cando/units/quantity.h>
#include <cando/units/unit.h>
#include <clasp/core/wrappers.h>
#include <clasp/core/num_co.h> // cl:numerator, cl:denominator

namespace units {

// ----------------------------------------------------------------------
//

Unit_sp Unit_O::create(Unit_sp orig, int power) {
  Unit_sp unit = Unit_O::create();
  unit->incorporateUnit(orig, 1.0, power);
  return unit;
}

Unit_sp Unit_O::createSquareRoot(Unit_sp orig) {
  Unit_sp unit = Unit_O::create();
  unit->_Amount = ::sqrt(orig->_Amount);
  for (int i = 0; i < NumBaseDimensions; i++) {
    if (unit->_Powers[i]._Num % 2 == 0) {
      unit->_Powers[i]._Num /= 2;
    } else {
      unit->_Powers[i]._Denom *= 2;
    }
    unit->_Powers[i].maybe_simplify();
  }
  return unit;
}

/*! Create a new unit by combining other units with powers.
eg: (make-unit (list units:meters 1 units:seconds -1))
*/
CL_LISPIFY_NAME(make-unit);
DOCGROUP(cando);
CL_DEFUN Unit_sp Unit_O::make(core::List_sp args) {
  auto u = gctools::GC<Unit_O>::allocate_with_default_constructor();
  Unit_O::parseUnitList(u->_Amount, u->_Powers, args);
  return u;
}

void Unit_O::parseUnitList(double &amountScale, SimpleRational powers[NumBaseDimensions], core::List_sp list) {
  //	core::Binder_sp unitDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
  core::List_sp dimCur = list;
  while (dimCur.notnilp()) {
    core::T_sp odim = core::oCar(dimCur);
    if (odim.isA<Unit_O>()) {
      Unit_sp u = odim.as<Unit_O>();
      int power = 1;
      if (core::oCadr(dimCur).isA<core::Rational_O>()) {
        power = core::clasp_to_fixnum(core::oCadr(dimCur));
        LOG("Got power[{}]", power);
        dimCur = core::oCdr(dimCur);
      }
      u->adjustPowersAndAmountScale(power, powers, amountScale);
    } else {
      SIMPLE_ERROR("Unknown unit[{}]", _rep_(odim));
    }
    dimCur = core::oCdr(dimCur);
  }
  LOG("The amountScale[{}]", amountScale);
}

void Unit_O::fields(core::Record_sp node) {
#define UFIELDS() \
  node->field(INTERN_(units, amount), amt);\
  node->field(INTERN_(units, meters), meters);\
  node->field(INTERN_(units, kilograms), kilograms);\
  node->field(INTERN_(units, seconds), seconds);\
  node->field(INTERN_(units, amperes), amperes);\
  node->field(INTERN_(units, kelvin), kelvin);\
  node->field(INTERN_(units, candelas), candelas);\
  node->field(INTERN_(units, moles), moles);\
  node->field(INTERN_(units, radians), radians);
  switch (node->stage()) {
  case core::Record_O::initializing: [[fallthrough]];
  case core::Record_O::loading: {
    core::DoubleFloat_sp amt;
    core::Rational_sp meters, kilograms, seconds, amperes, kelvin, candelas, moles, radians;
    UFIELDS();
    this->_Amount = amt->as_double_();
#define RAT(rat) SimpleRational(core::clasp_to_integral<int>(cl__numerator(rat)), core::clasp_to_integral<int>(cl__denominator(rat)))
    this->_Powers[0] = RAT(meters);
    this->_Powers[1] = RAT(kilograms);
    this->_Powers[2] = RAT(seconds);
    this->_Powers[3] = RAT(amperes);
    this->_Powers[4] = RAT(kelvin);
    this->_Powers[5] = RAT(candelas);
    this->_Powers[6] = RAT(moles);
    this->_Powers[7] = RAT(radians);
#undef RAT
    break;
  }
  case core::Record_O::saving: {
#define RAT(dim) core::Rational_O::create(core::Integer_O::create((dim)._Num), core::Integer_O::create((dim)._Denom))
    core::DoubleFloat_sp amt = core::DoubleFloat_O::create(this->_Amount);
    core::Rational_sp meters = RAT(this->_Powers[0]);
    core::Rational_sp kilograms = RAT(this->_Powers[1]);
    core::Rational_sp seconds = RAT(this->_Powers[2]);
    core::Rational_sp amperes = RAT(this->_Powers[3]);
    core::Rational_sp kelvin = RAT(this->_Powers[4]);
    core::Rational_sp candelas = RAT(this->_Powers[5]);
    core::Rational_sp moles = RAT(this->_Powers[6]);
    core::Rational_sp radians = RAT(this->_Powers[7]);
    UFIELDS();
#undef RAT
    break;
  }
  default:
      IMPLEMENT_MEF("Add support to patch units");
  }
#undef UFIELDS
}

void Unit_O::initialize() {
  this->Base::initialize();
  Dimension_O::zeroPowers(this->_Powers);
  this->_Amount = 1.0;
}

void Unit_O::incorporateUnit(Unit_sp u, double amountScale, int power) {
  /*Set the default system as units:*SI* */
  //	core::Binder_sp unitDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
  // Repeat this block for multiple symbols
  {
    u->adjustPowersAndAmountScale(power, this->_Powers, amountScale);
    this->_Amount *= amountScale;
  }
}

void Unit_O::adjustPowersAndAmountScale(SimpleRational power, SimpleRational powers[], double &amountScale) const {
  for (int id = 0; id < NumBaseDimensions; id++) {
    powers[id] = powers[id] + this->_Powers[id] * power;
  }
  double scale = ::pow(this->_Amount, (double)power._Num/(double)power._Denom);
  amountScale *= scale;
}

string Unit_O::__repr__() const {
  stringstream ss;
  ss << this->_Amount << "*" << this->unitsOnlyAsString();
  return ss.str();
}

void Unit_O::__write__(core::T_sp stream) const {
  core::clasp_write_string("#<UNIT ", stream);
  core::clasp_write_string(this->__repr__(), stream);
  core::clasp_write_char('>', stream);
}

CL_LISPIFY_NAME("test_set_amount");
CL_DEFMETHOD void Unit_O::test_set_amount(double amount) { this->_Amount = amount; }

string Unit_O::unitsOnlyAsString() const {
  stringstream units;
  for (int i = 0; i < NumBaseDimensions; i++) {
    if (this->_Powers[i]._Num != 0) {
      units << Dimension_O::baseDimensionUnitName(i);
      if (this->_Powers[i]._Num != this->_Powers[i]._Denom) {
        units << "^" << this->_Powers[i]._Num;
        if (this->_Powers[i]._Denom>1) {
          units<< "/" << this->_Powers[i]._Denom;
        }
      }
      units << "*";
    }
  }
  return units.str().substr(0, units.str().size() - 1);
};

CL_LISPIFY_NAME("unit_is_compatible");
CL_DEFMETHOD bool Unit_O::is_compatible(Unit_sp other, int power) const {
  for (int i = 0; i < NumBaseDimensions; i++) {
    if (this->_Powers[i] != other->_Powers[i] * power)
      return false;
  }
  return true;
}

CL_LISPIFY_NAME("conversion_factor_to");
CL_DEFMETHOD double Unit_O::conversion_factor_to(Unit_sp other, int power) const {
  if (this->is_compatible(other, power)) {
    double conversion = this->_Amount / ::pow(other->_Amount, power);
    return conversion;
  }
  SIMPLE_ERROR("Units[{}] are not compatible with Unit[{}]^{}", this->__repr__(), _rep_(other), power);
}

Unit_sp Unit_O::copyWithoutName() const {
  Unit_sp result = Unit_O::create();
  Dimension_O::copyPowers(result->_Powers, this->_Powers);
  result->_Amount = this->_Amount;
  return result;
}

CL_LISPIFY_NAME("*");
CL_DEFMETHOD core::T_sp Unit_O::operator*(core::T_sp obj) const {
  core::T_sp result;
  if (obj.isA<Unit_O>()) {
    Unit_sp other = obj.as<Unit_O>();
    Unit_sp uresult = this->copyWithoutName();
    for (int i = 0; i < NumBaseDimensions; i++) {
      uresult->_Powers[i] = uresult->_Powers[i] + other->_Powers[i];
    }
    uresult->_Amount *= other->_Amount;
    result = uresult;
  } else if (obj.isA<core::Number_O>() || obj.isA<geom::OVector3_O>() || obj.isA<core::Array_O>()) {
    result = Quantity_O::create(obj, this->const_sharedThis<Unit_O>());
  } else if (obj.isA<Quantity_O>()) {
    SIMPLE_ERROR("Handle Unit*Quantity");
  } else {
    SIMPLE_ERROR("Handle Unit*XXX where XXX={}", core::cl__class_of(obj)->_classNameAsString());
  }
  return result;
}

CL_LISPIFY_NAME("/");
CL_DEFMETHOD core::T_sp Unit_O::operator/(core::T_sp obj) const {
  if (obj.isA<Unit_O>()) {
    Unit_sp other = obj.as<Unit_O>();
    Unit_sp result = this->copyWithoutName();
    for (int i = 0; i < NumBaseDimensions; i++) {
      result->_Powers[i] = result->_Powers[i] - other->_Powers[i];
    }
    result->_Amount /= other->_Amount;
    return result;
  } else {
    SIMPLE_ERROR("Handle Unit/XXX where XXX={}", core::cl__class_of(obj)->_classNameAsString());
  }
}

}; // namespace units
