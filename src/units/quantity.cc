/*
    File: quantity.cc
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
#include <cando/geom/coordinateArray.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/array.h>
#include <cando/units/quantity.h>
#include <clasp/core/wrappers.h>
namespace units {

// ----------------------------------------------------------------------
//

void Quantity_O::throwOnInvalidValue(core::T_sp obj) {
  if (obj.isA<core::Number_O>())
    return;
  if (obj.isA<geom::OVector3_O>())
    return;
  if (obj.isA<core::Array_O>())
    return;
  TYPE_ERROR(obj, core::Cons_O::createList(cl::_sym_or, cl::_sym_Number_O, geom::_sym_OVector3_O, cl::_sym_Array_O));
}

core::T_sp Quantity_O::copyAndScaleValue(core::T_sp obj, double conversion) {
  if (obj.isA<core::Number_O>()) {
    return obj.as_unsafe<core::Number_O>() * core::DoubleFloat_O::create(conversion);
  } else if (obj.isA<geom::OVector3_O>()) {
    geom::OVector3_sp oval = obj.as<geom::OVector3_O>();
    geom::OVector3_sp nval = geom::OVector3_O::create(oval->get().multiplyByScalar(conversion));
    return nval;
  } else if (obj.isA<core::Array_O>()) {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
    core::Array_sp array = obj.as<core::Array_O>();
#pragma clang diagnostic pop
    // core::Array_sp narray = array->deepCopy().as<core::Array_O>();
    FIX_ME();
#if 0
	    narray->multiplyByScalar(conversion);
	    return narray;
#endif
  }
  SIMPLE_ERROR("Illegal value type[{}] for copyAndScaleValue", core::cl__class_of(obj)->_classNameAsString());
}

core::T_sp Quantity_O::copyAndScaleValueElement(core::T_sp obj, uint index, double conversion) {
  if (obj.isA<core::Vector_O>()) {
    core::Vector_sp vec = gc::As_unsafe<core::Vector_sp>(obj);
    core::T_sp element = vec->rowMajorAref(index);
    return Quantity_O::copyAndScaleValue(element, conversion);
  }
  SIMPLE_ERROR("Illegal value type[{}] for copyAndScaleValueElement", core::cl__class_of(obj)->_classNameAsString());
}

core::T_sp Quantity_O::copyValueElement(core::T_sp obj, uint index) {
  if (obj.isA<core::Vector_O>()) {
    //    core::Vector_sp vec = gc::As_unsafe<core::Vector_sp>(obj);
    FIX_ME(); // implement copy of object?
#if 0
          // Why am I copying objects?
          //core::T_sp element = vec->rowMajorAref(index).as<core::General_O>()->deepCopy();
          return element;
#endif
  }
  SIMPLE_ERROR("Illegal value type[{}] for copyAndScaleValueElement", core::cl__class_of(obj)->_classNameAsString());
}

void Quantity_O::setValueElement(core::T_sp obj, uint index, core::T_sp newVal) {
  if (obj.isA<core::Vector_O>()) {
    core::Vector_sp vec = obj.as<core::Vector_O>();
    vec->rowMajorAset(index, newVal);
  } else {
    SIMPLE_ERROR("Illegal value type[{}] for Quantity_O::setValueElement", core::cl__class_of(obj)->_classNameAsString());
  }
}

uint Quantity_O::sizeOfValue(core::T_sp obj) {
  if (obj.isA<core::DoubleFloat_O>()) {
    return 1;
  } else if (obj.isA<geom::OVector3_O>()) {
    return 1;
  } else if (obj.isA<core::Vector_O>()) {
    core::Vector_sp vector = obj.as<core::Vector_O>();
    return vector->length();
  }
  SIMPLE_ERROR("Illegal value type[{}] for sizeOfValue", core::cl__class_of(obj)->_classNameAsString());
}

bool Quantity_O::isnanValue(core::T_sp obj) {
  return obj.isA<core::Float_O>() && core::Float_O::isnan(obj.as_unsafe<core::Float_O>());
}

Quantity_sp Quantity_O::create(core::T_sp obj, Unit_sp unit) {
  Quantity_sp q = Quantity_O::create();
  Quantity_O::throwOnInvalidValue(obj);
  q->_Value = obj;
  q->_Unit = unit;
  return q;
}

Quantity_sp Quantity_O::create(double dbl, Unit_sp unit) {
  Quantity_sp q = Quantity_O::create();
  q->_Value = core::DoubleFloat_O::create(dbl);
  q->_Unit = unit;
  return q;
}

CL_LISPIFY_NAME(make-quantity);
DOCGROUP(cando);
CL_DEFUN Quantity_sp Quantity_O::make(core::T_sp val, Unit_sp unit) {
  // If val is a Number_sp then coerce it to a DoubleFloat_sp
  auto quant = gctools::GC<Quantity_O>::allocate_with_default_constructor();
  if (gc::IsA<core::DoubleFloat_sp>(val)) {
    quant->_Value = gc::As_unsafe<core::DoubleFloat_sp>(val);
  } else if (gc::IsA<core::Number_sp>(val)) {
    core::DoubleFloat_sp dbl = core::clasp_make_double_float(clasp_to_double(val));
    quant->_Value = dbl;
  } else {
    quant->_Value = val;
  }
  quant->_Unit = unit;
  return quant;
}

#ifdef XML_ARCHIVE
void Quantity_O::archiveBase(::core::ArchiveP node) {
  // Archive other instance variables here
  node->archiveObject("value", this->_Value);
  node->archiveObject("unit", this->_Unit);
}
#endif
#ifdef OLD_SERIALIZE
void Quantity_O::serialize(serialize::SNode node) {
  // Archive other instance variables here
  node->archiveObject("value", this->_Value);
  node->archiveObject("unit", this->_Unit);
}
#endif

void Quantity_O::initialize() { this->Base::initialize(); }

CL_LISPIFY_NAME("rawAsString");
CL_DEFMETHOD string Quantity_O::rawAsString() const {
  stringstream ss;
  if (this->_Value.isA<core::DoubleFloat_O>()) {
    ss << this->_Value.as<core::DoubleFloat_O>()->get() * this->_Unit->_Amount;
    ss << "->";
  }
  ss << _rep_(this->_Value);
  ss << "*(";
  ss << _rep_(this->_Unit);
  ss << ")";
  return ss.str();
}

string Quantity_O::__repr__() const {
  stringstream ss;
  core::T_sp valueCopy = Quantity_O::copyAndScaleValue(this->_Value, this->_Unit->_Amount);
  ss << _rep_(valueCopy);
  string unitstr = this->_Unit->unitsOnlyAsString();
  if (unitstr != "") {
    ss << "*" << unitstr;
  }
  return ss.str();
}

void Quantity_O::fields(core::Record_sp node) {
  node->field(INTERN_(units, value), this->_Value);
  node->field(INTERN_(units, unit), this->_Unit);
}

CL_LISPIFY_NAME("quantity_is_compatible");
CL_DEFMETHOD bool Quantity_O::is_compatible(core::T_sp other) const {
  if (other.isA<Unit_O>()) {
    Unit_sp uother = other.as<Unit_O>();
    return this->_Unit->is_compatible(uother);
  } else if (other.isA<Quantity_O>()) {
    Unit_sp uother = other.as<Quantity_O>()->_Unit;
    return this->_Unit->is_compatible(uother);
  }
  SIMPLE_ERROR("You can only check compatibility of a Quantity with another Quantity or Unit - you passed a {}",
               core::cl__class_of(other)->_classNameAsString());
}

CL_LISPIFY_NAME("value_in_unit");
CL_DEFMETHOD core::T_sp Quantity_O::value_in_unit(Unit_sp other, int power) const {
  double conversion = this->_Unit->conversion_factor_to(other, power);
  core::T_sp valueCopy = Quantity_O::copyAndScaleValue(this->_Value, conversion);
  return valueCopy;
}

double Quantity_O::value_in_unit_asReal(Unit_sp other, int power) const {
  core::DoubleFloat_sp val = this->value_in_unit(other, power).as<core::DoubleFloat_O>();
  double dval = val->get();
  return dval;
}

Quantity_sp Quantity_O::getElement(uint i) const {
  core::T_sp val = Quantity_O::copyValueElement(this->_Value, i);
  return Quantity_O::create(val, this->_Unit);
}

CL_LISPIFY_NAME("getElement_in_unit");
CL_DEFMETHOD core::T_sp Quantity_O::getElement_in_unit(uint index, Unit_sp other) const {
  double conversion = this->_Unit->conversion_factor_to(other);
  core::T_sp valueCopy = Quantity_O::copyAndScaleValueElement(this->_Value, index, conversion);
  return valueCopy;
}

CL_LISPIFY_NAME("setElement");
CL_DEFMETHOD core::T_sp Quantity_O::setElement(uint index, Quantity_sp other) {
  // At this point (this) is a Quantity in a particular set of units
  // and other is a Quantity in (hopefully) compatible units
  // I want to convert other into my units and put it into
  // (this) at index (index).
  if (this->_Unit->is_compatible(other->_Unit)) {
    double conversion = other->_Unit->conversion_factor_to(this->_Unit);
    core::T_sp newVal = Quantity_O::copyAndScaleValue(other->_Value, conversion);
    Quantity_O::setValueElement(this->_Value, index, newVal);
  } else {
    SIMPLE_ERROR("Incompatible units - you cannot put an {} with units[{}] into a Quantity with units[{}]",
                 core::cl__class_of(other->_Value)->_classNameAsString(), _rep_(other->_Unit), _rep_(this->_Unit));
  }
  return other;
}

Vector3 Quantity_O::getElement_in_unit_asVector3(uint index, Unit_sp other) const {
  geom::OVector3_sp element = this->getElement_in_unit(index, other).as<geom::OVector3_O>();
  return element->get();
}

CL_LISPIFY_NAME("size");
CL_DEFMETHOD int Quantity_O::size() const { return Quantity_O::sizeOfValue(this->_Value); }

CL_LISPIFY_NAME("*");
CL_DEFMETHOD core::T_sp Quantity_O::operator*(core::T_sp other) const {
  if (gc::IsA<core::Number_sp>(other)) {
    double otherdbl = core::clasp_to_double(other);
    core::T_sp newVal = Quantity_O::copyAndScaleValue(this->_Value, otherdbl);
    return Quantity_O::create(newVal, this->_Unit);
  } else if (other.isA<geom::OVector3_O>()) {
    if (this->_Value.isA<core::Number_O>()) {
      geom::OVector3_sp vother = other.as<geom::OVector3_O>();
      double mydbl = core::clasp_to_double(this->_Value.as<core::Number_O>());
      core::T_sp newVal = Quantity_O::copyAndScaleValue(vother, mydbl);
      return Quantity_O::create(newVal, this->_Unit);
    } else {
      SIMPLE_ERROR(
          ("Currently a Quantity has to be a Number if you are going to multiply it by an OVector3 - the Quantity is a %s"),
          core::cl__class_of(this->_Value)->_classNameAsString());
    }
  } else if (other.isA<Unit_O>()) {
    Unit_sp otherUnit = other.as<Unit_O>();
    double otherdbl = otherUnit->_Amount;
    Quantity_sp q = Quantity_O::create();
    q->_Value = Quantity_O::copyAndScaleValue(this->_Value, otherdbl);
    q->_Unit = Unit_O::create(this->_Unit, 1);
    q->_Unit->incorporateUnit(otherUnit, 1.0, 1);
    return q;
  } else if (other.isA<Quantity_O>()) { // Currently only handle other being a Scalar
    Quantity_sp otherQuant = other.as<Quantity_O>();
    Unit_sp otherUnit = otherQuant->_Unit;
    double otherval = core::clasp_to_double(otherQuant->_Value);
    Quantity_sp q = Quantity_O::create();
    q->_Value = Quantity_O::copyAndScaleValue(this->_Value, otherval);
    q->_Unit = Unit_O::create(this->_Unit, 1);
    q->_Unit->incorporateUnit(otherUnit, 1.0, 1);
    return q;
  }
  TYPE_ERROR(other, core::Cons_O::createList(_sym_Quantity_O, cl::_sym_Number_O));
}

CL_LISPIFY_NAME("/");
CL_DEFMETHOD core::T_sp Quantity_O::operator/(core::T_sp other) const {
  if (gc::IsA<core::Number_sp>(other)) {
    double otherdbl = core::clasp_to_double(other);
    core::T_sp newVal = Quantity_O::copyAndScaleValue(this->_Value, 1.0 / otherdbl);
    return Quantity_O::create(newVal, this->_Unit);
  } else if (other.isA<Unit_O>()) {
    Unit_sp otherUnit = other.as<Unit_O>();
    double otherdbl = otherUnit->_Amount;
    Quantity_sp q = Quantity_O::create();
    q->_Value = Quantity_O::copyAndScaleValue(this->_Value, 1.0 / otherdbl);
    q->_Unit = Unit_O::create(this->_Unit, 1);
    q->_Unit->incorporateUnit(otherUnit, 1.0, -1);
    return q;
  } else if (other.isA<Quantity_O>()) { // Currently only handle other being a Scalar
    Quantity_sp otherQuant = other.as<Quantity_O>();
    Unit_sp otherUnit = otherQuant->_Unit;
    double otherval = core::clasp_to_double(otherQuant->_Value);
    Quantity_sp q = Quantity_O::create();
    q->_Value = Quantity_O::copyAndScaleValue(this->_Value, 1.0 / otherval);
    q->_Unit = Unit_O::create(this->_Unit, 1);
    q->_Unit->incorporateUnit(otherUnit, 1.0, -1);
    return q;
  }
  TYPE_ERROR(other, core::Cons_O::createList(_sym_Quantity_O, cl::_sym_Number_O));
}

CL_LISPIFY_NAME("+");
CL_DEFMETHOD core::T_sp Quantity_O::operator+(core::T_sp other) const {
  if (other.isA<Quantity_O>()) {
    Quantity_sp qother = other.as<Quantity_O>();
    if (this->_Value.isA<core::Number_O>() && qother->_Value.isA<core::Number_O>()) {
      double myValue = core::clasp_to_double(this->_Value.as<core::Number_O>());
      double otherValue = qother->value_in_unit_asReal(this->_Unit);
      double newValue = myValue + otherValue;
      return Quantity_O::create(core::DoubleFloat_O::create(newValue), this->_Unit);
    } else {
      SIMPLE_ERROR("Currently you can only add scalar Quantities - you tried to add a {} to a {}",
                   core::cl__class_of(this->_Value)->_classNameAsString(),
                   core::cl__class_of(qother->_Value)->_classNameAsString());
    }
  }
  TYPE_ERROR(other, core::Cons_O::createList(_sym_Quantity_O, cl::_sym_Number_O));
}

CL_LISPIFY_NAME("-");
CL_DEFMETHOD core::T_sp Quantity_O::operator-(core::T_sp other) const {
  if (other.isA<Quantity_O>()) {
    Quantity_sp qother = other.as<Quantity_O>();
    if (this->_Value.isA<core::Number_O>() && qother->_Value.isA<core::Number_O>()) {
      double myValue = core::clasp_to_double(this->_Value.as<core::Number_O>());
      double otherValue = qother->value_in_unit_asReal(this->_Unit);
      double newValue = myValue - otherValue;
      return Quantity_O::create(core::DoubleFloat_O::create(newValue), this->_Unit);
    } else {
      SIMPLE_ERROR("Currently you can only sub scalar Quantities - you tried to add a {} to a {}",
                   core::cl__class_of(this->_Value)->_classNameAsString(),
                   core::cl__class_of(qother->_Value)->_classNameAsString());
    }
  }
  TYPE_ERROR(other, core::Cons_O::createList(_sym_Quantity_O, cl::_sym_Number_O));
}

CL_LISPIFY_NAME("power");
CL_DEFMETHOD core::T_sp Quantity_O::power(int pwr) const {
  if (this->_Value.isA<core::Number_O>()) {
    double myValue = ::pow(core::clasp_to_double(this->_Value.as<core::Number_O>()), (double)pwr);
    Unit_sp newUnit = Unit_O::create(this->_Unit, pwr);
    return Quantity_O::create(core::DoubleFloat_O::create(myValue), newUnit);
  }
  SIMPLE_ERROR("Currently I can only take powers of Numbers - you tried to take a power of a {}",
               core::cl__class_of(this->_Value)->_classNameAsString());
}

CL_LISPIFY_NAME("sqrt");
CL_DEFMETHOD core::T_sp Quantity_O::sqrt() const {
  if (this->_Value.isA<core::Number_O>()) {
    double myValue = ::sqrt(core::clasp_to_double(this->_Value.as<core::Number_O>()));
    Unit_sp newUnit = Unit_O::createSquareRoot(this->_Unit);
    return Quantity_O::create(core::DoubleFloat_O::create(myValue), newUnit);
  }
  SIMPLE_ERROR("Currently I can only take square root of Numbers - you tried to take square root of a {}",
               core::cl__class_of(this->_Value)->_classNameAsString());
}

CL_LISPIFY_NAME("isnan");
CL_DEFMETHOD bool Quantity_O::isnan() const { return Quantity_O::isnanValue(this->_Value); }

core::T_sp Quantity_O::deepCopy() const {
  /* units are immutable so we don't need to copy them */
  IMPLEMENT_ME();
  // Implement deepCopy
  // return Quantity_O::create(this->_Value.as<core::General_O>()->deepCopy(),this->_Unit);
}

CL_LISPIFY_NAME("<");
CL_DEFMETHOD bool Quantity_O::operator<(core::T_sp other) const {
  if (other.isA<Quantity_O>()) {
    Quantity_sp qother = other.as<Quantity_O>();
    return core::Number_O::compare(gc::As<core::Real_sp>(this->_Value),
                                   gc::As<core::Real_sp>(qother->value_in_unit(this->_Unit, 1))) < 0;
    //	    return this->_Value.as<core::General_O>()->operator<(qother->value_in_unit(this->_Unit,1));
  }
  SIMPLE_ERROR("You cannot compare a Quantity with value of class[{}] to an object of class[{}]",
               core::cl__class_of(this->_Value)->_classNameAsString(), core::cl__class_of(other)->_classNameAsString());
}

CL_LISPIFY_NAME("<=");
CL_DEFMETHOD bool Quantity_O::operator<=(core::T_sp other) const {
  if (other.isA<Quantity_O>()) {
    Quantity_sp qother = other.as<Quantity_O>();
    return core::Number_O::compare(gc::As<core::Real_sp>(this->_Value),
                                   gc::As<core::Real_sp>(qother->value_in_unit(this->_Unit, 1))) <= 0;
    //	    return this->_Value.as<core::General_O>()->operator<=(qother->value_in_unit(this->_Unit,1));
  }
  SIMPLE_ERROR("You cannot compare a Quantity with value of class[{}] to an object of class[{}]",
               core::cl__class_of(this->_Value)->_classNameAsString(), core::cl__class_of(other)->_classNameAsString());
}

CL_LISPIFY_NAME(">");
CL_DEFMETHOD bool Quantity_O::operator>(core::T_sp other) const {
  if (other.isA<Quantity_O>()) {
    Quantity_sp qother = other.as<Quantity_O>();
    return core::Number_O::compare(gc::As<core::Real_sp>(this->_Value),
                                   gc::As<core::Real_sp>(qother->value_in_unit(this->_Unit, 1))) > 0;
    //	    return this->_Value.as<core::General_O>()->operator>(qother->value_in_unit(this->_Unit,1));
  }
  SIMPLE_ERROR("You cannot compare a Quantity with value of class[{}] to an object of class[{}]",
               core::cl__class_of(this->_Value)->_classNameAsString(), core::cl__class_of(other)->_classNameAsString());
}

CL_LISPIFY_NAME(">=");
CL_DEFMETHOD bool Quantity_O::operator>=(core::T_sp other) const {
  if (other.isA<Quantity_O>()) {
    Quantity_sp qother = other.as<Quantity_O>();
    return core::Number_O::compare(gc::As<core::Real_sp>(this->_Value),
                                   gc::As<core::Real_sp>(qother->value_in_unit(this->_Unit, 1))) >= 0;
    //	    return this->_Value.as<core::General_O>()->operator>=(qother->value_in_unit(this->_Unit,1));
  }
  SIMPLE_ERROR("You cannot compare a Quantity with value of class[{}] to an object of class[{}]",
               core::cl__class_of(this->_Value)->_classNameAsString(), core::cl__class_of(other)->_classNameAsString());
}

}; // namespace units
