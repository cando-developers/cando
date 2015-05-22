#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
//#include "core/serialize.h"
#include <cando/candoBase/coordinateArray.h>
#include <cando/candoBase/ovector3.h>
#include <clasp/core/array.h>
#include <cando/units/quantity.h>
#include <clasp/core/wrappers.h>
namespace units
{

// ----------------------------------------------------------------------
//

    

    void Quantity_O::throwOnInvalidValue(core::T_sp obj,core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::Number_O>() ) return;
	if ( obj->isAssignableTo<candoBase::OVector3_O>() ) return;
	if ( obj->isAssignableTo<core::Array_O>() ) return;
	SIMPLE_ERROR(BF("Illegal value type[%s] for Quantity") % obj->__class()->classNameAsString() );
    }


    core::T_sp Quantity_O::copyAndScaleValue(core::T_sp obj, double conversion, core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::DoubleFloat_O>() )
	{
	    core::DoubleFloat_sp rval = core::DoubleFloat_O::create(obj.as<core::DoubleFloat_O>()->get()*conversion);
	    return rval;
	} else if ( obj->isAssignableTo<candoBase::OVector3_O>() )
	{
	    candoBase::OVector3_sp oval = obj.as<candoBase::OVector3_O>();
	    candoBase::OVector3_sp nval = lisp->create<candoBase::OVector3_O>(oval->get().multiplyByScalar(conversion));
	    return nval;
	} else if ( obj->isAssignableTo<core::Array_O>() )
	{
	    core::Array_sp array = obj.as<core::Array_O>();
	    core::Array_sp narray = array->deepCopy().as<core::Array_O>();
	    narray->multiplyByScalar(conversion);
	    return narray;
	}
	SIMPLE_ERROR(BF("Illegal value type[%s] for copyAndScaleValue") % obj->__class()->classNameAsString() );
    }
	

    core::T_sp Quantity_O::copyAndScaleValueElement(core::T_sp obj, uint index,
						   double conversion, core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::Vector_O>() )
	{
	    core::Vector_sp vec = obj.as<core::Vector_O>();
	    core::T_sp element = vec->elt(index);
	    return Quantity_O::copyAndScaleValue(element,conversion,lisp);
	}
	SIMPLE_ERROR(BF("Illegal value type[%s] for copyAndScaleValueElement") % obj->__class()->classNameAsString() );
    }


    core::T_sp Quantity_O::copyValueElement(core::T_sp obj, uint index,
					   core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::Vector_O>() )
	{
	    core::Vector_sp vec = obj.as<core::Vector_O>();
	    core::T_sp element = vec->elt(index)->deepCopy();
	    return element;
	}
	SIMPLE_ERROR(BF("Illegal value type[%s] for copyAndScaleValueElement") % obj->__class()->classNameAsString() );
    }



    void Quantity_O::setValueElement(core::T_sp obj, uint index,
				     core::T_sp newVal, core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::Vector_O>() )
	{
	    core::Vector_sp vec = obj.as<core::Vector_O>();
	    vec->setf_elt(index,newVal);
	} else
	{
	    SIMPLE_ERROR(BF("Illegal value type[%s] for Quantity_O::setValueElement") % obj->__class()->classNameAsString() );
	}
    }


    uint Quantity_O::sizeOfValue(core::T_sp obj, core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::DoubleFloat_O>() )
	{
	    return 1;
	} else if ( obj->isAssignableTo<candoBase::OVector3_O>() )
	{
	    return 1;
	} else if ( obj->isAssignableTo<core::Vector_O>() )
	{
	    core::Vector_sp vector = obj.as<core::Vector_O>();
	    return vector->length();
	}
	SIMPLE_ERROR(BF("Illegal value type[%s] for sizeOfValue") % obj->__class()->classNameAsString() );
    }
	


    bool Quantity_O::isnanValue(core::T_sp obj, core::Lisp_sp lisp)
    {_G();
	if ( obj->isAssignableTo<core::DoubleFloat_O>() )
	{
	    return obj.as<core::DoubleFloat_O>()->isnan();
#if 0
	} else if ( obj->isAssignableTo<candoBase::OVector3_O>() )
	{
	    return
	} else if ( obj->isAssignableTo<core::Vector_O>() )
	{
	    core::Vector_sp vector = obj.as<core::Vector_O>();
	    return vector->length();
#endif
	}
	SIMPLE_ERROR(BF("Illegal value type[%s] for isnanValue") % obj->__class()->classNameAsString() );
    }
	


    Quantity_sp Quantity_O::create( core::T_sp obj, Unit_sp unit, core::Lisp_sp lisp)
    {_G();
	Quantity_sp q = lisp->create<Quantity_O>();
	Quantity_O::throwOnInvalidValue(obj,_lisp);
	q->_Value = obj;
	q->_Unit = unit;
	return q;
    }

    Quantity_sp Quantity_O::create( double dbl, Unit_sp unit, core::Lisp_sp lisp)
    {_G();
	Quantity_sp q = lisp->create<Quantity_O>();
	q->_Value = core::DoubleFloat_O::create(dbl);
	q->_Unit = unit;
	return q;
    }




    Quantity_sp Quantity_O::make(core::T_sp val, Unit_sp unit)
    {_G();
	GC_ALLOCATE(Quantity_O,quant);
	quant->_Value = val;
	quant->_Unit = unit;
	return quant;
    }

#ifdef XML_ARCHIVE
    void Quantity_O::archiveBase(::core::ArchiveP node)
    {
	// Archive other instance variables here
	node->archiveObject("value",this->_Value);
	node->archiveObject("unit",this->_Unit);
    }
#endif
#ifdef OLD_SERIALIZE
    void Quantity_O::serialize(serialize::SNode node)
    {
	// Archive other instance variables here
	node->archiveObject("value",this->_Value);
	node->archiveObject("unit",this->_Unit);
    }
#endif

    void Quantity_O::initialize()
    {_OF();
        this->Base::initialize();
    }


    string Quantity_O::rawAsString() const
    {_OF();
	stringstream ss;
	if ( this->_Value->isAssignableTo<core::DoubleFloat_O>() )
	{
	    ss << this->_Value.as<core::DoubleFloat_O>()->get()*this->_Unit->_Amount;
	    ss << "->";
	}
	ss << this->_Value->__repr__();
	ss << "*(";
	ss << this->_Unit->__repr__();
	ss << ")";
	return ss.str();
    }


    string Quantity_O::__repr__() const
    {_OF();
	stringstream ss;
	core::T_sp valueCopy = Quantity_O::copyAndScaleValue(this->_Value,this->_Unit->_Amount,_lisp);
	ss << valueCopy->__repr__();
	string unitstr = this->_Unit->unitsOnlyAsString();
	if ( unitstr != "" )
	{
	    ss << "*" << unitstr;
	}
	return ss.str();
    }


    bool Quantity_O::is_compatible(core::T_sp other) const
    {_OF();
	if ( other->isAssignableTo<Unit_O>() )
	{
	    Unit_sp uother = other.as<Unit_O>();
	    return this->_Unit->is_compatible(uother);
	} else if ( other->isAssignableTo<Quantity_O>() )
	{
	    Unit_sp uother = other.as<Quantity_O>()->_Unit;
	    return this->_Unit->is_compatible(uother);
	}
	SIMPLE_ERROR(BF("You can only check compatibility of a Quantity with another Quantity or Unit - you passed a %s") % other->__class()->classNameAsString() );
    }

    core::T_sp Quantity_O::value_in_unit(Unit_sp other, int power) const
    {_OF();
	double conversion = this->_Unit->conversion_factor_to(other,power);
	core::T_sp valueCopy = Quantity_O::copyAndScaleValue(this->_Value,conversion,_lisp);
	return valueCopy;
    }


    double Quantity_O::value_in_unit_asReal(Unit_sp other, int power) const
    {_OF();
	core::DoubleFloat_sp val = this->value_in_unit(other,power).as<core::DoubleFloat_O>();
	double dval = val->get();
	return dval;
    }


    Quantity_sp Quantity_O::getElement(uint i) const
    {_OF();
	core::T_sp val = Quantity_O::copyValueElement(this->_Value,i,_lisp);
	return Quantity_O::create(val,this->_Unit,_lisp);
    }

    core::T_sp Quantity_O::getElement_in_unit(uint index, Unit_sp other ) const
    {_OF();
	double conversion = this->_Unit->conversion_factor_to(other);
	core::T_sp valueCopy = Quantity_O::copyAndScaleValueElement(this->_Value,index,conversion,_lisp);
	return valueCopy;
    }

    core::T_sp Quantity_O::setElement(uint index, Quantity_sp other)
    {_OF();
	// At this point (this) is a Quantity in a particular set of units
	// and other is a Quantity in (hopefully) compatible units
	// I want to convert other into my units and put it into
	// (this) at index (index).
	if ( this->_Unit->is_compatible(other->_Unit) )
	{
	    double conversion = other->_Unit->conversion_factor_to(this->_Unit);
	    core::T_sp newVal = Quantity_O::copyAndScaleValue(other->_Value,conversion,_lisp);
	    Quantity_O::setValueElement(this->_Value,index,newVal,_lisp);
	} else
	{
	    SIMPLE_ERROR(BF("Incompatible units - you cannot put an %s with units[%s] into a Quantity with units[%s]") % other->_Value->__class()->classNameAsString() % other->_Unit->__repr__() % this->_Unit->__repr__() );
	}
	return other;
    }

    Vector3 Quantity_O::getElement_in_unit_asVector3(uint index, Unit_sp other) const
    {_OF();
	candoBase::OVector3_sp element = this->getElement_in_unit(index,other).as<candoBase::OVector3_O>();
	return element->get();
    }


    int Quantity_O::size() const
    {_OF();
	return Quantity_O::sizeOfValue(this->_Value,_lisp);
    }




    core::T_sp Quantity_O::operator*(core::T_sp other) const
    {_OF();
	if ( other.nilp() ) return this->const_sharedThis<Quantity_O>();
	if ( other->isAssignableTo<core::DoubleFloat_O>() )
	{
	    double otherdbl = other.as<core::DoubleFloat_O>()->get();
	    core::T_sp newVal = Quantity_O::copyAndScaleValue(this->_Value,otherdbl,_lisp);
	    return Quantity_O::create(newVal,this->_Unit,_lisp);
	} else if ( other->isAssignableTo<candoBase::OVector3_O>() )
	{
	    if ( this->_Value->isAssignableTo<core::Number_O>() )
	    {
		candoBase::OVector3_sp vother = other.as<candoBase::OVector3_O>();
		double mydbl = this->_Value.as<core::Number_O>()->as_double();
		core::T_sp newVal = Quantity_O::copyAndScaleValue(vother,mydbl,_lisp);
		return Quantity_O::create(newVal,this->_Unit,_lisp);
	    } else
	    {
		SIMPLE_ERROR(BF("Currently a Quantity has to be a Number if you are going to multiply it by an OVector3 - the Quantity is a %s") % this->_Value->__class()->classNameAsString() );
	    }
	} else if ( other->isAssignableTo<Unit_O>() )
	{
	    Unit_sp otherUnit = other.as<Unit_O>();
	    double otherdbl = otherUnit->_Amount;
	    Quantity_sp q = Quantity_O::create();
	    q->_Value = Quantity_O::copyAndScaleValue(this->_Value,otherdbl,_lisp);
	    q->_Unit = Unit_O::create(this->_Unit,1,_lisp);
	    q->_Unit->incorporateUnit(otherUnit,1.0,1);
	    return q;
	} else if ( other->isAssignableTo<Quantity_O>() )
	{	// Currently only handle other being a Scalar
	    Quantity_sp otherQuant = other.as<Quantity_O>();
	    Unit_sp otherUnit = otherQuant->_Unit;
	    double otherval = otherQuant->_Value.as<core::DoubleFloat_O>()->get();
	    Quantity_sp q = Quantity_O::create();
	    q->_Value = Quantity_O::copyAndScaleValue(this->_Value,otherval,_lisp);
	    q->_Unit = Unit_O::create(this->_Unit,1,_lisp);
	    q->_Unit->incorporateUnit(otherUnit,1.0,1);
	    return q;
	}
	SIMPLE_ERROR(BF("I cannot yet handle * of a Quantity with a %s") % other->__class()->classNameAsString() );
    }




    core::T_sp Quantity_O::operator/(core::T_sp other) const
    {_OF();
	if ( other.nilp() ) return this->const_sharedThis<Quantity_O>();
	if ( other->isAssignableTo<core::DoubleFloat_O>() )
	{
	    double otherdbl = other.as<core::DoubleFloat_O>()->get();
	    core::T_sp newVal = Quantity_O::copyAndScaleValue(this->_Value,1.0/otherdbl,_lisp);
	    return Quantity_O::create(newVal,this->_Unit,_lisp);
	} else if ( other->isAssignableTo<Unit_O>() )
	{
	    Unit_sp otherUnit = other.as<Unit_O>();
	    double otherdbl = otherUnit->_Amount;
	    Quantity_sp q = Quantity_O::create();
	    q->_Value = Quantity_O::copyAndScaleValue(this->_Value,1.0/otherdbl,_lisp);
	    q->_Unit = Unit_O::create(this->_Unit,1,_lisp);
	    q->_Unit->incorporateUnit(otherUnit,1.0,-1);
	    return q;
	} else if ( other->isAssignableTo<Quantity_O>() )
	{	// Currently only handle other being a Scalar
	    Quantity_sp otherQuant = other.as<Quantity_O>();
	    Unit_sp otherUnit = otherQuant->_Unit;
	    double otherval = otherQuant->_Value.as<core::DoubleFloat_O>()->get();
	    Quantity_sp q = Quantity_O::create();
	    q->_Value = Quantity_O::copyAndScaleValue(this->_Value,1.0/otherval,_lisp);
	    q->_Unit = Unit_O::create(this->_Unit,1,_lisp);
	    q->_Unit->incorporateUnit(otherUnit,1.0,-1);
	    return q;
	}
	SIMPLE_ERROR(BF("I cannot yet handle * of a Quantity with a %s") % other->__class()->classNameAsString() );
    }


    core::T_sp Quantity_O::operator+(core::T_sp other) const
    {_OF();
	if ( other.nilp() ) return this->const_sharedThis<Quantity_O>();
	if ( other->isAssignableTo<Quantity_O>() )
	{
	    Quantity_sp qother = other.as<Quantity_O>();
	    if ( this->_Value->isAssignableTo<core::Number_O>()
		 && qother->_Value->isAssignableTo<core::Number_O>() )
	    {
		double myValue = this->_Value.as<core::Number_O>()->as_double();
		double otherValue = qother->value_in_unit_asReal(this->_Unit);
		double newValue = myValue + otherValue;
		return Quantity_O::create(core::DoubleFloat_O::create(newValue),
					  this->_Unit,_lisp);
	    }
	    else
	    {
		SIMPLE_ERROR(BF("Currently you can only add scalar Quantities - you tried to add a %s to a %s" ) % this->_Value->__class()->classNameAsString() % qother->_Value->__class()->classNameAsString() );
	    }
	} else
	{
	    SIMPLE_ERROR(BF("You tried to add a Quantity to a %s") % other->__class()->classNameAsString() );
	}
    }


    core::T_sp Quantity_O::operator-(core::T_sp other) const
    {_OF();
	if ( other.nilp() ) return this->const_sharedThis<Quantity_O>();
	if ( other->isAssignableTo<Quantity_O>() )
	{
	    Quantity_sp qother = other.as<Quantity_O>();
	    if ( this->_Value->isAssignableTo<core::Number_O>()
		 && qother->_Value->isAssignableTo<core::Number_O>() )
	    {
		double myValue = this->_Value.as<core::Number_O>()->as_double();
		double otherValue = qother->value_in_unit_asReal(this->_Unit);
		double newValue = myValue - otherValue;
		return Quantity_O::create(core::DoubleFloat_O::create(newValue),
					  this->_Unit,_lisp);
	    }
	    else
	    {
		SIMPLE_ERROR(BF("Currently you can only sub scalar Quantities - you tried to add a %s to a %s" ) % this->_Value->__class()->classNameAsString() % qother->_Value->__class()->classNameAsString() );
	    }
	} else
	{
	    SIMPLE_ERROR(BF("You tried to subtract from a Quantity a %s - not currently allowed") % other->__class()->classNameAsString() );
	}
    }




    core::T_sp Quantity_O::power(int pwr) const
    {_OF();
	if ( this->_Value->isAssignableTo<core::Number_O>() )
	{
	    double myValue = ::pow(this->_Value.as<core::Number_O>()->as_double(),(double)pwr);
	    Unit_sp newUnit = Unit_O::create(this->_Unit,pwr,_lisp);
	    return Quantity_O::create(core::DoubleFloat_O::create(myValue),newUnit,_lisp);
	}
	SIMPLE_ERROR(BF("Currently I can only take powers of Numbers - you tried to take a power of a %s") % this->_Value->__class()->classNameAsString());
    }

    core::T_sp Quantity_O::sqrt() const
    {_OF();
	if ( this->_Value->isAssignableTo<core::Number_O>() )
	{
	    double myValue = ::sqrt(this->_Value.as<core::Number_O>()->as_double());
	    Unit_sp newUnit = Unit_O::createSquareRoot(this->_Unit,_lisp);
	    return Quantity_O::create(core::DoubleFloat_O::create(myValue),newUnit,_lisp);
	}
	SIMPLE_ERROR(BF("Currently I can only take square root of Numbers - you tried to take square root of a %s") % this->_Value->__class()->classNameAsString());
    }


    bool Quantity_O::isnan() const
    {_OF();
	return Quantity_O::isnanValue(this->_Value,_lisp);
    }


    core::T_sp Quantity_O::deepCopy() const
    {_OF();
	/* units are immutable so we don't need to copy them */
	return Quantity_O::create(this->_Value->deepCopy(),this->_Unit,_lisp);
    }





    bool Quantity_O::operator<(core::T_sp other) const
    {_OF();
	if ( other->isAssignableTo<Quantity_O>() )
	{
	    Quantity_sp qother = other.as<Quantity_O>();
	    return this->_Value->operator<(qother->value_in_unit(this->_Unit,1));
	}
	SIMPLE_ERROR(BF("You cannot compare a Quantity with value of class[%s] to an object of class[%s]") % this->_Value->__class()->classNameAsString() % other->__class()->classNameAsString() );
    }

    bool Quantity_O::operator<=(core::T_sp other) const
    {_OF();
	if ( other->isAssignableTo<Quantity_O>() )
	{
	    Quantity_sp qother = other.as<Quantity_O>();
	    return this->_Value->operator<=(qother->value_in_unit(this->_Unit,1));
	}
	SIMPLE_ERROR(BF("You cannot compare a Quantity with value of class[%s] to an object of class[%s]") % this->_Value->__class()->classNameAsString() % other->__class()->classNameAsString() );
    }


    bool Quantity_O::operator>(core::T_sp other) const
    {_OF();
	if ( other->isAssignableTo<Quantity_O>() )
	{
	    Quantity_sp qother = other.as<Quantity_O>();
	    return this->_Value->operator>(qother->value_in_unit(this->_Unit,1));
	}
	SIMPLE_ERROR(BF("You cannot compare a Quantity with value of class[%s] to an object of class[%s]") % this->_Value->__class()->classNameAsString() % other->__class()->classNameAsString() );
    }

    bool Quantity_O::operator>=(core::T_sp other) const
    {_OF();
	if ( other->isAssignableTo<Quantity_O>() )
	{
	    Quantity_sp qother = other.as<Quantity_O>();
	    return this->_Value->operator>=(qother->value_in_unit(this->_Unit,1));
	}
	SIMPLE_ERROR(BF("You cannot compare a Quantity with value of class[%s] to an object of class[%s]") % this->_Value->__class()->classNameAsString() % other->__class()->classNameAsString() );
    }





    EXPOSE_CLASS(units,Quantity_O);

    void Quantity_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<Quantity_O>()
//	.initArgs("(self)")
//	    .def_raw("__init__",&Quantity_O::__init__,"(self units::value units::unit)")
	    .def("is_compatible",&Quantity_O::is_compatible)
	    .def("value_in_unit",&Quantity_O::value_in_unit)
	    .def("setElement",&Quantity_O::setElement)
	    .def("getElement_in_unit",&Quantity_O::getElement_in_unit)
	    .def("rawValue",&Quantity_O::rawValue)
	    .def("rawUnit",&Quantity_O::rawUnit)
	    .def("rawAsString",&Quantity_O::rawAsString)
	    .def("size",&Quantity_O::size)
	    .def("*",&Quantity_O::operator*)
	    .def("/",&Quantity_O::operator/)
	    .def("+",&Quantity_O::operator+)
	    .def("-",&Quantity_O::operator-)
	    .def("power",&Quantity_O::power)
	    .def("sqrt",&Quantity_O::sqrt)
	    .def("isnan",&Quantity_O::isnan)
	    .def("<",&Quantity_O::operator<)
	    .def("<=",&Quantity_O::operator<=)
	    .def(">",&Quantity_O::operator>)
	    .def(">=",&Quantity_O::operator>=)
	;
    }



    void Quantity_O::exposePython(::core::Lisp_sp lisp)
    {
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(UnitsPkg,Quantity,"","",lisp)
//	.initArgs("(self)")
	    .def("is_compatible",&Quantity_O::is_compatible)
	    .def("value_in_unit",&Quantity_O::value_in_unit)
	    .def("getElement_in_unit",&Quantity_O::getElement_in_unit)
	    .def("setElement",&Quantity_O::setElement)
	    .def("rawValue",&Quantity_O::rawValue)
	    .def("rawUnit",&Quantity_O::rawUnit)
	    .def("rawAsString",&Quantity_O::rawAsString)
	    .def("size",&Quantity_O::size)
	    .def("*",&Quantity_O::operator*)
	    .def("/",&Quantity_O::operator/)
	    .def("+",&Quantity_O::operator+)
	    .def("-",&Quantity_O::operator-)
	    .def("power",&Quantity_O::power)
	    .def("sqrt",&Quantity_O::sqrt)
	    .def("isnan",&Quantity_O::isnan)
	    .def("<",&Quantity_O::operator<)
	    .def("<=",&Quantity_O::operator<=)
	    .def(">",&Quantity_O::operator>)
	    .def(">=",&Quantity_O::operator>=)
	;
#endif
    }







}; /* units */
