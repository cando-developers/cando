#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
//#include "core/serialize.h"
#include <cando/units/dimension.fwd.h>
#include <cando/geom/ovector3.h>
#include <clasp/core/array.h>
#include <cando/units/quantity.h>
#include <cando/units/unit.h>
#include <clasp/core/wrappers.h>


namespace units
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(units,Unit_O);

    void Unit_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<Unit_O>()
//	    .def_raw("__init__",&Unit_O::__init__,"(self &rest args)")
	    .def("test_set_amount",&Unit_O::test_set_amount)
	    .def("unit_is_compatible",&Unit_O::is_compatible)
	    .def("conversion_factor_to",&Unit_O::conversion_factor_to)
	    .def("*",&Unit_O::operator*)
	    .def("/",&Unit_O::operator/)
	;
      core::af_def(UnitsPkg,"make-unit",&Unit_O::make);
    }

    void Unit_O::exposePython(core::Lisp_sp lisp)
    {_G();
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(UnitsPkg,Unit,"","",_lisp)
//	.def_raw("__init__",&initArgs("(self)")
	;
#endif
    }



Unit_sp Unit_O::create(Unit_sp orig, int power)
    {_G();
      Unit_sp unit = Unit_O::create();
	unit->incorporateUnit(orig,1.0,power);
	return unit;
    }

Unit_sp Unit_O::createSquareRoot(Unit_sp orig)
    {_G();
      Unit_sp unit = Unit_O::create();
	unit->_Amount = ::sqrt(orig->_Amount);
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    if ( orig->_Powers[i] % 2 != 0 )
	    {
		SIMPLE_ERROR(BF("One of the powers is not divisible by 2"));
	    }
	    unit->_Powers[i] = orig->_Powers[i]/2;
	}
	return unit;
    }


/*! Create a new unit by combining other units with powers.
eg: (make-unit (list units:meters 1 units:seconds -1))
*/
    Unit_sp Unit_O::make(core::List_sp args)
    {_G();
	GC_ALLOCATE(Unit_O,u);
	Unit_O::parseUnitList(u->_Amount,u->_Powers,args);
	return u;
    }


    void Unit_O::parseUnitList(double& amountScale, int powers[NumBaseDimensions], core::List_sp list)
    {_G();
	LOG(BF("Creating unit with arguments[%s]") % lisp->__repr__() );
//	core::Binder_sp unitDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
	core::List_sp dimCur = list;
	while (dimCur.notnilp())
	{
	    core::T_sp odim = core::oCar(dimCur);
	    if ( odim.isA<Unit_O>() )
	    {
		Unit_sp u = odim.as<Unit_O>();
		int power = 1;
		if ( core::oCadr(dimCur).isA<core::Rational_O>() )
		{
                  power = core::clasp_to_fixnum(core::oCadr(dimCur).as<core::Rational_O>());
		    LOG(BF("Got power[%d]") % power);
		    dimCur = core::oCdr(dimCur);
		}
		u->adjustPowersAndAmountScale(power,powers,amountScale);
	    } else
	    {
		SIMPLE_ERROR(BF("Unknown unit[%s]") % odim->__repr__());
	    }
	    dimCur = core::oCdr(dimCur);
	}
	LOG(BF("The amountScale[%lf]") % amountScale );
    }



#ifdef OLD_SERIALIZE
    void Unit_O::serialize(serialize::SNode node)
    {_G();
	// Archive other instance variables here
	ASSERTF(NumBaseDimensions==8,BF("You changed NumBaseDimensions but didnt change Unit_O::archiveBase"));
	node->attribute("amt",this->_Amount);
	node->attributeIfNotDefault("d0",this->_Powers[0],0);
	node->attributeIfNotDefault("w1",this->_Powers[1],0);
	node->attributeIfNotDefault("t2",this->_Powers[2],0);
	node->attributeIfNotDefault("c3",this->_Powers[3],0);
	node->attributeIfNotDefault("t4",this->_Powers[4],0);
	node->attributeIfNotDefault("l5",this->_Powers[5],0);
	node->attributeIfNotDefault("a6",this->_Powers[6],0);
	node->attributeIfNotDefault("r7",this->_Powers[7],0);
    }
#endif

#ifdef XML_ARCHIVE
    void Unit_O::archiveBase(::core::ArchiveP node)
    {_G();
	// Archive other instance variables here
	ASSERTF(NumBaseDimensions==8,BF("You changed NumBaseDimensions but didnt change Unit_O::archiveBase"));
	node->attribute("amt",this->_Amount);
	node->attribute("d0",this->_Powers[0]);
	node->attribute("w1",this->_Powers[1]);
	node->attribute("t2",this->_Powers[2]);
	node->attribute("c3",this->_Powers[3]);
	node->attribute("t4",this->_Powers[4]);
	node->attribute("l5",this->_Powers[5]);
	node->attribute("a6",this->_Powers[6]);
	node->attribute("r7",this->_Powers[7]);
    }
#endif

    void Unit_O::initialize()
    {_OF();
        this->Base::initialize();
	Dimension_O::zeroPowers(this->_Powers);
	this->_Amount = 1.0;
    }



    void Unit_O::incorporateUnit(Unit_sp u, double amountScale, int power)
    {_OF();
	/*Set the default system as units:*SI* */
//	core::Binder_sp unitDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
	// Repeat this block for multiple symbols
	{ _BLOCK_TRACEF(BF("Processing unit[%s]") % u->__repr__() );
	    u->adjustPowersAndAmountScale(power,this->_Powers,amountScale);
	    this->_Amount *= amountScale;
	}
    }


    void Unit_O::adjustPowersAndAmountScale(int power, int powers[], double& amountScale ) const
    {_OF();
	for ( int id=0; id<NumBaseDimensions; id++ )
	{
	    powers[id] += this->_Powers[id]*power;
	}
	double scale = ::pow(this->_Amount,power);
	amountScale *= scale;
    }

    string Unit_O::__repr__() const
    {
	stringstream ss;
	ss << this->_Amount << "*" << this->unitsOnlyAsString();
	return ss.str();
    }


    void Unit_O::test_set_amount(double amount)
    {_OF();
	this->_Amount = amount;
    }

    string Unit_O::unitsOnlyAsString() const
    {_OF();
	stringstream units;
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    if ( this->_Powers[i] != 0 )
	    {
		units << Dimension_O::baseDimensionUnitName(i);
		if ( this->_Powers[i] != 1)
		{
		    units << "^" << this->_Powers[i];
		}
		units << "*";
	    }
	}
	return units.str().substr(0,units.str().size()-1);
    };


    bool Unit_O::is_compatible(Unit_sp other,int power) const
    {_OF();
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    if ( this->_Powers[i] != other->_Powers[i]*power ) return false;
	}
	return true;
    }

    double Unit_O::conversion_factor_to(Unit_sp other, int power) const
    {_OF();
	if ( this->is_compatible(other,power) )
	{
	    double conversion = this->_Amount/::pow(other->_Amount,power);
	    return conversion;
	}
	SIMPLE_ERROR(BF("Units[%s] are not compatible with Unit[%s]^%d") 
		     % this->__repr__() % other->__repr__() % power );
    }

    Unit_sp Unit_O::copyWithoutName() const
    {_OF();
	Unit_sp result = Unit_O::create();
	Dimension_O::copyPowers(result->_Powers,this->_Powers);
        result->_Amount = this->_Amount;
	return result;
    }

    core::T_sp Unit_O::operator*(core::T_sp obj) const
    {_OF();
	core::T_sp result;
	if ( obj.isA<Unit_O>() )
	{
	    Unit_sp other = obj.as<Unit_O>();
	    Unit_sp uresult = this->copyWithoutName();
	    for ( int i=0; i<NumBaseDimensions; i++ )
	    {
		uresult->_Powers[i] += other->_Powers[i];
	    }
	    uresult->_Amount *= other->_Amount;
	    result = uresult;
	} else if ( obj.isA<core::Number_O>()
		    || obj.isA<geom::OVector3_O>()
		    || obj.isA<core::Array_O>() )
	{
	    core::T_sp val = obj->deepCopy();
	    result = Quantity_O::create(val,this->const_sharedThis<Unit_O>());
	} else if ( obj.isA<Quantity_O>() )
	{
	    SIMPLE_ERROR(BF("Handle Unit*Quantity"));
	} else
	{
	    SIMPLE_ERROR(BF("Handle Unit*XXX where XXX=%s") % obj->__class()->classNameAsString());
	}
	return result;
    }


    core::T_sp Unit_O::operator/(core::T_sp obj) const
    {_OF();
	if ( obj.isA<Unit_O>() )
	{
	    Unit_sp other = obj.as<Unit_O>();
	    Unit_sp result = this->copyWithoutName();
	    for ( int i=0; i<NumBaseDimensions; i++ )
	    {
		result->_Powers[i] -= other->_Powers[i];
	    }
	    result->_Amount /= other->_Amount;
	    return result;
	} else
	{
	    SIMPLE_ERROR(BF("Handle Unit/XXX where XXX=%s") % obj->__class()->classNameAsString());
	}
    }


}; /* units */
