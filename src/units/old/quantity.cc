#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "quantity.h"
#include "system.h"
#include "float.h"
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//


    Object_sp prim_quantity(core::Executable_sp exec, core::Cons_sp args,
			    core::Environment_sp env, core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	LOG(BF("Initializing Quantity_O with args[%s]") % args->__repr__() );
	Object_sp oamount = args->ocar();
	if ( oamount->isAssignableTo<Float_O>() )
	{
	    Float_sp quant = _lisp->create<Float_O>();
	    quant->__init__(exec,args,env,lisp);
	    return quant;
	}
	THROW(_lisp->error(BF("Illegal quantity")));
    }




    EXPOSE_CLASS(Quantity_O);

    void Quantity_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<Quantity_O>()
	    .def("dimensions-match",&Quantity_O::dimensionsMatch)
	    .def_raw("amount",&Quantity_O::amount_matchUnits)
	    .def("+",&Quantity_O::add)
	    .def("-",&Quantity_O::sub)
	    .def("*",&Quantity_O::mul)
	    .def("/",&Quantity_O::div)
	    .def("add",&Quantity_O::add)
	    .def("sub",&Quantity_O::sub)
	    .def("mul",&Quantity_O::mul)
	    .def("div",&Quantity_O::div)
	    .def("reciprocal",&Quantity_O::reciprocal)
	    .def("negate",&Quantity_O::negate)
	    .def("sqrt",&Quantity_O::sqrt)
	    .def("mul-by-scalar",&Quantity_O::mulByScalar)
	    .def("div-by-scalar",&Quantity_O::divByScalar)
	    .def("pow",&Quantity_O::pow)
	    .def("unsafe-amount",&Quantity_O::unsafeAmount)
	    .def("isnan",&Quantity_O::isnan)
	    .def("dump",&Quantity_O::dump)
	;
    }


    void Quantity_O::exposePython(core::Lisp_sp lisp)
    {_G();
	PYTHON_CLASS(UnitsPkg,Quantity,"","",_lisp)
	    .def("dimensions-match",&Quantity_O::dimensionsMatch)
	    .def("amount",boost::python::raw_method<Quantity_O>("amount",&Quantity_O::amount_matchUnits,_lisp))
	    .def("+",&Quantity_O::add)
	    .def("-",&Quantity_O::sub)
	    .def("*",&Quantity_O::mul)
	    .def("/",&Quantity_O::div)
	    .def("add",&Quantity_O::add)
	    .def("sub",&Quantity_O::sub)
	    .def("mul",&Quantity_O::mul)
	    .def("div",&Quantity_O::div)
	    .def("reciprocal",&Quantity_O::reciprocal)
	    .def("negate",&Quantity_O::negate)
	    .def("sqrt",&Quantity_O::sqrt)
	    .def("mul-by-scalar",&Quantity_O::mulByScalar)
	    .def("div-by-scalar",&Quantity_O::divByScalar)
	    .def("pow",&Quantity_O::pow)
	    .def("unsafeAmount",&Quantity_O::unsafeAmount)
	    .def("isnan",&Quantity_O::isnan)
	    .def("dump",&Quantity_O::dump)
	;
    }



    void Quantity_O::parseDimensionSymbols(int powers[NumBaseDimensions], double& amountScale,
							  core::Symbol_sp unit, int power)
    {_OF();
	/*Set the default system as units:*SI* */
	core::Symbol_sp system = _lisp->symbol(_sym_UnitsPkg_StarSISystemStar);
	core::Binder_sp dimensionDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
	// Repeat this block for multiple symbols
	{ _BLOCK_TRACEF(BF("Processing unit[%s]") % unit->__repr__() );
	    core::Binder_O::iterator it = dimensionDatabase->find(unit);
	    ASSERTF(it!=dimensionDatabase->end(),BF("Could not find BaseUnit[%s]") % unit->__repr__() );
	    Dimension_sp dimension = (*it).second.as<Dimension_O>();
	    dimension->adjustPowersAndAmountScale(power,powers,amountScale);
	}
    }

    void Quantity_O::parseDimensionSymbols(int powers[NumBaseDimensions], double& amountScale,
					       core::Symbol_sp unit1, int power1,
					       core::Symbol_sp unit2, int power2)
    {_OF();
	this->parseDimensionSymbols(powers,amountScale,unit1,power1);
	this->parseDimensionSymbols(powers,amountScale,unit2,power2);
	// repeat if there are more
    }


    core::Symbol_sp Quantity_O::parseDimensionSymbolList(int powers[NumBaseDimensions], double& amountScale,
						       core::Cons_sp list, bool allowFinalSystemArg) const
    {_OF();
	/*Set the default system as units:*SI* */
	core::Symbol_sp system = _lisp->symbol(_sym_UnitsPkg_StarSISystemStar);
	core::Binder_sp dimensionDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
	core::Cons_sp dimCur = list;
	while (dimCur->notNil())
	{
	    Object_sp odim = dimCur->ocar();
	    if ( odim->isAssignableTo<core::Symbol_O>() )
	    {
		core::Symbol_sp sdim = odim.as<core::Symbol_O>();
		if ( sdim == _lisp->symbol(_sym_kw_system) )
		{
		    if ( dimCur->ocadr()->isAssignableTo<core::Symbol_O>() )
		    {
			system = dimCur->ocadr().as<core::Symbol_O>();
			dimCur = dimCur->cdr();
		    }
		} else
		{
		    core::Binder_O::iterator it = dimensionDatabase->find(sdim);
		    if ( it!=dimensionDatabase->end() )
		    {
			Dimension_sp dimension = (*it).second.as<Dimension_O>();
			int power = 1;
			if ( dimCur->ocadr()->isAssignableTo<core::Rational_O>() )
			{
			    power = dimCur->ocadr().as<core::Rational_O>()->asInt();
			    dimCur = dimCur->cdr();
			}
			dimension->adjustPowersAndAmountScale(power,powers,amountScale);
		    } else
		    {
			THROW(_lisp->error(BF("Unknown dimension[%s]") % sdim->__repr__()));
		    }
		}
	    } else
	    {
		THROW(_lisp->error(BF("Units definition expected symbol or :system but got[%s]") % dimCur->ocar()->__repr__() ));
	    }
	    dimCur = dimCur->cdr();
	}
	return system;
    }

#if 0
    /*! If last argument is a System then use that system, otherwise use the default units:*SI* */
    Object_sp Quantity_O::__init__(::core::Executable_sp exec, ::core::Cons_sp allargs, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	core::Cons_sp args = allargs->cdr();
	LOG(BF("Initializing Quantity_O with args[%s]") % args->__repr__() );
	double amount = 1.0;
	Object_sp oamount = args->ocar();
	if ( oamount->notNil() )
	{
	    amount = oamount.as<core::Number_O>()->asReal();
	}
	/*Set the default system as units:*SI* */
	this->_System = this->parseDimensionSymbolList(this->_Powers,amount,args->cdr().as<core::Cons_O>(),true);
	this->_Amount = amount;
	return _lisp->onil();
    }
#endif

    void Quantity_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
        node->archiveObjectIfDefined("system",this->_System);
	node->attribute("d0",this->_Powers[0]);
	node->attribute("w1",this->_Powers[1]);
	node->attribute("t2",this->_Powers[2]);
	node->attribute("c3",this->_Powers[3]);
	node->attribute("t4",this->_Powers[4]);
	node->attribute("l5",this->_Powers[5]);
	node->attribute("a6",this->_Powers[6]);
    }


    void Quantity_O::initialize()
    {_OF();
        this->Base::initialize();
	this->_Powers[distanceMeter] = 0;
	this->_Powers[weightKilogram] = 0;
	this->_Powers[timeSecond] = 0;
	this->_Powers[currentAmpere] = 0;
	this->_Powers[temperatureKelvin] = 0;
	this->_Powers[luminosityCandela] = 0;
	this->_Powers[amountMole] = 0;
	this->_System = _lisp->symbol(_sym_UnitsPkg_StarSISystemStar);
    }


    void Quantity_O::throwOnPowersMisMatch(int matchPowers[]) const
    {_OF();
	if ( !this->powersMatch(matchPowers) )
	{
	    Quantity_sp temp = RP_Create<Quantity_O>(_lisp);
	    for (int ii=0; ii<NumBaseDimensions; ii++) temp->_Powers[ii] = matchPowers[ii];
	    temp->_System = this->_System;
	    THROW(_lisp->error(BF("You tried to get the amount of the %s[%s] in terms of units[%s]"
				  "- cannot get the amount of a Quantity unless you match the units")
			       % this->className()
			       % this->__repr__()
			       % temp->__repr__() ));
	}
    }


    bool Quantity_O::powersMatch(const int otherPowers[NumBaseDimensions]) const
    {_OF();
	return ((this->_Powers[distanceMeter] == otherPowers[distanceMeter]) &&
		(this->_Powers[weightKilogram] == otherPowers[weightKilogram]) &&
		(this->_Powers[timeSecond] == otherPowers[timeSecond]) &&
		(this->_Powers[currentAmpere] == otherPowers[currentAmpere]) &&
		(this->_Powers[temperatureKelvin] == otherPowers[temperatureKelvin]) &&
		(this->_Powers[luminosityCandela] == otherPowers[luminosityCandela]) &&
		(this->_Powers[amountMole] == otherPowers[amountMole]));
    }



    double Quantity_O::setDimension(BaseUnit_sp dim, int power)
    {_OF();
	IMPLEMENT_ME();
#if 0
	this->_Powers[dim->_BaseUnit] = power;
	double amountMultiplier = ::pow(dim->_CanonicalToUnitScale,-power);
	return amountMultiplier;
#endif
    }


    Object_sp Quantity_O::amount() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    Object_sp Quantity_O::amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	IMPLEMENT_ME();
    }


    Object_sp Quantity_O::unsafeAmount() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    string Quantity_O::__repr__() const
    {_OF();
	stringstream ss;
	double amount = 1.0;
	ss << this->_class()->classNameSymbol()->__repr__()<< "*" <<this->system()->unitsAsString(amount,this->sharedThis<Quantity_O>());
	return ss.str();
    }


    bool Quantity_O::dimensionsMatch(Quantity_sp other) const
    {_OF();
	if ( this->powersMatch(other->_Powers) ) return true;
	return false;
    }







    Quantity_sp Quantity_O::add(Quantity_sp other) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }

	
    Quantity_sp Quantity_O::sub(Quantity_sp other) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }

    Quantity_sp Quantity_O::mul(Quantity_sp other) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }

    Quantity_sp Quantity_O::div(Quantity_sp other) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }



    Quantity_sp Quantity_O::reciprocal() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    Quantity_sp Quantity_O::negate() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    Quantity_sp Quantity_O::sqrt() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }





    Quantity_sp Quantity_O::mulByScalar(double scalar) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }



    Quantity_sp Quantity_O::divByScalar(double scalar) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }



    Quantity_sp Quantity_O::pow(int power) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    bool Quantity_O::isnan() const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


    System_sp Quantity_O::system() const
    {_OF();
	System_sp system = _lisp->symbol(_sym_UnitsPkg_StarSISystemStar)->symbolValue().as<System_O>();
	ASSERTF(this->_System->notNil(),BF("The system symbol should never be nil"));
	if (!this->_System->hasDynamicValue())
	{
	    system = this->_System->symbolValue().as<System_O>();
	}
	return system;
    }




}; /* units */
