#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "system.h"
#include "quantity.h"
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(System_O);

    void System_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<System_O>()
//	.def_raw("__init__",&System_O::__init__,"(self)")
	;
    }


    void System_O::exposePython(core::Lisp_sp lisp)
    {_G();
	PYTHON_CLASS(UnitsPkg,System,"","",_lisp)
//	.def_raw("__init__",&initArgs("(self)")
	;
    }


    Object_sp System_O::__init__(::core::Executable_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _lisp->onil();
    }

    void System_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }


    void System_O::initialize()
    {_OF();
        this->Base::initialize();
    }



    System_sp System_O::create(string const& systemName,
				uint distanceSymbolIndex,
				uint weightSymbolIndex,
				uint timeSymbolIndex,
				uint currentSymbolIndex,
				uint temperatureSymbolIndex,
				uint luminositySymbolIndex,
				uint amountSymbolIndex,
				core::Lisp_sp lisp)
    {_G();
	System_sp system = RP_Create<System_O>(lisp);
	core::Binder_sp dimensionDatabase = _lisp->symbol(_sym_UnitsPkg_StarUnitDatabaseStar)->symbolValue().as<core::Binder_O>();
	system->_Name = systemName;
	system->_BaseUnits[distanceMeter] = dimensionDatabase->lookup(_lisp->symbol(distanceSymbolIndex)).as<BaseUnit_O>();
	system->_BaseUnits[weightKilogram] = dimensionDatabase->lookup(_lisp->symbol(weightSymbolIndex)).as<BaseUnit_O>();
	system->_BaseUnits[timeSecond] = dimensionDatabase->lookup(_lisp->symbol(timeSymbolIndex)).as<BaseUnit_O>();
	system->_BaseUnits[currentAmpere] = dimensionDatabase->lookup(_lisp->symbol(currentSymbolIndex)).as<BaseUnit_O>();
	system->_BaseUnits[temperatureKelvin] = dimensionDatabase->lookup(_lisp->symbol(temperatureSymbolIndex)).as<BaseUnit_O>();
	system->_BaseUnits[luminosityCandela] = dimensionDatabase->lookup(_lisp->symbol(luminositySymbolIndex)).as<BaseUnit_O>();
	system->_BaseUnits[amountMole] = dimensionDatabase->lookup(_lisp->symbol(amountSymbolIndex)).as<BaseUnit_O>();
	return system;
    }


    string System_O::unitsAsString(double& amount, Quantity_sp quantity) const
    {_OF();
	if ( this.nilp() ) return "-nil-quantity-";
	stringstream units;
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    ASSERT_NOT_NULL(this->_BaseUnits[i]);
	    amount /= this->_BaseUnits[i]->canonicalToUnitScale();
	    if ( quantity->_Powers[i] != 0 )
	    {
		units << this->_BaseUnits[i]->label();
		if ( quantity->_Powers[i] != 1)
		{
		    units << "^" << quantity->_Powers[i];
		}
		units << "*";
	    }
	}
	return units.str().substr(0,units.str().size()-1);
    };



    string System_O::quantityAsString(double amount, Quantity_sp quantity) const
    {_OF();
	if ( this.nilp() ) return "-nil-quantity-";
	stringstream units;
	units << this->unitsAsString(amount,quantity);
	stringstream result;
	result << amount << "*"<< units.str();
	return result.str();
    };

    

}; /* units */
