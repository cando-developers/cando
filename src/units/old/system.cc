/*
    File: system.cc
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
 
This is an open source license for the CANDO software from Temple University, but it is not the only one. Contact Temple University at mailto:techtransfer@temple.edu if you would like a different license.
*/
/* -^- */
#define	DEBUG_LEVEL_NONE

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






    Object_sp System_O::__init__(::core::Executable_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {
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
    {
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
