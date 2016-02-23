/*
    File: float.cc
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
#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "float.h"
#include "system.h"
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//





    Float_sp Float_O::create(double amount, core::Symbol_sp dimension, int power, core::Symbol_sp system,core::Lisp_sp lisp)
    {
	Float_sp f = _lisp->create<Float_O>();
	f->_System = system;
	f->parseDimensionSymbols(f->_Powers,amount,dimension,power);
	f->_Amount = amount;
	return f;
    }



    /*! If last argument is a System then use that system, otherwise use the default units:*SI* */
    Object_sp Float_O::__init__(::core::Executable_sp exec, ::core::Cons_sp allargs, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	core::Cons_sp args = allargs->cdr();
	LOG(BF("Initializing Float_O with args[%s]") % args->__repr__() );
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

    void Float_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->attribute("amt",this->_Amount);
    }


    void Float_O::initialize()
    {_OF();
        this->Base::initialize();
	this->_Amount = 0.0;
    }


    Object_sp Float_O::amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	this->parseDimensionSymbolList(matchPowers,amountScale,args->cdr(),false);
	this->throwOnPowersMisMatch(matchPowers);
	return core::Double_float_O::create(this->_Amount/amountScale,_lisp);
    }

    double Float_O::amount_matchUnit(core::Symbol_sp dim, int power)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	this->parseDimensionSymbols(matchPowers,amountScale,dim,power);
	this->throwOnPowersMisMatch(matchPowers);
	return this->_Amount/amountScale;
    }




    Object_sp Float_O::amount() const
    {_OF();
	int zeroPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	if ( this->powersMatch(zeroPowers) )
	{
	    return core::Double_float_O::create(this->_Amount,_lisp);
	}
	THROW(_lisp->error(BF("You cannot get the amount of a Quantity unless it is unit-less!")));
    }







    Object_sp Float_O::unsafeAmount() const
    {_OF();
	return core::Double_float_O::create(this->_Amount,_lisp);
    }


    string Float_O::__repr__() const
    {_OF();
	stringstream ss;
	return this->system()->quantityAsString(this->_Amount,this->sharedThis<Float_O>());
    }







    Quantity_sp Float_O::add(Quantity_sp other) const
    {_OF();
	Float_sp result = RP_Copy(this);
	if ( other->notNil() )
	{
	    if ( !this->dimensionsMatch(other) )
	    {
		THROW(_lisp->error(BF("Mismatch in dimensions in add - you tried to add[%s] to [%s]") % this->__repr__() % other->__repr__()));
	    }
	    result->_Amount += other.as<Float_O>()->_Amount;
	}
	return result;
    }

	
    Quantity_sp Float_O::sub(Quantity_sp other) const
    {_OF();
	Float_sp result = RP_Copy(this);
	if ( other->notNil() )
	{
	    if ( !this->dimensionsMatch(other) )
	    {
		THROW(_lisp->error(BF("Mismatch in dimensions in sub - you tried to sub[%s] to [%s]") % this->__repr__() % other->__repr__()));
	    }
	    result->_Amount -= other.as<Float_O>()->_Amount;
	}
	return result;
    }

    Quantity_sp Float_O::negate() const
    {_OF();
	Float_sp result = RP_Copy(this);
	result->_Amount = - result->_Amount;
	return result;
    }


    Quantity_sp Float_O::mul(Quantity_sp other) const
    {_OF();
	Float_sp result = RP_Copy(this);
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    result->_Powers[i] += other->_Powers[i];
	}
	result->_Amount *= other.as<Float_O>()->_Amount;
	return result;
    }

    Quantity_sp Float_O::div(Quantity_sp other) const
    {_OF();
	Float_sp result = RP_Copy(this);
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    result->_Powers[i] -= other->_Powers[i];
	}
	result->_Amount /= other.as<Float_O>()->_Amount;
	return result;
    }



    Quantity_sp Float_O::reciprocal() const
    {_OF();
	Float_sp result = RP_Copy(this);
	if ( this->_Amount == 0.0 )
	{
	    THROW(core::MathException_O::create("reciprocal of 0.0",_lisp));
	}
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    result->_Powers[i] = -result->_Powers[i];
	}
	result->_Amount = 1.0 / result->_Amount;
	return result;
    }

    Quantity_sp Float_O::sqrt() const
    {_OF();
	Float_sp result = RP_Copy(this);
	if ( this->_Amount < 0.0  )
	{
	    THROW(core::MathException_O::create("sqrt of negative number",_lisp));
	}
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    if ( (result->_Powers[i] % 2) == 1 )
	    {
		THROW(core::MathException_O::create("sqrt will result in fractional unit",_lisp));
	    }
	    result->_Powers[i] /= 2;
	}
	result->_Amount = ::sqrt(result->_Amount);
	return result;
    }





    Quantity_sp Float_O::mulByScalar(double scalar) const
    {_OF();
	Float_sp result = RP_Copy(this);
	result->_Amount *= scalar;
	return result;
    }



    Quantity_sp Float_O::divByScalar(double scalar) const
    {_OF();
	Float_sp result = RP_Copy(this);
	if ( scalar == 0.0 )
	{
	    THROW(_lisp->create<core::MathException_O>("Divide by zero in units::Float_O::divByScalar"));
	}
	result->_Amount /= scalar;
	return result;
    }



    Quantity_sp Float_O::pow(int power) const
    {_OF();
	Float_sp result = RP_Copy(this);
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    result->_Powers[i] *= power;
	}
	result->_Amount = ::pow(result->_Amount,power);
	return result;
    }


    bool Float_O::isnan() const
    {_OF();
	return this->_Amount != this->_Amount;
    }


    string Float_O::dump() const
    {_OF();
	return this->__repr__();
    }



}; /* units */
