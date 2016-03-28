/*
    File: xyz.cc
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
#include "vector3.h"
#include "core/ovector3.h"
#include "xyz.h"
#include "float.h"
#include "system.h"
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//






    Xyz_sp Xyz_O::create(Vector3 pos, Quantity_sp temp, core::Lisp_sp lisp)
    {	
	Xyz_sp x = _lisp->create<Xyz_O>();
	for ( int i=0; i<NumBaseDimensions; i++ ) x->_Powers[i] = temp->_Powers[i];
	x->_Pos = pos;
	x->_System = temp->_System;
	return x;
    }


    /*! If last argument is a System then use that system, otherwise use the default units:*SI* */
    Object_sp Xyz_O::__init__(::core::Executable_sp exec, ::core::Cons_sp allargs, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	core::Cons_sp args = allargs->cdr();
	LOG(BF("Initializing Xyz_O with args[%s]") % _rep_(args) );
	double amount = 1.0;
	core::OVector3_sp vamount = args->ocar().as<core::OVector3_O>();
	/*Set the default system as units:*SI* */
	this->_System = this->parseDimensionSymbolList(this->_Powers,amount,args->cdr().as<core::Cons_O>(),true);
	this->_Pos = vamount->timesScalar(amount);
	return _lisp->onil();
    }

    void Xyz_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->archivePlainObjectIfDefined<Vector3>( "pos","Vector3",
						    this->_Pos.isDefined(), this->_Pos );
    }


    void Xyz_O::initialize()
    {_OF();
        this->Base::initialize();
	this->_Pos.set(0.0,0.0,0.0);
    }


    Object_sp Xyz_O::amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	this->parseDimensionSymbolList(matchPowers,amountScale,args->cdr(),false);
	if ( !this->powersMatch(matchPowers) )
	{
	    Xyz_sp temp = RP_Create<Xyz_O>(_lisp);
	    for (int ii=0; ii<NumBaseDimensions; ii++) temp->_Powers[ii] = matchPowers[ii];
	    temp->_System = this->_System;
	    THROW(_lisp->error(BF("You tried to get the amount of the quantity[%s] in terms of units[%s]"
				  "- cannot get the amount of a Quantity unless you match the units")
			       % this->__repr__()
			       % _rep_(temp) ));
	}
	return core::OVector3_O::create(this->_Pos.multiplyByScalar(1.0/amountScale),_lisp);
    }




    Object_sp Xyz_O::amount() const
    {_OF();
	int zeroPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	if ( this->powersMatch(zeroPowers) )
	{
	    return core::OVector3_O::create(this->_Pos,_lisp);
	}
	THROW(_lisp->error(BF("You cannot get the amount of a Quantity unless it is unit-less!")));
    }







    Object_sp Xyz_O::unsafeAmount() const
    {_OF();
	return core::OVector3_O::create(this->_Pos,_lisp);
    }


    string Xyz_O::__repr__() const
    {_OF();
	double amount = 1.0;
	Vector3 pos = this->_Pos;
	stringstream units;
	units << "*" << this->system()->unitsAsString(amount,this->sharedThis<Xyz_O>());
	pos = pos.multiplyByScalar(amount);
	stringstream result;
	result << pos.asString() << units.str();
	return result.str();
    }







    Quantity_sp Xyz_O::add(Quantity_sp other) const
    {_OF();
	Xyz_sp result = RP_Copy(this);
	if ( other->notNil() )
	{
	    if ( !this->dimensionsMatch(other) )
	    {
              THROW(_lisp->error(BF("Mismatch in dimensions in add - you tried to add[%s] to [%s]") % this->__repr__() % _rep_(other)));
	    }
	    result->_Pos = result->_Pos +  (other.as<Xyz_O>()->_Pos);
	}
	return result;
    }

	
    Quantity_sp Xyz_O::sub(Quantity_sp other) const
    {_OF();
	Xyz_sp result = RP_Copy(this);
	if ( other->notNil() )
	{
	    if ( !this->dimensionsMatch(other) )
	    {
              THROW(_lisp->error(BF("Mismatch in dimensions in sub - you tried to sub[%s] to [%s]") % this->__repr__() % _rep_(other)));
	    }
	    result->_Pos = result->_Pos - other.as<Xyz_O>()->_Pos;
	}
	return result;
    }





    Quantity_sp Xyz_O::mulByFloat(Float_sp f) const
    {_OF();
	Xyz_sp result = RP_Copy(this);
	for ( int i=0; i<NumBaseDimensions; i++ )
	{
	    result->_Powers[i] += f->_Powers[i];
	}
	result->_Pos = result->_Pos.multiplyByScalar(f->_Amount);
	return result;
    }





    bool Xyz_O::isnan() const
    {_OF();
	return (
	    (this->_Pos.getX() != this->_Pos.getX())||
	    (this->_Pos.getY() != this->_Pos.getY())||
	    (this->_Pos.getZ() != this->_Pos.getZ()));
    }





}; /* units */
