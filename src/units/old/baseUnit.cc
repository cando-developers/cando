/*
    File: baseUnit.cc
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
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//





    Object_sp BaseUnit_O::__init__(::core::Executable_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _lisp->onil();
    }

    void BaseUnit_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }


    void BaseUnit_O::initialize()
    {_OF();
        this->Base::initialize();
    }

    BaseUnit_sp BaseUnit_O::create(uint symbolIndex, double scale, string const& label,
					     BaseDimensionEnum baseDim, core::Lisp_sp lisp)
    {
	BaseUnit_sp bd = RP_Create<BaseUnit_O>(lisp);
	bd->_BaseUnit = baseDim;
	bd->_CanonicalToUnitScale = 1.0/scale;
//	bd->_UnitToCanonicalScale = 1.0/scale;
	bd->_Label = label;
	return bd;
    }



    void BaseUnit_O::adjustPowersAndAmountScale(int power, int powers[], double& amountScale ) const
    {_OF();
	ASSERTF(this->_BaseUnit<NumBaseDimensions,BF("Bad base dimension index[%d]") % this->_BaseUnit);
	powers[this->_BaseUnit] = power;
	double scale = ::pow(this->_CanonicalToUnitScale,power);
	amountScale /= scale;
    }

}; /* units */
