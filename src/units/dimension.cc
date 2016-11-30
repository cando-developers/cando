/*
    File: dimension.cc
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

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
#include <cando/units/dimension.h>
#include <clasp/core/wrappers.h>
namespace units
{

// ----------------------------------------------------------------------
//





#if 0
    core::T_sp Dimension_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _lisp->onil();
    }
#endif


#ifdef XML_ARCHIVE
    void Dimension_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	// Archive other instance variables here
    }
#endif

    void Dimension_O::initialize()
    {_OF();
        this->Base::initialize();
    }



    void Dimension_O::zeroPowers(int powers[])
    {
	for (int i=0; i<NumBaseDimensions; i++ ) powers[i] = 0;
    }

    void Dimension_O::copyPowers(int dest[], int const src[])
    {
	for (int i=0; i<NumBaseDimensions; i++ ) dest[i] = src[i];
    }

    void Dimension_O::adjustPowersAndAmountScale(int power, int powers[], double& amountScale ) const
    {_OF();
	SUBCLASS_MUST_IMPLEMENT();
    }


string Dimension_O::baseDimensionUnitName(int i)
    {
	ASSERTF(NumBaseDimensions==8,BF("When you change NumBaseDimensions you must update Dimension_O::baseDimensionUnitName"));
	ASSERTF(i>=0 && i < NumBaseDimensions,BF("Illegal base dimension index[%d]") % i);
	switch (i)
	{
	case 0: return "m";
	case 1: return "kg";
	case 2: return "s";
	case 3: return "A";
	case 4:	return "K";
	case 5: return "C";
	case 6: return "mole";
	case 7: return "rad";
	};
	SIMPLE_ERROR(BF("Illegal dimension index[%d]") % i);
    }
    

}; /* units */
