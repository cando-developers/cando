/*
    File: namedUnit.cc
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
#include <cando/units/namedUnit.h>
#include <clasp/core/wrappers.h>
namespace units
{

// ----------------------------------------------------------------------
//


    NamedUnit_sp NamedUnit_O::create(core::Symbol_sp symbolId, double scale, Unit_sp unit, string const& label, string const& code)
    {
	NamedUnit_sp n = _lisp->create<NamedUnit_O>();
	n->_Symbol = symbolId;
	n->_Label = label;
	n->_Code = code;
	Dimension_O::copyPowers(n->_Powers,unit->_Powers);
	n->_Amount = scale*unit->_Amount;
	return n;
    }



    NamedUnit_sp NamedUnit_O::create(core::Symbol_sp symbolId, double scale, int baseDimension, int power, string const& label, string const& code)
    {
	NamedUnit_sp n = _lisp->create<NamedUnit_O>();
	n->_Symbol = symbolId;
	n->_Label = label;
	n->_Code = code;
	n->_Powers[baseDimension] = power;
	n->_Amount = scale;
	return n;
    }



    NamedUnit_sp NamedUnit_O::create(core::Symbol_sp symbolId, string const& label, string const& code)
    {
	NamedUnit_sp n = _lisp->create<NamedUnit_O>();
	n->_Symbol = symbolId;
	n->_Label = label;
	n->_Code = code;
	n->_Amount = 1.0;
	return n;
    }






#if 0
    core::T_sp NamedUnit_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env)
    {
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	return _lisp->onil();
    }
#endif



#ifdef XML_ARCHIVE
    void NamedUnit_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->attributeSymbolIfNotNil("symbol",this->_Symbol);
	node->attribute("label",this->_Label);
	node->attribute("code",this->_Code);
	// Archive other instance variables here
    }
#endif

#ifdef OLD_SERIALIZE
    void NamedUnit_O::serialize(serialize::SNode node)
    {
        this->Base::serialize(node);
	node->attributeSymbolIfNotNil("symbol",this->_Symbol);
	node->attribute("label",this->_Label);
	node->attribute("code",this->_Code);
    }
#endif

    void NamedUnit_O::initialize()
    {
        this->Base::initialize();
    }

    

}; /* units */
