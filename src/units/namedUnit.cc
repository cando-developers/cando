#define	DEBUG_LEVEL_FULL

#include <clasp/core/common.h>
#include <clasp/core/environment.h>
//#include "core/serialize.h"
#include <cando/units/namedUnit.h>
#include <clasp/core/wrappers.h>
namespace units
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(units,NamedUnit_O);

    NamedUnit_sp NamedUnit_O::create(core::Symbol_sp symbolId, double scale, Unit_sp unit, string const& label, string const& code)
    {_G();
	NamedUnit_sp n = _lisp->create<NamedUnit_O>();
	n->_Symbol = symbolId;
	n->_Label = label;
	n->_Code = code;
	Dimension_O::copyPowers(n->_Powers,unit->_Powers);
	n->_Amount = scale*unit->_Amount;
	return n;
    }



    NamedUnit_sp NamedUnit_O::create(core::Symbol_sp symbolId, double scale, int baseDimension, int power, string const& label, string const& code)
    {_G();
	NamedUnit_sp n = _lisp->create<NamedUnit_O>();
	n->_Symbol = symbolId;
	n->_Label = label;
	n->_Code = code;
	n->_Powers[baseDimension] = power;
	n->_Amount = scale;
	return n;
    }



    NamedUnit_sp NamedUnit_O::create(core::Symbol_sp symbolId, string const& label, string const& code)
    {_G();
	NamedUnit_sp n = _lisp->create<NamedUnit_O>();
	n->_Symbol = symbolId;
	n->_Label = label;
	n->_Code = code;
	n->_Amount = 1.0;
	return n;
    }



    void NamedUnit_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<NamedUnit_O>()
//	.initArgs("(self)")
	;
    }


    void NamedUnit_O::exposePython(::core::Lisp_sp lisp)
    {
#ifdef USEBOOSTPYTHON
	PYTHON_CLASS(UnitsPkg,NamedUnit,"","",lisp);
	;
#endif
    }

#if 0
    core::T_sp NamedUnit_O::__init__(::core::Function_sp exec, ::core::Cons_sp args, ::core::Environment_sp env)
    {_G();
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
    {_OF();
        this->Base::initialize();
    }

    

}; /* units */
