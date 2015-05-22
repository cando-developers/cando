#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(BaseUnit_O);

    void BaseUnit_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<BaseUnit_O>()
//	.def_raw("__init__",&BaseUnit_O::__init__,"(self)")
	;
    }

    void BaseUnit_O::exposePython(core::Lisp_sp lisp)
    {_G();
	PYTHON_CLASS(UnitsPkg,BaseUnit,"","",_lisp)
	;
    }


    Object_sp BaseUnit_O::__init__(::core::Executable_sp exec, ::core::Cons_sp args, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {_G();
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
    {_G();
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
