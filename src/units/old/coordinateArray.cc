#define	DEBUG_LEVEL_FULL

#include "core/common.h"
#include "core/environment.h"
#include "coordinateArray.h"
#include "core/ovector3.h"
#include "system.h"
#include "xyz.h"
#include "singleDimension.h"
#include "core/wrappers.h"
namespace units
{

// ----------------------------------------------------------------------
//

    EXPOSE_CLASS(CoordinateArray_O);

    void CoordinateArray_O::exposeCando(::core::Lisp_sp lisp)
    {
	::core::class_<CoordinateArray_O>()
	    .def_raw("__init__",&CoordinateArray_O::__init__,"(self &rest args)")
	    .def_raw("getElement",&CoordinateArray_O::getElement_matchUnits)
	    .def("coordinate-array-append",&CoordinateArray_O::append)
	    .def("copy-CoordinateArray",&CoordinateArray_O::copy)
	;
    }


    void CoordinateArray_O::exposePython(core::Lisp_sp lisp)
    {_G();
	PYTHON_CLASS(UnitsPkg,CoordinateArray,"","",_lisp)
	    .def("coordinate-array-append",&CoordinateArray_O::append)
	    .def("copy-CoordinateArray",&CoordinateArray_O::copy)
	;
    }



    CoordinateArray_sp CoordinateArray_O::create(core::CoordinateArray_sp coords, core::Symbol_sp dimension, int power, core::Symbol_sp system,core::Lisp_sp lisp)
    {_G();
	CoordinateArray_sp uc = _lisp->create<CoordinateArray_O>();
	/*Set the default system as units:*SI* */
	core::CoordinateArray_sp camount = coords->copy();
	uc->_System = system;
	double amountScale = 1.0;
	uc->parseDimensionSymbols(uc->_Powers,amountScale,dimension,power);
	camount->multiplyBy(amountScale);
	uc->_Coordinates = camount;
	return uc;
    }


    /*! If last argument is a System then use that system, otherwise use the default units:*SI* */
    Object_sp CoordinateArray_O::__init__(::core::Executable_sp exec, ::core::Cons_sp allargs, ::core::Environment_sp env, ::core::Lisp_sp lisp)
    {_G();
//      this->Base::__init__(exec,args,env,lisp);
//      arg = from_object<XXXX>::convert(env->lookup(this->Package(),"YYY"));
	core::Cons_sp args = allargs->cdr();
	LOG(BF("Initializing CoordinateArray_O with args[%s]") % args->__repr__() );
	double amount = 1.0;
	core::CoordinateArray_sp camount = args->ocar().as<core::CoordinateArray_O>();
	if ( camount.nilp() )
	{
	    camount = _lisp->create<core::CoordinateArray_O>();
	}
	/*Set the default system as units:*SI* */
	this->_System = this->parseDimensionSymbolList(this->_Powers,amount,args->cdr().as<core::Cons_O>(),true);
	camount->multiplyBy(amount);
	this->_Coordinates = camount;
	return _lisp->onil();
    }

    void CoordinateArray_O::archiveBase(::core::ArchiveP node)
    {
        this->Base::archiveBase(node);
	node->archiveObject("coords",this->_Coordinates);
    }


    void CoordinateArray_O::initialize()
    {_OF();
        this->Base::initialize();
	this->_Coordinates = core::CoordinateArray_O::create(_lisp);
    }


    Object_sp CoordinateArray_O::amount_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	this->parseDimensionSymbolList(matchPowers,amountScale,args->cdr(),false);
	this->throwOnPowersMisMatch(matchPowers);
	core::CoordinateArray_sp ccopy = this->_Coordinates->copy();
	ccopy->multiplyBy(1.0/amountScale);
	return ccopy;
    }



    Object_sp CoordinateArray_O::getElement_matchUnits(core::Executable_sp exec, core::Cons_sp args, core::Environment_sp env, core::Lisp_sp lisp)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	int i = args->ocadr().as<core::Rational_O>()->asInt();
	this->parseDimensionSymbolList(matchPowers,amountScale,args->cddr().as<core::Cons_O>(),false);
	this->throwOnPowersMisMatch(matchPowers);
	Vector3 vc = this->_Coordinates->getElement(i).multiplyByScalar(1.0/amountScale);
	return core::OVector3_O::create(vc,_lisp);
    }


    Xyz_sp CoordinateArray_O::getElement(int i) const
    {_OF();
	return Xyz_O::create(this->_Coordinates->getElement(i),this->sharedThis<CoordinateArray_O>(),_lisp);
    }



    Vector3 CoordinateArray_O::getElement(int i, core::Symbol_sp unit1, int power1)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	this->parseDimensionSymbols(matchPowers,amountScale,unit1,power1);
	this->throwOnPowersMisMatch(matchPowers);
	Vector3 vc = this->_Coordinates->getElement(i).multiplyByScalar(1.0/amountScale);
	return vc;
    }

    Vector3 CoordinateArray_O::getElement(int i, core::Symbol_sp unit1, int power1, core::Symbol_sp unit2, int power2)
    {_OF();
	int matchPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	double amountScale = 1.0;
	this->parseDimensionSymbols(matchPowers,amountScale,unit1,power1);
	this->parseDimensionSymbols(matchPowers,amountScale,unit2,power2);
	this->throwOnPowersMisMatch(matchPowers);
	Vector3 vc = this->_Coordinates->getElement(i).multiplyByScalar(1.0/amountScale);
	return vc;
    }




    Object_sp CoordinateArray_O::amount() const
    {_OF();
	int zeroPowers[NumBaseDimensions] = {0,0,0,0,0,0,0};
	if ( this->powersMatch(zeroPowers) )
	{
	    return this->_Coordinates->copy();
	}
	THROW(_lisp->error(BF("You cannot get the amount of a Quantity unless it is unit-less!")));
    }




    void CoordinateArray_O::append(Xyz_sp xyz)
    {_OF();
	ASSERTF(this->_System==xyz->_System,BF("Systems must match"));
	this->_Coordinates->append(xyz->_Pos);
    }




    Object_sp CoordinateArray_O::unsafeAmount() const
    {_OF();
	return this->_Coordinates->copy();
    }


    string CoordinateArray_O::__repr__() const
    {_OF();
	stringstream ss;
	double amount = 1.0;
	System_sp system = this->system();
	ss << this->_class()->classNameSymbol()->__repr__() << "*"<< this->system()->unitsAsString(amount,this->sharedThis<Quantity_O>());
	return ss.str();
    }


    CoordinateArray_sp CoordinateArray_O::copy() const
    {_OF();
	CoordinateArray_sp ca = RP_Copy(this);
	ca->_Coordinates = this->_Coordinates->copy();
	return ca;
    }



    string CoordinateArray_O::dump() const
    {_OF();
	stringstream ss;
	for ( int i=0; i<this->size(); i++ )
	{
	    ss << Xyz_O::create(this->_Coordinates->getElement(i),this->sharedThis<CoordinateArray_O>(),_lisp)->__repr__() << endl;
	}
	return ss.str();
    }



}; /* units */
