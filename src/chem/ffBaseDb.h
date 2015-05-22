       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffBase.h
 *
 *	Maintains a database of stretch types
 */

#ifndef FFBASE_H
#define	FFBASE_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "core/common.h"
#include "bond.h"
#include "addon/vector3.h"
#include "atom.h"
#include "residue.h"


#include "chemPackage.h"


namespace chem {

SMART(ForceField);

};

namespace chem {
    typedef	enum { parameterized, estimated, rough } ParameterizationLevel;

};


namespace kw {
    extern core::Symbol_sp _sym_parameterized;
    extern core::Symbol_sp _sym_estimated;
    extern core::Symbol_sp _sym_rough;
};

namespace translate {
    template <> struct to_object<chem::ParameterizationLevel>
    {								
	typedef	chem::ParameterizationLevel	GivenType;	
	static core::T_sp convert(const GivenType& val)
	{_G();
	    SYMBOL_EXPORT_SC_(ChemKwPkg,parameterized);
	    SYMBOL_EXPORT_SC_(ChemKwPkg,estimated);
	    SYMBOL_EXPORT_SC_(ChemKwPkg,rough);
	    switch (val) {
	    case chem::parameterized:
		return (kw::_sym_parameterized);
	    case chem::estimated:
		return (kw::_sym_estimated);
	    case chem::rough:
		return (kw::_sym_rough);
	    };
	    SIMPLE_ERROR(BF("Illegal ParameterizationLevel"));
	}
    };

    template <>
    struct from_object<chem::ParameterizationLevel>
    {								
	typedef	chem::ParameterizationLevel 	ExpectedType;
	typedef	chem::ParameterizationLevel 	DeclareType;
	DeclareType _v;
	from_object(gctools::smart_ptr<core::T_O> o) {
	    if ( o == kw::_sym_parameterized ) {
		_v = chem::parameterized;
	    } else if ( o == kw::_sym_estimated ) {
		_v = chem::estimated;
	    } else if ( o == kw::_sym_rough ) {
		_v = chem::rough;
	    } else {
		SIMPLE_ERROR(BF("Illegal value for ParameterizationLevel %s") % _rep_(o) );
	    }
	}
    };


};


namespace chem {




typedef	enum { stretch, angle, ptor, itor, nonbond, vdw } ParameterType;


class FFParameter_O;
    typedef	gctools::smart_ptr<FFParameter_O>	FFParameter_sp;
SMART(FFParameter);
class FFParameter_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,FFParameter_O,"FFParameter");

public:
    void	archiveBase(core::ArchiveP node);
public:
	ParameterizationLevel	_Level;



/*! Return the parameterization level of this parameter.
 * In most cases if it isn't "excellent" then it will be noted and the
 * user can go and calculate a good parameter and put it into the force field.
 */
	ParameterizationLevel	level() { return this->_Level;};

/*! Describe this parameter in terms of its ParameterizationLevel and the atom types involved.
 * This will be used when collecting missing or poorly parameterized parameters.
 */
virtual	string levelDescription();
    virtual	ParameterType type() {SUBIMP();};

DEFAULT_CTOR_DTOR(FFParameter_O);
};



SMART(FFBaseDb);
class FFBaseDb_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,FFBaseDb_O,"FFBaseDb");

public:
    void	archiveBase(core::ArchiveP node);

public:
	ForceField_wp	_ForceField;

public:

	void	setForceField(ForceField_sp ff);

	DEFAULT_CTOR_DTOR(FFBaseDb_O);
};

};

TRANSLATE(chem::FFParameter_O);
TRANSLATE(chem::FFBaseDb_O);
#endif
