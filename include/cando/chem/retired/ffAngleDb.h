       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffAngle.h
 *
 *	Maintains a database of angle terms
 */

#ifndef FFANGLE_H
#define	FFANGLE_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/candoBase/vector3.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>

#include <cando/chem/chemPackage.h>


namespace chem {

#define		FFAngleWildCardCode	-1

inline	string	XmlTag_FFAngle() { return "FFAngle";};

class	FFAngle_O;
    typedef	gctools::smart_ptr<FFAngle_O>	FFAngle_sp;
SMART(FFAngle);
class FFAngle_O : public FFParameter_O
{
    LISP_BASE1(FFParameter_O);
    LISP_CLASS(chem,ChemPkg,FFAngle_O,"FFAngle");

public:
	void initialize();
public:
	void	archiveBase(core::ArchiveP node);

public:
    core::Symbol_sp 		_Type1;
    core::Symbol_sp 		_Type2;
    core::Symbol_sp 		_Type3;
	double		_AngRad;
	double		_K2__kJPerRadianSquared;
//	double		_K3;
//	double		_K4;
//	double		_Ub_k;
//	double		_Ub_len;

	double getAngle_Radian() const { return this->_AngRad;};

	void setK2_kCalPerRadianSquared(double val);
	void setK2_kJPerRadianSquared(double val);
	
	double getK2_kCalPerRadianSquared() const;
	double getK2_kJPerRadianSquared() const;

virtual	string	levelDescription();
virtual ParameterType type() { return angle; };

DEFAULT_CTOR_DTOR(FFAngle_O);
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


inline	string	XmlTag_FFAngleDb() { return "FFAngleDb";};
class	FFAngleDb_O;
    typedef	gctools::smart_ptr<FFAngleDb_O>	FFAngleDb_sp;
SMART(FFAngleDb);
class FFAngleDb_O : public FFBaseDb_O
{
    LISP_BASE1(FFBaseDb_O);
    LISP_CLASS(chem,ChemPkg,FFAngleDb_O,"FFAngleDb");

public:
    void	archiveBase(core::ArchiveP node);
public:
    typedef	gctools::Vec0<FFAngle_sp>::iterator FFAngle_spIterator;
	string			_AngleFunction;
    gctools::Vec0<FFAngle_sp>	_Terms;
	core::SymbolMap<FFAngle_O>          _Lookup;
	core::SymbolMap<core::DoubleFloat_O>		_ZConstants;	//<! Angle parameter constants from Wang et al. J. Comput. Chem 25, 1157-1174 (2004)
	core::SymbolMap<core::DoubleFloat_O>		_CConstants;	//<! Angle parameter constants from Wang et al. J. Comput. Chem 25, 1157-1174 (2004)

public:
	void	setAngleFunction(const string& st) {this->_AngleFunction = st;};
	void	add( FFAngle_sp ang );

    void	addZConstant(core::Symbol_sp element, double value ) { this->_ZConstants.set(element,core::DoubleFloat_O::create(value));};
    void	addCConstant(core::Symbol_sp element, double value ) { this->_CConstants.set(element,core::DoubleFloat_O::create(value));};

	FFAngle_spIterator	begin() {return this->_Terms.begin();};
	FFAngle_spIterator	end() {return this->_Terms.end();};


	FFAngle_sp	findTerm( chem::Atom_sp a1, chem::Atom_sp a2, chem::Atom_sp a3);
	FFAngle_sp	estimateTerm( chem::Atom_sp a1, chem::Atom_sp a2, chem::Atom_sp a3 );

	void	initialize();

	DEFAULT_CTOR_DTOR(FFAngleDb_O);
};

};

TRANSLATE(chem::FFAngle_O);
TRANSLATE(chem::FFAngleDb_O);
#endif
