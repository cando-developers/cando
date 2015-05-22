      
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	forceField.h
 *
 *	Maintains a parameter set for the Moe general force field
 *	from these parameters we can derive any force field parameters
 *
 */

#ifndef FORCEFIELD_H
#define	FORCEFIELD_H

#include "core/common.h"
#include "stdio.h"
#include "string"
#include "vector"
#include "set"
#include "bond.h"
#include "addon/vector3.h"
//#include	"conformation.h"
#include "atom.h"
#include "residue.h"
#include "ffBaseDb.h"
#include "ffTypesDb.h"
#include "ffStretchDb.h"
#include "ffAngleDb.h"
#include "ffItorDb.h"
#include "ffPtorDb.h"
#include "ffNonbondDb.h"
#include "ffVdwDb.h"



namespace       chem {


extern string	ForceField_XmlName;



typedef	enum	{ ffStr, ffAng, ffStb, ffOop, ffItor, ffTor, ffVdw, ffEle, ffSol } DisableEnum;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
class	InfoDb_O;
    typedef	gctools::smart_ptr<InfoDb_O>	InfoDb_sp;
SMART(InfoDb);
class InfoDb_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,InfoDb_O,"InfoDb");

public:
    void	archiveBase(core::ArchiveP node);

public:
    core::SymbolMap<core::Str_O>	_database;
public:

    void	addInfo( core::Symbol_sp key, const string& data );
	DEFAULT_CTOR_DTOR(InfoDb_O);
};



////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

SMART(ForceField);
class ForceField_O : public core::T_O
{
    LISP_BASE1(core::T_O);
    LISP_CLASS(chem,ChemPkg,ForceField_O,"ForceField");

public:
	void initialize();
public:
	void archiveBase(core::ArchiveP node);

public:
//		/*! Read the forceField from a file*/
//	static ForceField_sp	open_ForceField(const string& fn);

public:
	string				_Title;
	string				_Ref;
    gctools::Vec0<core::Symbol_sp>		_SingleBondMultiBondDistinctions;
	InfoDb_sp			_Info;
	FFTypesDb_sp			_Types;
	FFStretchDb_sp			_Stretches;
	FFAngleDb_sp			_Angles;
	FFItorDb_sp			_Itors;
	FFPtorDb_sp			_Ptors;
	FFNonbondDb_sp			_Nonbonds;
	FFVdwDb_sp			_Vdws;
//	bool				_Disable[sizeof(DisableEnum)];
#if 0
	RPFFEquivalenceDb		Equivalences;
	RPFFOutOfPlaneDb		OutOfPlanes;
	FFStretch_spBendDb		StretchBends;
#endif

    void	addSingleBondMultiBondDistinctionType(core::Symbol_sp s) {
			this->_SingleBondMultiBondDistinctions.push_back(s);
		};

	FFTypesDb_sp	getTypes() { return this->_Types; };

	FFStretchDb_sp getStretchDb() { return this->_Stretches;};
	FFAngleDb_sp getAngleDb() { return this->_Angles;};
	FFItorDb_sp getItorDb() { return this->_Itors;};
	FFPtorDb_sp getPtorDb() { return this->_Ptors;};
	FFNonbondDb_sp	getNonbondDb() { return this->_Nonbonds; };
	FFVdwDb_sp getVdwDb() { return this->_Vdws;};
	
	void	assignTypes(Matter_sp matter);

	void	setTitle(const string& title);
	void	setInfoDb( InfoDb_sp Info );
	void	setFFTypeDb( FFTypesDb_sp Types);
	void	setFFStretchDb( FFStretchDb_sp Stretches);
	void	setFFAngleDb( FFAngleDb_sp Angles);
	void	setFFItorDb( FFItorDb_sp Itors);
	void	setFFPtorDb( FFPtorDb_sp Ptors);
	void	setFFNonbondDb(FFNonbondDb_sp Nonbonds );
	void	setFFVdwDb(FFVdwDb_sp Vdws );

//	void		parseFromMoeStream(istream in);
	void		parseFromMoeFormatFileName(string name);


	void	pointAllPartsToForceField(ForceField_sp ff);

	void	saveAs(const string& fileName);

	DEFAULT_CTOR_DTOR(ForceField_O);
};



//
// Some forcefield values may be undefined
// they will be set to UNDEFINED_NUMBER
// test this with isUndefined(x)
//
#define	UNDEFINED_NUMBER	9.0e99
#define	UNDEFINED_CUTOFF	9.0e98
inline	bool	isUndefined(double x) {return (x)>=UNDEFINED_CUTOFF;}



};


TRANSLATE(chem::InfoDb_O);
TRANSLATE(chem::ForceField_O);
#endif
