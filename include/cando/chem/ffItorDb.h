       
       
//
// (C) 2004 Christian E. Schafmeister
//


/*
 *	ffItor.h
 *
 *	Maintains a database of itor terms
 */

#ifndef FFITOR_H
#define	FFITOR_H
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include <clasp/core/common.h>
#include <cando/chem/bond.h>
#include <cando/geom/vector3.h>
#include <cando/adapt/symbolMap.h>
#include <cando/chem/atom.h>
#include <cando/chem/residue.h>
#include <cando/chem/ffBaseDb.h>
#include <clasp/core/smallMap.h>

namespace       chem
{

#define	FFItorWildCardCode	-1


    inline	string	XmlTag_FFItor() { return "FFItor";};

    SMART(FFItor);
    class FFItor_O : public FFParameter_O
    {
	LISP_BASE1(FFParameter_O);
	LISP_CLASS(chem,ChemPkg,FFItor_O,"FFItor");
    public:
	static int const IMaxPeriodicity = 6;
    public:
	void initialize();
#ifdef XML_ARCHIVE
	void	archiveBase(core::ArchiveP node);
#endif

    public:
	core::Symbol_sp         _T1;
	core::Symbol_sp 	_T2;
	core::Symbol_sp 	_T3;
	core::Symbol_sp 	_T4;
	bool      	_hasPeriodicity[IMaxPeriodicity];
	double      _Vs_kJ[IMaxPeriodicity];
	double      _PhaseRads[IMaxPeriodicity];

	void	setTypes( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4);

	bool hasPeriodicity(int period) const;
	double  getV_kJ(int period) const;
	double  getV_kCal(int period) const;
	void    setV_kJ(int period, double v);
	void    setV_kCal(int period, double v);
	double  getPhaseRad(int period) const;
	void    setPhaseRad(int period, double v);

	void    mergeWith(FFItor_sp itor);


	virtual	string	levelDescription();
	virtual ParameterType type() { return itor; };
	DEFAULT_CTOR_DTOR(FFItor_O);
    };




////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
    SMART(FFItorDb);
    class FFItorDb_O : public FFBaseDb_O
    {
	LISP_BASE1(FFBaseDb_O);
	LISP_CLASS(chem,ChemPkg,FFItorDb_O,"FFItorDb");

    public:
#ifdef XML_ARCHIVE
	void	archiveBase(core::ArchiveP node);
#endif

    public:
	typedef	gctools::Vec0<FFItor_sp>::iterator	FFItor_spIterator;
        gctools::Vec0<FFItor_sp>		_Terms;
	adapt::SymbolMap<FFItor_O>		_Lookup;

	FFItor_spIterator	begin() { return this->_Terms.begin();};
	FFItor_spIterator	end() { return this->_Terms.end();};

	void	add( FFItor_sp itor );

	FFItor_sp findExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
	bool    hasExactTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

	//! Look for exact term then general one
        FFItor_sp findBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );
        bool    hasBestTerm( core::Symbol_sp a1, core::Symbol_sp a2, core::Symbol_sp a3, core::Symbol_sp a4 );

        void    cantFind(core::Symbol_sp t1, core::Symbol_sp t2, core::Symbol_sp t3, core::Symbol_sp t4 );

	//! Dump all ptors to stdout that match these types
	DEFAULT_CTOR_DTOR(FFItorDb_O);
    };


};

TRANSLATE(chem::FFItor_O);
TRANSLATE(chem::FFItorDb_O);
#endif
